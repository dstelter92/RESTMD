/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Notes:
    The stmd.f::stmdcntrl() subroutine is translated in the FixSTMD constructor and init().
    The stmd.f::stmdinitrans() subroutine is translated in init().
    The stmd.f::stmdinitoutu() subroutine can be added to the bottom of init() if desired.
    If stmd.f::stmdAssignVel() is just velocity-rescaling, then this is not needed and a LAMMPS thermostat fix can be used.
    The stmd.f::stmddig() subroutine is translated to dig().
    ...
    The stmd.f::stmdMain() subroutine is translated to Main().

    The Verlet::run() function shows the flow of calculation at each
    MD step. The post_force() function is the current location for the
    STMD update via Main() and subsequent scaling of forces.

    In several spots, ".eq." was used in the Fortran code for testing reals. That is 
    repeated here with "==", but probably should be testing similarity against some tolerance.

    Many of the smaller subroutines/functions could probably be inlined.
    
------------------------------------------------------------------------- */

#include <cstring>
#include <cstdlib>
#include "fix_stmd.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "comm.h"
#include "group.h"
#include "compute.h"
#include "output.h"
#include "universe.h"
#include <fstream>


using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define INVOKED_SCALAR 1

/* ---------------------------------------------------------------------- */

FixStmd::FixStmd(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 15 || narg > 16) error->all(FLERR,"Illegal fix stmd command");

  // DEBUG FLAG
  //stmd_debug = 1;

  scalar_flag = 1;
  vector_flag = 1;
  array_flag = 1;
  
  extscalar = 0;
  extvector = 0;
  extarray = 0;
  global_freq = 1;
  restart_file = 1;
 
  // This is the subset of variables explicitly given in the charmm.inp file
  // If the full set is expected to be modified by a user, then reading 
  // a stmd.inp file is probably the best mechanism for input.
  //
  // fix fxstmd all stmd RSTFRQ f_style init_f final_f Tlo Thi Elo Ehi binsize TSC1 TSC2 ST OREST

  // Probably a good idea to set this equal to restart value in input
  RSTFRQ = atoi(arg[3]);      

  // Setup type of f-reduction scheme
  f_flag = -1;
  if (strcmp(arg[4],"none") == 0)
    f_flag = 0;
  else if (strcmp(arg[4],"hchk") == 0) 
    f_flag = 1;
  else if (strcmp(arg[4],"sqrt") == 0)
    f_flag = 2;
  else if (strcmp(arg[4],"constant_f") == 0)
    f_flag = 3;
  else if (strcmp(arg[4],"constant_df") == 0)
    f_flag = 4;
  else
    error->all(FLERR,"STMD: invalid f-reduction scheme");
  if (f_flag == -1)
    error->all(FLERR,"STMD: invalid f-reduction scheme");
  
  // Only used initially, controlled by restart
  initf  = atof(arg[5]); // Delta F, init
  if (initf > 1.) 
    error->all(FLERR,"STMD: initial deltaF value too large");
  
  // Delta-f tolerances
  dFval3 = 0.000020;
  if (dFval3 < 0.00001)
    error->all(FLERR,"STMD: final deltaF value too small");
  dFval4 = dFval3 / 10.;
  pfinFval = exp(dFval3 * 2 * bin);
  finFval = exp(dFval4 * 2 * bin);

  // Used once, controlled by restart file
  TL     = atof(arg[6]);
  TH     = atof(arg[7]);

  // Controlled by input file, not restart file
  Emin   = atof(arg[8]);
  Emax   = atof(arg[9]);
  bin    = atof(arg[10]);
  TSC1   = atoi(arg[11]);
  TSC2   = atoi(arg[12]);

  // find current thermostating fix
  int n = strlen(arg[13])+1;
  id_nh = new char[n];
  strcpy(id_nh,arg[13]);
  
  // check if fix exists
  int ifix = modify->find_fix(id_nh);
  if (ifix < 0)
    error->all(FLERR,"Fix id for nvt or npt fix does not exist");
  Fix *nh = modify->fix[ifix];

  // 0 for new run, 1 for restart
  OREST = -1;
  if (strcmp(arg[14],"yes") == 0)
    OREST = 1;
  else if (strcmp(arg[14],"no") == 0)
    OREST = 0;
  else 
    error->all(FLERR,"STMD: invalid restart option");
  if (OREST == -1)
    error->all(FLERR,"STMD: invalid restart option");

  // Make dir_output default to local dir
  if (narg == 16)
    strcpy(dir_output,arg[15]);
  else
    strcpy(dir_output,"./");
  
  // Init arrays
  Y2 = Prob = NULL;
  Hist = Htot = PROH = NULL;

  // STMD_specific flags
  hist_flag = 0; // 0=read from restart, 1=reset
  freset_flag = 0; // 0=read from restart, 1=reset

  // Setup communication flags
  stmd_logfile = stmd_debug = stmd_screen = 0;
  if ((comm->me == 0) && (logfile)) stmd_logfile = 1;
  if ((comm->me == 0) && (screen)) stmd_screen = 1;

  // Init file pointers
  fp_wtnm = fp_whnm = fp_whpnm = fp_orest = NULL;

  
  // Energy bin setup
  BinMin = round(Emin / bin);
  BinMax = round(Emax / bin);
  N = BinMax - BinMin + 1;

  if (N < 0)
    error->all(FLERR,"Emin > Emax, negative number of bins");
  if (N < 1)
    error->all(FLERR,"Invalid energy range");


  // ceate new compute temp style
  // id = fix-ID + temp
  n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;

  // do same for pressure/stmd
  // id = fix-ID + press
  n = strlen(id) + 6;
  id_press = new char[n];
  strcpy(id_press,id);
  strcat(id_press,"_press");

  newarg = new char*[5];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "PRESSURE/STMD";
  newarg[3] = id_temp;
  newarg[4] = id;
  modify->add_compute(5,newarg);
  delete [] newarg;

  // check if barostat exists
  pressflag = 0;
  int *p_flag = (int *)nh->extract("p_flag",ifix);
  if ((p_flag == NULL) || (ifix != 1) || (p_flag[0] == 0)
      || (p_flag[1] == 0) || (p_flag[2] == 0)) {
    pressflag = 0;
  } else if ((p_flag[0] == 1) && (p_flag[1] == 1) && (p_flag[2] == 1) && (ifix == 1)) {
    pressflag = 1;
    char *modargs[2];
    modargs[0] = (char*) "press";
    modargs[1] = id_press;
    nh->modify_param(2,modargs);
  }
  
  // Setup size of global vector/arrays
  size_vector = 9;
  size_array_cols = 4;
  size_array_rows = N;

}

/* ---------------------------------------------------------------------- */

FixStmd::~FixStmd()
{
  memory->destroy(Y2);
  memory->destroy(Hist);
  memory->destroy(Htot);
  memory->destroy(PROH);
  memory->destroy(Prob);
  modify->delete_compute(id_temp);
  modify->delete_compute(id_press);
  delete [] id_nh;
  delete [] id_temp;
  delete [] id_press;
}

/* ---------------------------------------------------------------------- */

int FixStmd::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStmd::init()
{
  // Get number of replicas (worlds) and walker number
  nworlds = universe->nworlds;
  iworld = universe->iworld;
  char walker[256]; 
  sprintf(walker,"%i",iworld);

  if (comm->me == 0) {
    char filename[256];
    if (!fp_wtnm) {
      strcpy(filename,dir_output);
      strcat(filename,"/WT.");
      strcat(filename,walker);
      strcat(filename,".d");
      strcpy(filename_wtnm,filename);
      fp_wtnm  = fopen(filename,"w");
    }
    if (!fp_whnm) {
      strcpy(filename,dir_output);
      strcat(filename,"/WH.");
      strcat(filename,walker);
      strcat(filename,".d");
      strcpy(filename_whnm,filename);
      fp_whnm  = fopen(filename,"w");
    }
    /*
    if (!fp_whpnm) {
      strcpy(filename,dir_output);
      strcat(filename,"/WHP.");
      strcat(filename,walker);
      strcat(filename,".d");
      strcpy(filename_whpnm,filename);
      fp_whpnm = fopen(filename,"w");
    }
    */
    if ((!fp_orest) && (!OREST)) {
      strcpy(filename,dir_output);
      strcat(filename,"/oREST.");
      strcat(filename,walker);
      strcat(filename,".d");
      strcpy(filename_orest,filename);
      fp_orest = fopen(filename,"w");
    }
    if ((!fp_orest) && (OREST)) {
      strcpy(filename,dir_output);
      strcat(filename,"/oREST.");
      strcat(filename,walker);
      strcat(filename,".d");
      strcpy(filename_orest,filename);

      // Check if file exists
      fp_orest = fopen(filename,"r");
      if (!fp_orest) {
        if (stmd_logfile)
          fprintf(logfile,"Restart file: %soREST.%s.d is empty\n",dir_output,walker);
        if (stmd_screen)
          fprintf(screen,"Restart file: %soREST.%s.d is empty\n",dir_output,walker);
        error->one(FLERR,"STMD: Restart file does not exist\n");
      }
    }
  }

  // set compute ptrs
  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature compute ID for fix grem does not exist");
  temperature = modify->compute[icompute];

  // check if thermostat fix exists
  int ifix = modify->find_fix(id_nh);
  if (ifix < 0)
    error->all(FLERR,"Fix id for nvt or npt fix does not exist");
  Fix *nh = modify->fix[ifix];

  // check for temperature ramp
  double *t_start = (double *)nh->extract("t_start",ifix);
  double *t_stop = (double *)nh->extract("t_stop",ifix);
  if ((t_start != NULL) && (t_stop != NULL) && (ifix == 0)) {
    ST = *t_start;
    if (*t_start != *t_stop)
      error->all(FLERR,"Thermostat temperature ramp not allowed");
  } else
    error->all(FLERR,"Problem extracting target temperature from fix nvt or npt");

  // check for supported barostat settings
  pressref = 0.0;
  if (pressflag) {
    int *p_flag = (int *)nh->extract("p_flag",ifix);
    double *p_start = (double *) nh->extract("p_start",ifix);
    double *p_stop = (double *) nh->extract("p_stop",ifix);
    if ((p_flag != NULL) && (p_start != NULL) && (p_stop != NULL)
        && (ifix == 1)) {
      ifix = 0;
      pressref = p_start[0];
      if ((p_start[0] != p_stop[0]) || (p_flag[0] != 1)) ++ ifix;
      if ((p_start[1] != p_stop[1]) || (p_flag[0] != 1)) ++ ifix;
      if ((p_start[2] != p_stop[2]) || (p_flag[0] != 1)) ++ ifix;
      if ((p_start[0] != p_start[1]) || (p_start[1] != p_start[2])) ++ifix;
      if ((p_flag[3] != 0) || (p_flag[4] != 0) || (p_flag[5] != 0)) ++ifix;
      if (ifix > 0)
        error->all(FLERR,"Unsupported pressure settings in fix npt");
    } else
      error->all(FLERR,"Problem extracting target pressure from fix npt");
  }
  
  // Defaults
  CutTmin  = 50.0;
  CutTmax  = 50.0;
  HCKtol   = 0.2;

  STG     = 1;
  SWf     = 1;
  SWfold  = 1;
  Gamma   = 1.0;
  Count   = 0;
  CountH  = 0;
  totC    = 0;
  totCi   = 0;
  SWchk   = 1;
  CountPH = 0;
  T = ST; // latest T_s at bin i

  /*
  // Delta-f tolerances
  dFval3 = 0.0001;
  dFval4 = 0.000002;
  */
  pfinFval = exp(dFval3 * 2 * bin);
  finFval = exp(dFval4 * 2 * bin);

  f = exp(initf * 2 * bin);
  df = log(f) * 0.5 / bin;

  T0 = ST;
  T1 = TL / ST;
  T2 = TH / ST;
  CTmin = (TL + CutTmin) / ST;
  CTmax = (TH - CutTmax) / ST;

  memory->grow(Y2, N, "FixSTMD:Y2");
  memory->grow(Hist, N, "FixSTMD:Hist");
  memory->grow(Htot, N, "FixSTMD:Htot");
  memory->grow(PROH, N, "FixSTMD:PROH");
  memory->grow(Prob, N, "FixSTMD:Prob");

  for (int i=0; i<N; i++) {
    Y2[i] = T2;
    Hist[i] = 0;
    Htot[i] = 0;
    PROH[i] = 0;
    Prob[i] = 0.0;
  }

  // Search for pe compute, otherwise create a new one
  pe_compute_id = -1;
  for (int i=0; i<modify->ncompute; i++) {
    if (strcmp(modify->compute[i]->style,"pe") == 0) {
      pe_compute_id = i;
      break;
    }
  }

  // Pretty sure a pe compute is always present, but check anyways.
  // Did we find a pe compute? If not, then create one.
  if (pe_compute_id < 0) {
    int n = strlen(id) + 4;
    id_pe = new char[n];
    strcpy(id_pe,id);
    strcat(id_pe,"_pe");

    char **newarg = new char*[3];
    newarg[0] = id_pe;
    newarg[1] = group->names[igroup];
    newarg[2] = (char *) "pe";

    modify->add_compute(3,newarg);
    delete [] newarg;

    pe_compute_id = modify->ncompute - 1;
  }

  if (domain->triclinic)
    error->all(FLERR,"Triclinic cells are not supported");

  if (OREST) { // Read oREST.d into variables
    if (comm->me == 0) {
      int k = 0;
      int numb = 13;
      int nsize = 3*N + numb;
      double *list;
      memory->create(list,nsize,"stmd:list");

      char filename[256];
      strcpy(filename,dir_output);
      strcat(filename,"/oREST.");
      strcat(filename,walker);
      strcat(filename,".d");
      std::ifstream file(filename);

      // Check if file empty before getting data
      file.seekg(0,file.end);
      int sz = file.tellg();
      file.seekg(0,file.beg);
      if (sz < (nsize*sizeof(double))) {
        if (stmd_logfile)
          fprintf(logfile,"Restart file: %soREST.%s.d is an invalid format\n",dir_output,walker);
        if (stmd_screen)
          fprintf(screen,"Restart file: %soREST.%s.d is an invalid format\n",dir_output,walker);
        error->one(FLERR,"STMD: Restart file is empty/invalid\n");
      }

      for (int i=0; i<nsize; i++) 
        file >> list[i];

      STG = static_cast<int> (list[k++]);
      if (!freset_flag)
        f = list[k++];
      CountH = static_cast<int> (list[k++]);
      SWf = static_cast<int> (list[k++]);
      SWfold = static_cast<int> (list[k++]);
      SWchk = static_cast<int> (list[k++]);
      Count = static_cast<int> (list[k++]);
      totCi = static_cast<int> (list[k++]);
      CountPH = static_cast<int> (list[k++]);
      //TSC1 = static_cast<int> (list[k++]);
      //TSC2 = static_cast<int> (list[k++]);
      T1 = list[k++];
      T2 = list[k++];
      CTmin = list[k++];
      CTmax = list[k++];

      for (int i=0; i<N; i++)
        Y2[i] = list[k++];
      for (int i=0; i<N; i++)
        Htot[i] = list[k++];
      if (!hist_flag) {
        for (int i=0; i<N; i++)
          PROH[i] = list[k++];
      }

      memory->destroy(list);
    }
    if (!freset_flag)
      df = log(f) * 0.5 / bin;
    OREST = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixStmd::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    post_force(vflag);

  // Write info to screen/log
  if ((stmd_logfile) && (nworlds > 1))
    fprintf(logfile,"RESTMD: #replicas=%i  walker=%i\n",nworlds,iworld);
  if ((stmd_screen) && (nworlds > 1))
    fprintf(screen,"RESTMD: #replicas=%i  walker=%i\n",nworlds,iworld);
  if (stmd_logfile) {
    fprintf(logfile,"STMD: STAGE=%i, #bins=%i  binsize=%f\n",STG,N,bin); 
    fprintf(logfile,"  Emin=%f Emax=%f f-value=%f df=%f\n",Emin,Emax,f,df); 
    fprintf(logfile,"  f-tolerances: STG3=%f STG4=%f\n",pfinFval,finFval);
  }
  if (stmd_screen) {
    fprintf(screen,"STMD: STAGE=%i, #bins=%i  binsize=%f\n",STG,N,bin);
    fprintf(screen,"  Emin=%f Emax=%f f-value=%f df=%f\n",Emin,Emax,f,df); 
    fprintf(screen,"  f-tolerances: STG3=%f STG4=%f\n",pfinFval,finFval);
  }

  // Write current Ts estimate to logfile
  if ((stmd_logfile) && (stmd_debug)) {
    fprintf(logfile,"STMD Temperature (Y2)= ");
    for (int i=0; i<N; i++) 
      fprintf(logfile," %f",Y2[i]);
    fprintf(logfile,"\n");
  }
    // Force computation of energies
    modify->compute[pe_compute_id]->invoked_flag |= INVOKED_SCALAR;
    modify->addstep_compute(update->ntimestep + 1);
  } else
    error->all(FLERR,"Currently expecting run_style verlet");
}

/* ---------------------------------------------------------------------- */

void FixStmd::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixStmd::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // Get current value of potential energy from compute/pe
  double tmp_pe = modify->compute[pe_compute_id]->compute_scalar();
  double tmp_vol = domain->xprd * domain->yprd * domain->zprd;

  sampledE = tmp_pe + (pressref*tmp_vol/(force->nktv2p));

  // Check if sampledE is outside of bounds before continuing
  if ((sampledE < Emin) || (sampledE > Emax))
    error->all(FLERR,"STMD: Sampled energy out of range\n");

  // Master rank will compute scaling factor and then Bcast to world
  MAIN(update->ntimestep,sampledE);

  // Gamma(U) = T_0 / T(U)
  MPI_Bcast(&Gamma, 1, MPI_DOUBLE, 0, world); 

  // Scale forces
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      f[i][0]*= Gamma;
      f[i][1]*= Gamma;
      f[i][2]*= Gamma;
    }
}

/* ---------------------------------------------------------------------- */

void FixStmd::end_of_step()
{
  // Force computation of energies on next step
  modify->compute[pe_compute_id]->invoked_flag |= INVOKED_SCALAR;
  modify->addstep_compute(update->ntimestep + 1);

  // If stmd, write output, otherwise let temper/stmd handle it
  if (universe->nworlds == 1) {
    write_temperature();
    write_orest();
  }
}

/* ---------------------------------------------------------------------- */

void FixStmd::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixStmd::memory_usage()
{
  double bytes = 0.0;
  bytes+= 7 * N * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   write temperature to external file
------------------------------------------------------------------------- */

void FixStmd::write_temperature()
{
  int istep = update->ntimestep;
  int m = istep % RSTFRQ;
  if ((m == 0) && (comm->me == 0)) {
    fprintf(fp_wtnm,"### STMD Step %i: bin E Ts(E)\n",istep);
    for (int i=0; i<N; i++) 
      fprintf(fp_wtnm,"%i %f %f\n", i,(i*bin)+Emin,Y2[i]*ST);
    fprintf(fp_wtnm,"\n\n");
  }
}

/* ----------------------------------------------------------------------
   write external restart file
------------------------------------------------------------------------- */

void FixStmd::write_orest()
{
  // Write restart info to external file
  int m = (update->ntimestep) % RSTFRQ;
  if ((m == 0) && (comm->me == 0)) {
    int k = 0;
    int numb = 13;
    int nsize = 3*N + numb;
    double *list;
    memory->create(list,nsize,"stmd:list");

    list[k++] = STG;
    list[k++] = f;
    list[k++] = CountH;
    list[k++] = SWf;
    list[k++] = SWfold;
    list[k++] = SWchk;
    list[k++] = Count;
    list[k++] = totCi;
    list[k++] = CountPH;

    // Control these by LAMMPS input file
    //list[k++] = TSC1; 
    //list[k++] = TSC2;

    list[k++] = T1;
    list[k++] = T2;

    // only used for chk_hist
    list[k++] = CTmin;
    list[k++] = CTmax;

    for (int i=0; i<N; i++) 
      list[k++] = Y2[i];
    for (int i=0; i<N; i++) 
      list[k++] = Htot[i];
    for (int i=0; i<N; i++) 
      list[k++] = PROH[i];

    // wipe file contents...
    char filename[256];
    char walker[256];
    sprintf(walker,"%i",universe->iworld);
    strcpy(filename,dir_output);
    strcat(filename,"/oREST.");
    strcat(filename,walker);
    strcat(filename,".d");
    freopen(filename,"w",fp_orest);

    if (fp_orest == NULL) 
      error->all(FLERR,"Cannot open STMD restart file");

    if (stmd_debug && stmd_logfile) {
      fprintf(screen,"%d\n",STG);
      fprintf(screen,"%f\n",f);
      fprintf(logfile,"%d\n",STG);
      fprintf(logfile,"%f\n",f);
    }
    fprintf(fp_orest,"%d\n",STG);
    fprintf(fp_orest,"%f\n",f);
    for (int i=2; i<numb; i++) {
      if (stmd_debug && stmd_logfile) {
        fprintf(screen,"%f\n",list[i]);
        fprintf(logfile,"%f\n",list[i]);
      }
      fprintf(fp_orest,"%f\n",list[i]);
    }
    for (int i=numb; i<N+numb; i++) {
      if (stmd_debug && stmd_logfile) {
        fprintf(screen,"%f ",list[i]);
        fprintf(logfile,"%f ",list[i]);
      }
      fprintf(fp_orest,"%f ",list[i]);
    }
    if (stmd_debug && stmd_logfile) {
      fprintf(screen,"\n");
      fprintf(logfile,"\n");
    }
    fprintf(fp_orest,"\n");
    for (int i=N+numb; i<(2*N)+numb; i++) {
      if (stmd_debug && stmd_logfile) {
        fprintf(screen,"%f ",list[i]);
        fprintf(logfile,"%f ",list[i]);
      }
      fprintf(fp_orest,"%f ",list[i]);
    }
    if (stmd_debug && stmd_logfile) {
      fprintf(screen,"\n");
      fprintf(logfile,"\n");
    }
    fprintf(fp_orest,"\n");
    for (int i=(2*N)+numb; i<nsize; i++) {
      if (stmd_debug && stmd_logfile) {
        fprintf(screen,"%f ",list[i]);
        fprintf(logfile,"%f ",list[i]);
      }
      fprintf(fp_orest,"%f ",list[i]);
    }
    if (stmd_debug && stmd_logfile) {
      fprintf(screen,"\n");
      fprintf(logfile,"\n");
    }
    fprintf(fp_orest,"\n");

    memory->destroy(list);
  }
}

/* ----------------------------------------------------------------------
   Translation of stmd.f subroutines
------------------------------------------------------------------------- */

void FixStmd::dig()
{
  int nkeepmin = 0;
  double keepmin = Y2[nkeepmin];

  for (int i=0; i<N; i++) {
    if (Y2[i] <= keepmin) {
      keepmin = Y2[i];
      nkeepmin = i;
    }
  }

  for (int i=0; i<nkeepmin; i++) 
    Y2[i] = keepmin;
}

/* ---------------------------------------------------------------------- */

int FixStmd::Yval(double sampledE)
{
  curbin = round(sampledE / double(bin)) - BinMin + 1;
  int i = curbin;

  if ((i<1) || (i>N-1)) {
    if ((stmd_logfile) && (comm->me == 0))
      fprintf(logfile,"Error in Yval: pe=%f  bin=%f  i=%i\n",sampledE,bin,i);
    if ((stmd_screen) && (comm->me == 0))
      fprintf(screen,"Error in Yval: pe=%f  bin=%f  i=%i\n",sampledE,bin,i);
    error->all(FLERR,"STMD: Histogram index out of range");
  }

  double Yhi = Y2[i+1];
  double Ylo = Y2[i-1];

  Y2[i+1] = Y2[i+1] / (1.0 - df * Y2[i+1]);
  Y2[i-1] = Y2[i-1] / (1.0 + df * Y2[i-1]);

  if (stmd_debug && stmd_logfile) {
    fprintf(screen,"  STMD T-UPDATE: sampledE= %f  sampledbin= %i  df=%f\n",sampledE,i,df);
    fprintf(logfile,"  STMD T-UPDATE: sampledE= %f  sampledbin= %i  df=%f\n",sampledE,i,df);
    fprintf(screen,"    bin %d+1: T'= %f  T=%f  delta= %f\n",i,Y2[i+1],Yhi,Y2[i+1]-Yhi);
    fprintf(logfile,"    bin %d+1: T'= %f  T=%f  delta= %f\n",i,Y2[i+1],Yhi,Y2[i+1]-Yhi);
    fprintf(screen,"    bin %d-1: T'=%f  T=%f  delta= %f\n",i,Y2[i-1],Ylo,Y2[i-1]-Ylo);
    fprintf(logfile,"    bin %d-1: T'=%f  T=%f  delta= %f\n",i,Y2[i-1],Ylo,Y2[i-1]-Ylo);
  }

  if (Y2[i-1] < T1) 
    Y2[i-1] = T1;
  if (Y2[i+1] > T2) 
    Y2[i+1] = T2;

  return i;
}

/* ---------------------------------------------------------------------- */

void FixStmd::GammaE(double sampledE, int indx)
{
  const int i  = indx;
  const int im = indx - 1;
  const int ip = indx + 1;

  const double e = sampledE - double( round(sampledE / double(bin)) * bin );

  //double T;  // made this public to share with RESTMD
  if (e > 0.0) {
    const double lam = (Y2[ip] - Y2[i]) / double(bin);
    T = Y2[i] + lam * e;
  } else if (e < 0.0) {
    const double lam = (Y2[i] - Y2[im]) / double(bin);
    T = Y2[i] + lam * e;
  } else T = Y2[i];

  Gamma = 1.0 / T;
}

/* ---------------------------------------------------------------------- */

void FixStmd::AddedEHis(int i)
{
  Hist[i] = Hist[i] + 1;
  Htot[i] = Htot[i] + 1;
}

/* ---------------------------------------------------------------------- */

void FixStmd::EPROB(int icycle)
{
  int sw, m;
  const int indx = icycle;
  m = indx % TSC1;
  if ((m == 0) && (indx != 0)) sw = 1;

  m = indx % TSC2;
  if ((m == 0) && (indx != 0)) sw = 2;

  if (sw == 1) 
    for (int i=0; i<N; i++) 
      Prob[i] = Prob[i] / double(TSC1);
  else if (sw == 2)
    for (int i=0; i<N; i++) 
      Prob[i] = Prob[i] / double(TSC2);
  sw = 0;
}

/* ---------------------------------------------------------------------- */

void FixStmd::ResetPH()
{
  for (int i=0; i<N; i++) Hist[i] = 0;
}

/* ---------------------------------------------------------------------- */

void FixStmd::TCHK()
{
  if ((stmd_logfile) && (stmd_debug)) {
    fprintf(logfile,"  STMD TCHK: T1= %f (%f K)  Y2[0]= %f (%f K)\n",T1,T1*ST,Y2[0],Y2[0]*ST);
    fprintf(screen,"  STMD TCHK: T1= %f (%f K)  Y2[0]= %f (%f K)\n",T1,T1*ST,Y2[0],Y2[0]*ST);
  }
  if (Y2[0] == T1) STG = 2;
}

/* ---------------------------------------------------------------------- */

void FixStmd::HCHK()
{
  SWfold = SWf;

  int ichk = 0;
  int icnt = 0;
  double aveH = 0.0;

  // check CTmin and CTmax
  // average histogram
  for (int i=0; i<N; i++) {
    if ((Y2[i] > CTmin) && (Y2[i] < CTmax)) {
      aveH = aveH + double(Hist[i]);
      icnt++;
    }
  }

  if ((stmd_logfile) && (stmd_debug)) {
    fprintf(logfile,"  STMD CHK HIST: icnt= %i  aveH= %f  N= %i\n",icnt,aveH,N);
    fprintf(screen,"  STMD CHK HIST: icnt= %i  aveH= %f  N= %i\n",icnt,aveH,N);
  }
  if (icnt==0) return;

  aveH = aveH / double(icnt);

  double eval;
  for (int i=0; i<N; i++) {
    if ((Y2[i] > CTmin) && (Y2[i] < CTmax) ) {
      eval = abs(double(Hist[i] - aveH) / aveH);
      if (eval > HCKtol) ichk++;
      if ((stmd_logfile) && (stmd_debug)) {
        fprintf(logfile,"  STMD CHK HIST: totCi= %i  i= %i  eval= %f  HCKtol= %f  " 
            "ichk= %i  Hist[i]= %i\n",totCi,i,eval,HCKtol,ichk,Hist[i]);
        fprintf(screen,"  STMD CHK HIST: totCi= %i  i= %i  eval= %f  HCKtol= %f  " 
            "ichk= %i  Hist[i]= %i\n",totCi,i,eval,HCKtol,ichk,Hist[i]);
      }
    }
  }

  if (ichk < 1) SWf = SWf + 1;
}

/* ---------------------------------------------------------------------- */

void FixStmd::MAIN(int istep, double sampledE)
{
  Count = istep;
  totCi++;

  if (STG >= 3) CountPH++;

  if ((stmd_logfile) && (stmd_debug)) {
    fprintf(logfile,"STMD DEBUG: STAGE %i\n",STG);
    fprintf(screen,"STMD DEBUG: STAGE %i\n",STG);
    fprintf(logfile,"  STMD: Count=%i, f=%f\n",Count,f);
    fprintf(screen,"  STMD: Count=%i, f=%f\n",Count,f);
  }

  // Statistical Temperature Update
  int stmdi = Yval(sampledE);

  // Gamma Update
  GammaE(sampledE,stmdi);

  if ((stmd_logfile) && (stmd_debug)) {
    fprintf(logfile,"  STMD: totCi= %i Gamma= %f Hist[%i]= %i "
        "T= %f\n",totCi,Gamma,stmdi,Hist[stmdi],T);
    fprintf(screen,"  STMD: totCi= %i Gamma= %f Hist[%i]= %i "
        "T= %f\n",totCi,Gamma,stmdi,Hist[stmdi],T);
  }

  // Histogram Update
  AddedEHis(stmdi);
  CountH++;

  // Add to Histogram for production run
  if (STG >= 3) {
    PROH[stmdi]++;
    CountPH++;
  }

  // Hist Output
  int o = istep % RSTFRQ;
  if ((o == 0) && (comm->me == 0)) {
    fprintf(fp_whnm,"### STMD Step=%d: bin E hist thist phist\n",istep);
    for (int i=0; i<N; i++) 
      fprintf(fp_whnm,"%i %f %i %i %i\n",i,(i*bin)+Emin,Hist[i],Htot[i],PROH[i]);
    fprintf(fp_whnm,"\n\n");
  }
  
  // Production Run if STG >= 3
  // STG3 START: Check histogram and further reduce f until cutoff
  if (STG >= 3) {
    int m = istep % TSC2;

    // STMD, reduce f based on histogram flatness (original form)
    if (m == 0) {
      if ((stmd_logfile) && (stmd_debug)) {
        fprintf(logfile,"  STMD: istep= %i  TSC2= %i\n",istep,TSC2);
        fprintf(screen,"  STMD: istep= %i  TSC2= %i\n",istep,TSC2);
      }

      // Reduction based on histogram flatness
      if (f_flag == 1) { 
        HCHK(); // Check flatness
        if ((stmd_logfile) && (stmd_debug)) {
          fprintf(logfile,"  STMD: SWfold= %i  SWf= %i\n",SWfold,SWf);
          fprintf(logfile,"  STMD: f= %f  SWchk= %i\n",f,SWchk);
          fprintf(screen,"  STMD: SWfold= %i  SWf= %i\n",SWfold,SWf);
          fprintf(screen,"  STMD: f= %f  SWchk= %i\n",f,SWchk);
        }
        if (SWfold != SWf) {
          if (STG == 3) // dont reduce if STG4
            f = sqrt(f); // reduce f
          df = log(f) * 0.5 / bin;
          if ((stmd_logfile) && (stmd_debug)) {
            fprintf(logfile,"  STMD f-UPDATE: f= %f  SWf= %i  df= %f\n",f,SWf,df);
            fprintf(screen,"  STMD f-UPDATE: f= %f  SWf= %i  df= %f\n",f,SWf,df);
          }
          SWchk = 1;
          // Histogram reset
          ResetPH();
          CountH = 0;
        } 
        else {
          SWchk++;
          if ((stmd_logfile) && (stmd_debug)) {
            fprintf(logfile,"  STMD: f= %f  Swchk= %i T= %f\n",f,SWchk,T);
            fprintf(screen,"  STMD: f= %f  Swchk= %i T= %f\n",f,SWchk,T);
          }
        }
      } // if (f_flag == 1)

      if (f_flag > 1) {
        if ((stmd_logfile) && (stmd_debug)) {
          fprintf(logfile,"  STMD: istep= %i  TSC2= %i\n",istep,TSC2);
          fprintf(screen,"  STMD: istep= %i  TSC2= %i\n",istep,TSC2);
        }
        if (STG == 3) // dont reduce if STG4
          f = sqrt(f);
        df = log(f) * 0.5 / bin;
        if ((stmd_logfile) && (stmd_debug)) {
          fprintf(logfile,"  STMD f-UPDATE: f= %f  df= %f\n",f,df);
          fprintf(screen,"  STMD f-UPDATE: f= %f  df= %f\n",f,df);
        }
      
        // Histogram reset
        ResetPH();
        CountH = 0;
      } // if (f_flag > 1)

      // Check stage 3
      if (f <= finFval) STG = 4;

      /*
      // Production run: Hist Output STMD
      if ((o == 0) && (comm->me == 0)) {
      fprintf(fp_whpnm,"### STMD Step=%d: bin E phist tot_hist\n",istep);
        for (int i=0; i<N; i++)
          fprintf(fp_whpnm,"%i %f %i %i\n", i, (i*bin)+Emin,PROH[i],Htot[i]);
        fprintf(fp_whpnm,"\n\n");
      }
      */
    } // if ((m == 0)
  } // if (STG >= 3)

  // STG2 START: Check histogram and modify f value on STG2
  // If STMD, run until histogram is flat, then reduce f value 
  // else if RESTMD, reduce every TSC2 steps
  if (STG == 2) {

    int m = istep % TSC2;
    if (m == 0) {
      if ((stmd_logfile) && (stmd_debug)) {
        fprintf(logfile,"  STMD: istep= %i  TSC2= %i\n",istep,TSC2);
        fprintf(screen,"  STMD: istep= %i  TSC2= %i\n",istep,TSC2);
      }

      // No f-reduction, simulate at initf only!
      if (f_flag == 0) { 
        ResetPH();
        CountH = 0;
      }

      // Standard f-reduction as HCHK() every m steps
      if (f_flag == 1) { 
        HCHK();
        if ((stmd_logfile) && (stmd_debug)) {
          fprintf(logfile,"  STMD: SWfold= %i SWf= %i\n",SWfold,SWf);
          fprintf(screen,"  STMD: SWfold= %i SWf= %i\n",SWfold,SWf);
        }
        // f value update
        if (SWfold != SWf) {
          f = sqrt(f);
          df = log(f) * 0.5 / bin;
          if ((stmd_logfile) && (stmd_debug)) {
            fprintf(logfile,"  STMD f-UPDATE: f= %f  SWf= %i  df= %f\n",f,SWf,df);
            fprintf(screen,"  STMD f-UPDATE: f= %f  SWf= %i  df= %f\n",f,SWf,df);
          }
          SWchk = 1;
          ResetPH();
          CountH = 0;
        }
        else SWchk++;

        if ((stmd_logfile) && (stmd_debug)) {
          fprintf(logfile,"  STMD RESULTS: totCi= %i  f= %f  SWf= %i  SWchk= %i  "
              "STG= %i\n",totCi,f,SWf,SWchk,STG);
          fprintf(screen,"  STMD RESULTS: totCi= %i  f= %f  SWf= %i  SWchk= %i  "
              "STG= %i\n",totCi,f,SWf,SWchk,STG);
        }
        if (f <= pfinFval) {
          STG = 3;
          CountPH = 0;
          SWchk = 1;
          ResetPH();
          CountH = 0;
        }
      } // if (f_flag == 1)

      if (f_flag == 2) { 
      // Reduce as sqrtf every m steps
        if (istep != 0) {
          f = sqrt(f);
          df = log(f) * 0.5 / bin;
        }

        ResetPH();
        CountH = 0;
      } // if (f_flag == 2)

      // Reduce f by constant every m steps
      // otherwise, reduce by sqrt(f) if too small
      if (f_flag == 3) { 
        double reduce_val = 0.1;
        if (istep != 0) {
          if (f > (1+(2*reduce_val)))
            f = f - (reduce_val*f);
          else f = sqrt(f);
        }
        df = log(f) * 0.5 / bin;
        ResetPH();
        CountH = 0;
      } // if (f_flag == 3)

      // Reduce df by constant every m steps
      if (f_flag == 4) {
        double reduce_val = 0.01; // 1% reduction
        if (istep != 0) {
          df = df - (df * reduce_val);
          f = exp(2 * bin * df);
        }
      } // if (f_flag == 4)

      if (f <= 1.0)
        error->all(FLERR,"f-value is less than unity");

      if ((stmd_logfile) && (stmd_debug) && (f_flag > 1)) {
	      fprintf(logfile,"  STMD f-UPDATE: f= %f  df= %f\n",f,df);
	      fprintf(screen,"  STMD f-UPDATE: f= %f  df= %f\n",f,df);
      }

      if ((f <= pfinFval) && (f_flag > 1)) {
        STG = 3;
        CountPH = 0;
      }

    } // if (m == 0)
  } // if (STG == 2)

  // STG1 START: Digging and chk stage on STG1
  // Run until lowest temperature sampled
  if (STG == 1) {
    int m = istep % TSC1;
    if ((m == 0) && (istep != 0)) {
      if (stmd_logfile)
        fprintf(logfile,"  STMD DIG: istep=%i  TSC1=%i Tlow=%f\n",istep,TSC1,T);
      if (stmd_screen)
        fprintf(screen,"  STMD DIG: istep=%i  TSC1=%i Tlow=%f\n",istep,TSC1,T);

      dig();
      TCHK();

      // Histogram reset
      if (STG > 1) {
        ResetPH();
        CountH = 0;
      }
    } // if (m == 0) 
  } // if (STG == 1) 

  if ((stmd_logfile) && (stmd_debug)) {
    fprintf(logfile,"STMD NEXT STG= %i\n",STG);
    fprintf(screen,"STMD NEXT STG= %i\n",STG);
  }
}

/* ---------------------------------------------------------------------- */

double FixStmd::compute_scalar()
{
  return T*ST;
}

/* ---------------------------------------------------------------------- */

double FixStmd::compute_vector(int i)
{
  // Returns useful info for STMD
  double xx;
  if      (i == 0) xx = static_cast<double>(STG);             // Current STG
  else if (i == 1) xx = static_cast<double>(N);               // Number of bins
  else if (i == 2) xx = static_cast<double>(BinMin);          // Lower limit of energy: Emin
  else if (i == 3) xx = static_cast<double>(BinMax);          // Upper limit of energy: Emax
  else if (i == 4) xx = static_cast<double>(curbin);          // Last sampled bin
  else if (i == 5) xx = bin;                                  // Bin spacing of energy: \Delta
  else if (i == 6) xx = df;                                   // df-value
  else if (i == 7) xx = Gamma;                                // force scalling factor
  else if (i == 8) xx = sampledE;                             // Energy/Enthalpy sampled in curbin

  return xx;
}

/* ---------------------------------------------------------------------- */

double FixStmd::compute_array(int i, int j)
{
  // Returns data from arrays
  double xx;
  if      (i == 0) xx = (j*bin)+Emin;    // Binned Energies
  else if (i == 1) xx = Y2[j];           // Histogram of temperature 1/T*j
  else if (i == 2) xx = Hist[j];         // Histogram of energies*j
  else if (i == 3) xx = PROH[j];         // Production histogram*j

  return xx;
}

/* ---------------------------------------------------------------------- */
/*     Enable control to reset things in restart file via fix_modify      */
/* ---------------------------------------------------------------------- */

int FixStmd::modify_param(int narg, char **arg)
{
  // Reset production histogram to 0
  if (strcmp(arg[0],"hist_reset") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");

    if (strcmp(arg[1],"yes") == 0)
      hist_flag = 1;
    else
      error->all(FLERR,"Illegal fix_modify command");
    return 2;
  }
  
  // Reset dfvalue, must be >=0. (=0 means Ts does not update)
  // df will take value from LAMMPS input, STG is NOT reset
  else if (strcmp(arg[0],"dfval") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (strcmp(arg[1],"yes") == 0)
      freset_flag = 1;
    else
      error->all(FLERR,"Illegal fix_modify command");

    f = exp(df * 2 * bin);
    return 2;
  }

  return 0;

}

/* ---------------------------------------------------------------------- */
/*     Extract scale factor                                               */
/* ---------------------------------------------------------------------- */

void *FixStmd::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"scale_stmd") == 0) {
    return &Gamma;
  }
  if (strcmp(str,"sampledE") == 0) {
    return &sampledE;
  }
  return NULL;
}
