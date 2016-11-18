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

#include "string.h"
#include "stdlib.h"
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
  if (narg < 16 || narg > 17) error->all(FLERR,"Illegal fix stmd command");

  global_freq = 1;
  scalar_flag = 1;
  restart_global = 0;
  restart_peratom = peratom_flag = 0;

  // This is the subset of variables explicitly given in the charmm.inp file
  // If the full set is expected to be modified by a user, then reading 
  //  a stmd.inp file is probably the best mechanism for input.
  //
  //  fix fxstmd all stmd RSTFRQ f Tlo Thi Plo Phi binsize 10000 40000 300 PRNFRQ OREST

  RSTFRQ = atoi(arg[3]);        // Probably a good idea to set this equal to restart value in input
  initf  = atof(arg[4]);
  TL     = atof(arg[5]);
  TH     = atof(arg[6]);
  Emin   = atof(arg[7]);
  Emax   = atof(arg[8]);
  bin    = atoi(arg[9]);
  TSC1   = atof(arg[10]);
  TSC2   = atof(arg[11]);
  ST     = atof(arg[12]);       // This value should be consistent with target temperature of thermostat fix
  PRNFRQ = atoi(arg[13]);
  f_flag = atoi(arg[14]);      // 0 for hckh() to reduce f-val, 1 for constant reduction
  OREST  = atoi(arg[15]);       // 0 for new run, 1 for restart
  
  // If RESTMD, check with temper_stmd to ensure walkers are the same.
  // Do this in temper_stmd...
  
  // Make dir_output hard coded to local dir
  strcpy(dir_output,"./");
  
  Y1 = Y2 = Prob = NULL;
  Hist = Htot = PROH = NULL;

  stmd_logfile = stmd_debug = 0;
  if ((comm->me == 0) && (logfile)) stmd_logfile = 1;

  // DEBUG FLAG
  //stmd_debug = 1;

  fp_wtnm = fp_whnm = fp_whpnm = fp_orest = NULL;
}

/* ---------------------------------------------------------------------- */

FixStmd::~FixStmd()
{
  memory->destroy(Y1);
  memory->destroy(Y2);
  memory->destroy(Hist);
  memory->destroy(Htot);
  memory->destroy(PROH);
  memory->destroy(Prob);
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

  if (OREST) {
    // Check if file exists/has data, otherwise exit!
    char filename[256];
    strcpy(filename,dir_output);
    strcat(filename,"/oREST.");
    strcat(filename,walker);
    strcat(filename,".d");
    strcpy(filename_orest,filename);

    if (fp_orest = fopen(filename, "r")) {
      fclose(fp_orest);
    } else {
      if (nworlds > 1) error->universe_all(FLERR,"RESTMD: Restart file does not exist\n");
      else error->all(FLERR,"STMD: Restart file does not exist\n");
    }
  }

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
    if (!fp_whpnm) {
      strcpy(filename,dir_output);
      strcat(filename,"/WHP.");
      strcat(filename,walker);
      strcat(filename,".d");
      strcpy(filename_whpnm,filename);
      fp_whpnm = fopen(filename,"w");
    }
    if ((!fp_orest) && (!OREST)) {
      strcpy(filename,dir_output);
      strcat(filename,"/oREST.");
      strcat(filename,walker);
      strcat(filename,".d");
      strcpy(filename_orest,filename);
      fp_orest  = fopen(filename,"w");
    }
    if ((!fp_orest) && (OREST)) {
      strcpy(filename,dir_output);
      strcat(filename,"/oREST.");
      strcat(filename,walker);
      strcat(filename,".d");
      strcpy(filename_orest,filename);
      fp_orest  = fopen(filename,"r");
    }
  }

  CutTmin  = 50.0;
  CutTmax  = 50.0;
  finFval  = 1.0000001;
  pfinFval = 1.000001;
  HCKtol   = 0.2;

  QEXPO = 0;

  BinMin = round(Emin / bin);
  BinMax = round(Emax / bin);

  // Exponential energy bin setup
  // if (QEXPO) {
  //   exv = exfB - exfA;
  //   exc = -exv / log(exal) / exbe;
  //   exb = (exeB - exeA) / (exp(exfB / exc) - exp(exfA / exc));
  //   exa = exeB - exb * exp(exfB / exc);
  //   if ((Emin - exa) / exb <= 0.0) error->warning(FLERR,"WARNNING!!! EXPONENT ENE BIN SETUP");
  //   BinMin = round(log((Emin - exa) / exb));
  //   BinMax = round(log((Emax - exa) / exb));
  // }
  
  N = BinMax - BinMin + 1;

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

  f = initf;
  df = log(f) * 0.5 / bin;
  T0 = ST;
  T1 = TL / ST;
  T2 = TH / ST;
  CTmin = (TL + CutTmin) / ST;
  CTmax = (TH - CutTmax) / ST;

  memory->grow(Y1, N, "FixSTMD:Y1");
  memory->grow(Y2, N, "FixSTMD:Y2");
  memory->grow(Hist, N, "FixSTMD:Hist");
  memory->grow(Htot, N, "FixSTMD:Htot");
  memory->grow(PROH, N, "FixSTMD:PROH");
  memory->grow(Prob, N, "FixSTMD:Prob");

  for (int i=0; i<N; i++) {
    Y1[i] = T2;
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
          
      for (int i=0; i<nsize; i++) 
        file >> list[i];

      //int nbins = static_cast<int> (list[k++]);
      STG = static_cast<int> (list[k++]);
      f = list[k++];
      CountH = static_cast<int> (list[k++]);
      //df = list[k++];
      SWf = static_cast<int> (list[k++]);
      SWfold = static_cast<int> (list[k++]);
      SWchk = static_cast<int> (list[k++]);
      Count = static_cast<int> (list[k++]);
      totCi = static_cast<int> (list[k++]);
      CountPH = static_cast<int> (list[k++]);
      //TSC1 = static_cast<int> (list[k++]);
      //TSC2 = static_cast<int> (list[k++]);
      //Gamma = list[k++];
      //T0 = list[k++];
      //ST = list[k++];
      T1 = list[k++];
      T2 = list[k++];
      CTmin = list[k++];
      CTmax = list[k++];
      
      for (int i=0; i<N; i++)
        Y2[i] = list[k++];
      for (int i=0; i<N; i++)
        Htot[i] = list[k++];
      for (int i=0; i<N; i++)
        PROH[i] = list[k++];
      
      memory->destroy(list);
    }
  }
  
  // Write values of all paramters to logfile
  if ((stmd_logfile) && (nworlds > 1)) {
    fprintf(logfile,"RESTMD: #replicas= %i  walker= %i\n",nworlds,iworld);
    fprintf(screen,"RESTMD: #replicas= %i  walker= %i\n",nworlds,iworld);
    }

  if ((stmd_logfile) && (nworlds == 1)) {
    fprintf(logfile,"STMD Check initial values...\n");
    fprintf(logfile,"STMD STAGE= %i #bins= %i  binsize %i\n",STG, N, bin); // diffE was included in stmd.f, but don't know what that is
    fprintf(screen,"STMD Check initial values...\n");
    fprintf(screen,"STMD STAGE= %i #bins= %i  binsize= %i\n",STG, N, bin);
  }
    
  if ((stmd_logfile) && (stmd_debug)) {
    //fprintf(logfile,"STMD Yold(Y1)= ");
    //for (int i=0; i<N; i++) fprintf(logfile," %f",Y1[i]);
    //fprintf(logfile,"\n");
    fprintf(logfile,"STMD Temperature (Y2)= ");
    for (int i=0; i<N; i++) 
      fprintf(logfile," %f",Y2[i]);
    fprintf(logfile,"\n");
  }
}

/* ---------------------------------------------------------------------- */

void FixStmd::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet")) {
    post_force(vflag);

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
  double potE = modify->compute[pe_compute_id]->compute_scalar(); 

  // Master rank will compute scaling factor and then Bcast to world
  MAIN(update->ntimestep,potE);

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

int FixStmd::Yval(double potE)
{
  int i = round(potE / double(bin)) - BinMin + 1;

  if ((i<1) || (i>N-1)) {
    fprintf(screen,"Error in Yval: potE= %f  bin= %i  i= %i\n",potE,bin,i);
    fprintf(logfile,"Error in Yval: potE= %f  bin= %i  i= %i\n",potE,bin,i);
    if (nworlds > 1) error->universe_all(FLERR,"RESTMD: Histogram index out of range");
    else error->one(FLERR,"STMD: Histogram index out of range");
  }

  Y2[i+1] = Y2[i+1] / (1.0 - df * Y2[i+1]);
  Y2[i-1] = Y2[i-1] / (1.0 + df * Y2[i-1]);

  if (Y2[i-1] < T1) 
    Y2[i-1] = T1;
  if (Y2[i+1] > T2) 
    Y2[i+1] = T2;

  return i;
}

/* ---------------------------------------------------------------------- */

void FixStmd::GammaE(double potE, int indx)
{
  const int i  = indx;
  const int im = indx - 1;
  const int ip = indx + 1;

  const double e = potE - double( round(potE / double(bin)) * bin );

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
  if ((stmd_logfile) && (stmd_debug))
    fprintf(logfile,"Stmd TCHK: T1= %f (%f K)  Y2[0]= %f (%f K)\n",T1,T1*ST,Y2[0],Y2[0]*ST);
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

  if ((stmd_logfile) && (stmd_debug)) 
    fprintf(logfile,"STMD CHK HIST: icnt= %i  aveH= %f  N= %i\n",icnt,aveH,N);
  if (icnt==0) return;

  aveH = aveH / double(icnt);
  
  double eval;
  for (int i=0; i<N; i++) {
    if ((Y2[i] > CTmin) && (Y2[i] < CTmax) ) {
      eval = abs(double(Hist[i] - aveH) / aveH);
      if (eval > HCKtol) ichk++;
      if ((stmd_logfile) && (stmd_debug))
        fprintf(logfile,"STMD CHK HIST: totCi= %i  i= %i  eval= %f  HCKtol= %f  " 
            "ichk= %i  Hist[i]= %i\n",totCi,i,eval,HCKtol,ichk,Hist[i]);
    }
  }

  if (ichk < 1) SWf = SWf + 1;
}

/* ---------------------------------------------------------------------- */

void FixStmd::MAIN(int istep, double potE)
{
  Count = istep;
  totCi++;

  if (STG >= 3) CountPH++;

  if ((stmd_logfile) && (stmd_debug))
    fprintf(logfile,"STMD STAGE %i\n",STG);

  // Statistical Temperature Update
  int stmdi = Yval(potE);
  
  // Gamma Update
  GammaE(potE,stmdi);

  if ((stmd_logfile) && (stmd_debug))
    fprintf(logfile,"STMD: totCi= %i Count= %i Gamma= %f stmdi= %i\n "
        "T= %f",totCi,Count,Gamma,stmdi,T);

  // Histogram Update
  AddedEHis(stmdi);
  CountH++;

  // Add to Histogram for production run
  if (STG >= 3) {
    PROH[stmdi]++;
    CountPH++;
  }

  // Hist Output
  int m = istep % PRNFRQ;
  if ((m == 0) && (comm->me == 0)) {
    for (int i=0; i<N; i++) fprintf(fp_whnm,"%i %f %i %i %f %i %i %f"
        "\n", i, (i*bin)+Emin, Hist[i], Htot[i], Y2[i], CountH, totCi, f);
    fprintf(fp_whnm,"\n\n");
  }

  // Production Run if STG >= 3
  // STG3 START: Check histogram and further reduce f until <= 1.0000001
  if (STG >= 3) {
    m = istep % TSC2;

    // STMD, reduce f based on histogram check
    if ((m == 0) && (f_flag == 0)) { 
      if ((stmd_logfile) && (stmd_debug))
        fprintf(logfile,"STMD STAGE 3\nSTMD CHK HIST: istep= %i  TSC2= %i\n",istep,TSC2);
      HCHK();
      if ((stmd_logfile) && (stmd_debug)) {
        fprintf(logfile,"STMD STG3 SWfold= %i  SWf= %i\n",SWfold,SWf);
        fprintf(logfile,"STMD STG3 f= %f  SWchk= %i\n",f,SWchk);
      }
      if (SWfold != SWf) {
        if ((stmd_logfile) && (stmd_debug))
          fprintf(logfile,"STMD STG f= %f  df= %f\n",f,df);
        f = sqrt(f);
        df = log(f) * 0.5 / double(bin);
        if ((stmd_logfile) && (stmd_debug)) {
          fprintf(logfile,"STMD STG3 f= %f  SWf= %i  df= %f\n",f,SWf,df);
          fprintf(logfile,"STMD STG3 NEXT STG= %i\n",STG);
        }
        SWchk = 1;
        // Histogram reset
        ResetPH();
        CountH = 0;
      } else {
        SWchk++;
        if ((stmd_logfile) && (stmd_debug))
          fprintf(logfile,"STMD STG3 f= %f  Swchk= %i T= %f\n",f,SWchk,T);
      }
      
      // Check stage 3
      if (f <= finFval) STG = 4;

      // Production run: Hist Output STMD
      m = istep % PRNFRQ;
      if ((m == 0) && (comm->me == 0)) {
        for (int i=0; i<N; i++)
          fprintf(fp_whpnm,"%i %f %i %i %i %f %i %i %f\n", i, (i*bin)+Emin,\
              Hist[i],PROH[i],Htot[i],Y2[i],CountH,CountPH,f);
        fprintf(fp_whpnm,"\n\n");
      }
    } // if ((m == 0) && (f_flag == 0))

    // RESTMD, reduce f to sqrt(f) every TSC2 steps
    if ((m == 0) && (f_flag != 0)) { 
      if ((stmd_logfile) && (stmd_debug))
        fprintf(logfile,"RESTMD STAGE 3\nRESTMD STG3 istep= %i  TSC2= %i\n",istep,TSC2);
      f = sqrt(f);
      df = log(f) * 0.5 / double(bin);
      if ((stmd_logfile) && (stmd_debug)) {
        fprintf(logfile,"RESTMD STG3 f= %f  df= %f\n",f,df);
        fprintf(logfile,"RESTMD STG3 NEXT STG= %i\n",STG);
      }

      // Histogram reset
      ResetPH();
      CountH = 0;

      // Check stage 3
      if (f <= finFval) STG = 4;

      // Production run: Hist Output RESTMD
      m = istep % PRNFRQ;
      if ((m == 0) && (comm->me == 0)) {
        for (int i=0; i<N; i++) 
          fprintf(fp_whpnm,"%i %f %i %i %i %f %i %i %f\n", i, (i*bin+Emin), Hist[i],\
              PROH[i], Htot[i], Y2[i], CountH, CountPH, f);
        fprintf(fp_whpnm,"\n\n");
      }
    } // if ((m == 0) && (f_flag == 1))
  } // if (STG >= 3)

  // STG2 START: Check histogram and modify f value on STG2
  // If STMD, run until histogram is flat, then reduce f value 
  // until <= 1.000001, else if RESTMD, reduce every TSC2 steps
  if (STG == 2) {
    m = istep % TSC2;
    
    // If STMD...
    if ((m == 0) && (f_flag == 0)) { 
      if ((stmd_logfile) && (stmd_debug))
        fprintf(logfile,"STMD STAGE 2\nSTMD STG2: CHK HIST istep= %i  "
            "TSC2= %i\n",istep,TSC2);
      HCHK();
      if ((stmd_logfile) && (stmd_debug))
        fprintf(logfile,"STMD STG2: SWfold= %i SWf= %i\n",SWfold,SWf);

      // F value update
      if (SWfold != SWf) {
        if ((stmd_logfile) && (stmd_debug))
          fprintf(logfile,"STMD STG2: f= %f  df= %f\n",f,df);
        f = sqrt(f);
        df = log(f) * 0.5 / double(bin);
        if ((stmd_logfile) && (stmd_debug)) {
          fprintf(logfile,"STMD STG2: f= %f  SWf= %i  df= %f\n",f,SWf,df);
          fprintf(logfile,"STMD STG2: STG= %i\n",STG);
        }
        
        SWchk = 1;
        ResetPH();
        CountH = 0;
      } else SWchk++;
      
      if ((stmd_logfile) && (stmd_debug))
        fprintf(logfile,"STMD SG2 RESULTS: totCi= %i  f= %f  SWf= %i  SWchk= %i  "
            "STG= %i\n",totCi,f,SWf,SWchk,STG);

      if (f <= pfinFval) {
        STG = 3;
        CountPH = 0;
        if ((stmd_logfile) && (stmd_debug)) {
          fprintf(logfile,"STMD STG2: f= %f  SWf= %i  df= %f\n",f,SWf,df);
          fprintf(logfile,"STMD STG2: STG= %i T= %f\n",STG,T);
        }
        SWchk = 1;
        ResetPH();
        CountH = 0;
      }      
    } // if ((m == 0) && (f_flag == 0))

    // If RESTMD...
    if ((m == 0) && (f_flag == 1)) { 
        if ((stmd_logfile) && (stmd_debug))
          fprintf(logfile,"RESTMD STAGE 2\nRESTMD STG2: istep= %i  TSC2= %i\n",istep,TSC2);
        if (istep != 0) f = sqrt(f);
        df = log(f) * 0.5 / double(bin);
        if (f <= pfinFval) {
          STG = 3;
          CountPH = 0;
        }
	    if ((stmd_logfile) && (stmd_debug)) {
	      fprintf(logfile,"RESTMD STG2: f= %f  df= %f\n",f,df);
	      fprintf(logfile,"RESTMD STG2: STG= %i\n",STG);
      }
      
      ResetPH();
      CountH = 0;
	  } // if ((m == 0) && (f_flag == 1)) 
    
    // RESTMD, alternative f-reduction scheme
    if ((m == 0) && (f_flag == 2)) { 
        if ((stmd_logfile) && (stmd_debug))
          fprintf(logfile,"RESTMD STAGE 2\nRESTMD STG2: istep= %i  TSC2= %i\n",istep,TSC2);
        // Reduce f by 10%, if above cutoff
        // otherwise, reduce by sqrt(f) 
        double reduce_val = 0.1;
        if (istep != 0) {
          if (f > (1+(2*reduce_val)))
            f = f - (reduce_val*f);
          else f = sqrt(f);
        }
        if (f <= 1.0)
          error->all(FLERR,"f-value is less than unity");

        df = log(f) * 0.5 / double(bin);
        if (f <= pfinFval) {
          STG = 3;
          CountPH = 0;
        }
	    if ((stmd_logfile) && (stmd_debug)) {
	      fprintf(logfile,"RESTMD STG2: f= %f  df= %f\n",f,df);
	      fprintf(logfile,"RESTMD STG2: STG= %i\n",STG);
      }
      
      ResetPH();
      CountH = 0;
	  } // if ((m == 0) && (f_flag == 2)) 

  } // if (STG == 2)

  // STG1 START: Digging and chk stage on STG1
  // Run until lowest temperature sampled
  if (STG == 1) {
    m = istep % TSC1;
    if (m == 0) {
      if (stmd_logfile) {
        fprintf(logfile,"STMD STAGE 1\n");
        fprintf(screen,"STMD STAGE 1\n");
        fprintf(logfile,"STMD STAGE1 DIG: istep= %i  TSC1= %i T= %f\n",istep,TSC1,T);
        fprintf(screen,"STMD STAGE1 DIG: istep= %i  TSC1= %i T= %f\n",istep,TSC1,T);
      }

      dig();
      TCHK();

      if ((stmd_logfile) && (stmd_debug))
        fprintf(logfile,"STMD STG1: NEXT STG= %i\n",STG);

      // Histogram reset
      if (STG > 1) {
        ResetPH();
        CountH = 0;
      }
    } // if (m == 0) 
  } // if (STG == 1) 

  // Yval output
  m = istep % PRNFRQ;
  if ((m == 0) && (comm->me == 0)) {
    for (int i=0; i<N; i++) 
      fprintf(fp_wtnm,"%i %f %f %f %i\n", i,(i*bin)+Emin,Y2[i]*ST,Y2[i],totCi);
    fprintf(fp_wtnm,"\n\n");
  }

  iworld = universe->iworld;
  // Write restart info to external file
  m = istep % RSTFRQ;
  if ((m == 0) && (comm->me == 0) && (Count > 0)) {
    int k = 0;
    int numb = 13;
    int nsize = 3*N + numb;
    double *list;
    memory->create(list,nsize,"stmd:list");

    //list[k++] = N;
    list[k++] = STG;
    list[k++] = f;
    list[k++] = CountH;
    //list[k++] = df;
    list[k++] = SWf;
    list[k++] = SWfold;
    list[k++] = SWchk;
    list[k++] = Count;
    list[k++] = totCi;
    list[k++] = CountPH;
    //list[k++] = TSC1; // Control these by LAMMPS input file
    //list[k++] = TSC2;
    //list[k++] = Gamma;
    //list[k++] = T0;
    //list[k++] = ST;
    list[k++] = T1;
    list[k++] = T2;
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
    sprintf(walker,"%i",iworld);
    strcpy(filename,dir_output);
    strcat(filename,"/oREST.");
    strcat(filename,walker);
    strcat(filename,".d");
    freopen(filename,"w",fp_orest);

    for (int i=0; i<numb; i++) 
      fprintf(fp_orest,"%f\n",list[i]);
    for (int i=numb; i<N+numb; i++) 
      fprintf(fp_orest,"%f ",list[i]);
    fprintf(fp_orest,"\n");
    for (int i=N+numb; i<2*N+numb; i++) 
      fprintf(fp_orest,"%f ",list[i]);
    fprintf(fp_orest,"\n");
    for (int i=2*N+numb; i<nsize; i++) 
      fprintf(fp_orest,"%f ",list[i]);
    fprintf(fp_orest,"\n");

    memory->destroy(list);
  }
}

/* ---------------------------------------------------------------------- */

double FixStmd::compute_scalar()
{
  return Gamma;
}

/* ---------------------------------------------------------------------- */

double FixStmd::compute_array(int i, int j)
{
  // Returns data from arrays
  double xx;
  if      (i == 0) xx = static_cast<double>(BinMax-BinMin+1); // Number of bins
  else if (i == 1) xx = static_cast<double>(BinMin);          // Lower limit of energy: Emin
  else if (i == 2) xx = static_cast<double>(BinMax);          // Upper limit of energy: Emax
  else if (i == 3) xx = static_cast<double>(bin);             // Bin spacing of energy: \Delta
  else if (i == 4) xx = Y2[j];                                // Histogram of temperature 1/T*j

  return xx;
}

/* ---------------------------------------------------------------------- */
/* --- Trigger reinitialization of key arrays after modify_fix called --- */
/* ---------------------------------------------------------------------- */

void FixStmd::modify_fix(int which, double *values, char *notused)
{
  // Sets a specified variable to the input value(s)
  if      (which == 0) BinMin = static_cast<int>(values[0] + 0.5);
  else if (which == 1) BinMax = static_cast<int>(values[0] + 0.5);
  else if (which == 2) bin    = static_cast<int>(values[0] + 0.5);
  else if (which == 3) {
    for (int i=0; i<N; i++) Y2[i] = values[i];
  }
}
