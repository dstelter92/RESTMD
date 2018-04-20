/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mark Sears (SNL)

   Modified for RESTMD by David Stelter (BU) and Chris Knight (ANL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "temper_stmd.h"
#include "universe.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "integrate.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "output.h"
#include "thermo.h"
#include "fix.h"
#include "random_park.h"
#include "finish.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include <fstream>
#include "fix_stmd.h"

using namespace LAMMPS_NS;

//#define TEMPER_DEBUG 1

/* ---------------------------------------------------------------------- */

TemperStmd::TemperStmd(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

TemperStmd::~TemperStmd()
{
  MPI_Comm_free(&roots);
  if (ranswap) delete ranswap;
  delete ranboltz;
  delete [] set_temp;
  delete [] temp2world;
  delete [] world2temp;
  delete [] world2root;
}

/* ----------------------------------------------------------------------
   perform tempering with inter-world swaps
------------------------------------------------------------------------- */

void TemperStmd::command(int narg, char **arg )
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition to temper");
  if (domain->box_exist == 0)
    error->all(FLERR,"Temper command before simulation box is defined");
  if (narg != 7 && narg != 8)
    error->universe_all(FLERR,"Illegal temper command");

  int nsteps = force->inumeric(FLERR,arg[0]);
  nevery = force->inumeric(FLERR,arg[1]);
  double temp = force->numeric(FLERR,arg[2]);

  for (whichfix = 0; whichfix < modify->nfix; whichfix++)
    if (strcmp(arg[3],modify->fix[whichfix]->id) == 0) break;
  if (whichfix == modify->nfix)
    error->universe_all(FLERR,"Tempering fix ID is not defined");

  // Set pointer to stmd fix
  fix_stmd = (FixStmd*)(modify->fix[whichfix]);

  seed_swap = force->inumeric(FLERR,arg[4]);
  seed_boltz = force->inumeric(FLERR,arg[5]);

  // Exchange flag, 0 if swap off (run many replicas without exchange)
  // 1 if swap on.
  if (strcmp(arg[6],"off") == 0)
    EX_flag = 0;
  else if (strcmp(arg[6],"on") == 0)
    EX_flag = 1;
  else
    error->all(FLERR,"RESTMD: illegal exchange option");

  if (fix_stmd->ST != temp)
    error->universe_all(FLERR,"Kinetic temperatures not "
        "the same, use homogeneous temperature control");

  my_set_temp = universe->iworld;
  if (narg == 8) my_set_temp = force->inumeric(FLERR,arg[7]);

  // swap frequency must evenly divide total # of timesteps
  if (nevery == 0)
    error->universe_all(FLERR,"Invalid frequency in temper command");
  nswaps = nsteps/nevery;
  if (nswaps*nevery != nsteps)
    error->universe_all(FLERR,"Non integer # of swaps in temper command");

  // fix style must be appropriate for temperature control
  if ((strcmp(modify->fix[whichfix]->style,"stmd") != 0)) 
    error->universe_all(FLERR,"Must use with fix STMD, fix is not valid");

  // setup for long tempering run
  update->whichflag = 1;
  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + nsteps;
  if (update->laststep < 0 || update->laststep > MAXBIGINT)
    error->all(FLERR,"Too many timesteps");

  lmp->init();

  // local storage
  me_universe = universe->me;
  MPI_Comm_rank(world,&me);
  nworlds = universe->nworlds;
  iworld = universe->iworld;
  boltz = force->boltz;

  // Setup Swap information
  int nlocal_values = fix_stmd->N + 2; // length of Y2 array + {TL and TH} + extra?
  int nglobal_values = nlocal_values * (nworlds);
  double *local_values = global_values = NULL;
  memory->create(local_values,nlocal_values,"temper/stmd:local_values");
  memory->create(global_values,nglobal_values,"temper/stmd:global_values");

  // pe_compute = ptr to thermo_pe compute
  // notify compute it will be called at first swap
  int id = modify->find_compute("thermo_pe");
  if (id < 0) error->all(FLERR,"Tempering could not find thermo_pe compute");
  Compute *pe_compute = modify->compute[id];
  pe_compute->addstep(update->ntimestep + nevery);

  // create MPI communicator for root proc from each world
  int color;
  if (me == 0) color = 0;
  else color = 1;
  MPI_Comm_split(universe->uworld,color,0,&roots);

  // RNGs for swaps and Boltzmann test
  // warm up Boltzmann RNG
  if (seed_swap) ranswap = new RanPark(lmp,seed_swap);
  else ranswap = NULL;
  ranboltz = new RanPark(lmp,seed_boltz + me_universe);
  for (int i = 0; i < 100; i++) ranboltz->uniform();

  // world2root[i] = global proc that is root proc of world i
  world2root = new int[nworlds];
  if (me == 0)
    MPI_Allgather(&me_universe,1,MPI_INT,world2root,1,MPI_INT,roots);
  MPI_Bcast(world2root,nworlds,MPI_INT,0,world);

  // create static list of set temperatures
  // allgather tempering arg "temp" across root procs
  // bcast from each root to other procs in world
  // leave this in place in the case of inhomogeous temperature
  // control for RESTMD
  set_temp = new double[nworlds];
  if (me == 0) MPI_Allgather(&temp,1,MPI_DOUBLE,set_temp,1,MPI_DOUBLE,roots);
  MPI_Bcast(set_temp,nworlds,MPI_DOUBLE,0,world);

  // create world2temp only on root procs from my_set_temp
  // create temp2world on root procs from world2temp,
  // then bcast to all procs within world
  world2temp = new int[nworlds];
  temp2world = new int[nworlds];
  if (me == 0) {
    MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
    for (int i = 0; i < nworlds; i++) temp2world[world2temp[i]] = i;
  }
  MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);

  // if restarting tempering, reset temp target of Fix to current my_set_temp
  // This should be handled by fix_stmd
  /*
  if (narg == 8) {
    double new_temp = set_temp[my_set_temp];
    modify->fix[whichfix]->reset_target(new_temp);
  }
  */

  // setup tempering runs
  int which,partner,swap,partner_set_temp,partner_world;
  double pe,pe_partner,boltz_factor;

  int stg_flag = 0;
  int stg_flag_me = 0;
  if (fix_stmd->STG == 1) stg_flag_me = 1;

  MPI_Reduce(&stg_flag_me,&stg_flag,1,MPI_INT,MPI_SUM,0,universe->uworld);

  if ((me_universe == 0) && (stg_flag > (universe->nprocs - nworlds)))
    error->universe_warn(FLERR,"RESTMD still in STAGE1, ensure exchanges "
        "turned off");

  if (me_universe == 0 && universe->uscreen)
    fprintf(universe->uscreen,"Setting up RESTMD ...\n");

  update->integrate->setup(1);

  if (me_universe == 0) {
    if (universe->uscreen) {
      fprintf(universe->uscreen,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->uscreen," T%d",i);
      fprintf(universe->uscreen,"\n");
    }
    if (universe->ulogfile) {
      fprintf(universe->ulogfile,"Step");
      for (int i = 0; i < nworlds; i++)
        fprintf(universe->ulogfile," T%d",i);
      fprintf(universe->ulogfile,"\n");
    }
    print_status();
  }

  timer->init();
  timer->barrier_start();

  for (int iswap = 0; iswap < nswaps; iswap++) {

    // run for nevery timesteps
    update->integrate->run(nevery);

    // compute PE
    // notify compute it will be called at next swap
    pe = pe_compute->compute_scalar();
    pe_compute->addstep(update->ntimestep + nevery);

    // Get fix stmd information
    current_STG = fix_stmd->STG;
    T_me = (fix_stmd->T)*(fix_stmd->ST);

    // which = which of 2 kinds of swaps to do (0,1)
    if (!ranswap) which = iswap % 2;
    else if (ranswap->uniform() < 0.5) which = 0;
    else which = 1;

    // partner_set_temp = which set temp I am partnering with for this swap
    if (which == 0) {
      if (my_set_temp % 2 == 0) partner_set_temp = my_set_temp + 1;
      else partner_set_temp = my_set_temp - 1;
    } else {
      if (my_set_temp % 2 == 1) partner_set_temp = my_set_temp + 1;
      else partner_set_temp = my_set_temp - 1;
    }

    // partner = proc ID to swap with
    // if partner = -1, then I am not a proc that swaps
    partner = -1;
    if (me == 0 && partner_set_temp >= 0 && partner_set_temp < nworlds) {
      partner_world = temp2world[partner_set_temp];
      partner = world2root[partner_world];
    }

    // swap with a partner, only root procs in each world participate
    // RESTMD Acceptance Criteria
    swap = 0;
    if (partner != -1) {
      if (me_universe > partner) {
        MPI_Send(&pe,1,MPI_DOUBLE,partner,0,universe->uworld);
        MPI_Send(&T_me,1,MPI_DOUBLE,partner,0,universe->uworld);
      }
      else {
        MPI_Recv(&pe_partner,1,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);
        MPI_Recv(&T_partner,1,MPI_DOUBLE,partner,0,universe->uworld,MPI_STATUS_IGNORE);
      }

      if (me_universe < partner) {
        boltz_factor = (pe_partner - pe)*(1.0/(boltz*T_partner) - 1.0/(boltz*T_me));
        if (boltz_factor >= 0.0) swap = 1;
        else if (ranboltz->uniform() < exp(boltz_factor)) swap = 1;
      }

      // Check what stage, if STG1, no swap
      //if (current_STG == 1) swap = 0; // warning instead...

      if (me_universe < partner)
        MPI_Send(&swap,1,MPI_INT,partner,0,universe->uworld);
      else
        MPI_Recv(&swap,1,MPI_INT,partner,0,universe->uworld,MPI_STATUS_IGNORE);

      if (EX_flag == 0) swap = 0; //If 0, exchanges turned off

#ifdef TEMPER_DEBUG
      if ((me_universe < partner) && (universe->uscreen)) {
        printf("SWAP %d & %d: yes = %d, T = %d %d, PEs = %g %g, Bz = %g %g rand = %g\n",me_universe,partner,swap,my_set_temp,partner_set_temp,pe,pe_partner,boltz_factor,exp(boltz_factor),ranboltz->uniform());
        printf("RESTMD SWAP: N = %d, STG = %d, T_s = %f %f, f = %f\n",fix_stmd->N,current_STG,T_me,T_partner,fix_stmd->f);
      }
#endif

    }

    // bcast swap result to other procs in my world
    MPI_Bcast(&swap,1,MPI_INT,0,world);

    // get information that is being swapped, pack into local_values, gather
    // then bcast to all worlds. All procs pack values for walker into local array
    for (int i=0; i<fix_stmd->N; i++) 
      local_values[i] = fix_stmd->Y2[i];
    local_values[fix_stmd->N] = fix_stmd->T1; //TLOW
    local_values[fix_stmd->N+1] = fix_stmd->T2; //THIGH

    // Gather all local_values from replicas
    if (me == 0) 
      MPI_Allgather(local_values,nlocal_values,MPI_DOUBLE,\
          global_values,nlocal_values,MPI_DOUBLE,roots);

    // Share global_values with universe
    MPI_Bcast(global_values,nglobal_values,MPI_DOUBLE,0,world);

    // if my world swapped, all procs in world reset variables in fix_stmd
    if (swap) {
      // Swap values and unpack from global_values
      const int indx = partner_world*nlocal_values;
      memcpy(&local_values[0],&global_values[indx],nlocal_values*sizeof(double));
      // update fix_stmd with swapped values
      for (int i=0; i<fix_stmd->N; i++) 
        fix_stmd->Y2[i] = local_values[i];
      fix_stmd->T1 = local_values[fix_stmd->N];
      fix_stmd->T2 = local_values[fix_stmd->N+1];
    } // if swap

    // update my_set_temp and temp2world on every proc
    // root procs update their value if swap took place
    // allgather across root procs
    // bcast within my world
    if (swap) my_set_temp = partner_set_temp;
    if (me == 0) {
      MPI_Allgather(&my_set_temp,1,MPI_INT,world2temp,1,MPI_INT,roots);
      for (int i=0; i<nworlds; i++) temp2world[world2temp[i]] = i;
    }
    MPI_Bcast(temp2world,nworlds,MPI_INT,0,world);

    // write stmd temperature files after swap
    fix_stmd->write_temperature();

    // write stmd restart files after swap
    fix_stmd->write_orest();

    // print out current swap status
    if (me_universe == 0) print_status();
  }

  timer->barrier_stop();

  update->integrate->cleanup();

  Finish finish(lmp);
  finish.end(1);

  update->whichflag = 0;
  update->firststep = update->laststep = 0;
  update->beginstep = update->endstep = 0;
}

/* ----------------------------------------------------------------------
   proc 0 prints current tempering status
------------------------------------------------------------------------- */

void TemperStmd::print_status()
{
  if (universe->uscreen) {
    fprintf(universe->uscreen,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->uscreen," %d",world2temp[i]);
    fprintf(universe->uscreen,"\n");
  }
  if (universe->ulogfile) {
    fprintf(universe->ulogfile,BIGINT_FORMAT,update->ntimestep);
    for (int i = 0; i < nworlds; i++)
      fprintf(universe->ulogfile," %d",world2temp[i]);
    fprintf(universe->ulogfile,"\n");
    fflush(universe->ulogfile);
  }
}

