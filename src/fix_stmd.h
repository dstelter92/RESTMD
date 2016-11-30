/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(stmd,FixStmd)

#else

#ifndef LMP_FIX_STMD_H
#define LMP_FIX_STMD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStmd : public Fix {
 public:
  FixStmd(class LAMMPS *, int, char **);
  ~FixStmd();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  void end_of_step();
  double memory_usage();

  double compute_scalar();
  double compute_array(int, int);
  void modify_fix(int, double *, char *);

  // Public for access by temper_grem
  double * Y2;              // statistical temperature array
  int STG;                  // stage flag
  int N;                    // number of bins
  double T;                 // latest sampled temperature
  double f;                 // current f-value
  double ST;                // kinetic temperature
  double T1, T2;            // scaled temperature cutoffs

 private:
  int bin;                  // binsize
  int RSTFRQ;               // restart and print frequency
  int f_flag;               // determines type of f-reduction
  int TSC1;                 // dig reduction frequency
  int TSC2;                 // hckh() or f-reduction frequency
  int OREST;                // restart flag, 1 to read restart
  int iworld,nworlds;       // world info
  int BinMin,BinMax;        // bin info
  int Count,CountH,CountPH; // histogram counts   
  int totC,totCi;           // total counts
  int SWf,SWchk,SWfold;     // histogram flatness checks
  int QEXPO;                // exponetial energy bins

  int stmd_logfile,stmd_debug;
  int pe_compute_id;

  double Emin,Emax;         // energy range
  double T0;                // kinetic temp
  double TL, TH;            // unscaled lower and upper T cutoff
  double CTmin,CTmax;       // temperature cutoffs
  double CutTmin,CutTmax;
  double finFval,pfinFval;  // f-tolerance for stg 3 and stg 4
  double initf,df;          // initial-f and delta-f
  double HCKtol;            // histogram tolerance when chk flatness
  double Gamma;             // force scaling factor

  char dir_output[256];     // output directory
  char filename_wtnm[256],filename_whnm[256];
  char filename_whpnm[256],filename_orest[256];

  char * id_pe;
  FILE * fp_wtnm, * fp_whnm, * fp_whpnm, * fp_orest;

  double * Y1, * Prob;
  int * Hist, * Htot, * PROH;

  void dig();               // Translation of stmd.f::stmddig()
  int Yval(double);         // Translation of stmd.f::stmdYval()
  void GammaE(double, int); // Translation of stmd.f::stmdGammaE()
  void AddedEHis(int);      // Translation of stmd.f::stmdAddedEHis()
  void EPROB(int);          // Translation of stmd.f::stmdEPROB()
  void ResetPH();           // Translation of stdm.f::stmdResetPH()
  void TCHK();              // Translation of stmd.f::stmdTCHK()
  void HCHK();              // Translation of stmd.f::stmdHCHK()
  void MAIN(int, double);   // Translation of stmd.f::stmdMAIN()

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid f-reduction scheme

Style provided for f-reduction is incorrect. Use hckh, sqrt, constant.

E: Restart file does not exist

Self-explanatory, change oREST flag.

E: Currently expecting run_style verlet

Fix must use run_style verlet

E: Histogram index out of range

Sampled enthalpy was outside of energy specified by input file.
This is usualy caused by either (1) improper inputs or (2) other
problems in the simulation which cause the energy -> infinity.

E: f-value is less than unity

Self-explanatory, f must be *at least* 1.

E: 
*/
