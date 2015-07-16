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

FixStyle(stmd,FixSTMD)

#else

#ifndef LMP_FIX_STMD_H
#define LMP_FIX_STMD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSTMD : public Fix {
 public:
  FixSTMD(class LAMMPS *, int, char **);
  ~FixSTMD();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  void end_of_step();
  double memory_usage();
  double * Y2;

  double compute_scalar();
  double compute_array(int, int);
  void modify_fix(int, double *, char *);
  //friend class Temper_STMD;

 private:
  // arguments of fix; these variables were initialized in charmm.inp
  int RSTFRQ, bin, PRNFRQ, TSC1, TSC2, OREST;
  double initf, TL, TH, Emin, Emax, ST;
  char dir_output[256]; // optional argument for output directory
  int iworld,nworlds; // world info
  // other variables initialized in stmd.f::stmdcntrl()
  int MODI;

  char filename_wtnm[256], filename_whnm[256], filename_whpnm[256], filename_wresnm[256], filename_iresnm[256], filename_wenm[256], filename_orest[256], filename_irest[256];
  FILE * fp_wtnm, * fp_whnm, * fp_whpnm, * fp_wresnm, * fp_iresnm, * fp_wenm, * fp_orest, * fp_irest;

  double CutTmin, CutTmax, finFval, pfinFval, HCKtol, multi, dymT;
  int QREST, QEXPO, QEXP1;

  int BinMin, BinMax, N;
  int STG, SWf, Count, CountH, totC, totCi, SWchk, CountPH, SWfold;
  double Gamma;

  double * Elist;
  double f, df, T0, T1, T2, CTmin, CTmax, scaleT;

  double * Y1, * Prob;
  int * Hist, * Htot, * PROH;

  int stmd_logfile, stmd_debug;

  void dig();               // Translation of stmd.f::stmddig()
  int Yval(double);         // Translation of stmd.f::stmdYval()
  void GammaE(double, int); // Translation of stmd.f::stmdGammaE()
  void AddedEHis(int);      // Translation of stmd.f::stmdAddedEHis()
  void EPROB(int);          // Translation of stmd.f::stmdEPROB()
  void ResetPH();           // Translation of stdm.f::stmdResetPH()
  void TCHK();              // Translation of stmd.f::stmdTCHK()
  void HCHK();              // Translation of stmd.f::stmdHCHK()
  //void OREST();             // Translation of stmd.f::stmdOREST()
  void MAIN(int, double);   // Translation of stmd.f::stmdMAIN()

  int pe_compute_id;
  char * id_pe;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix setforce does not exist

Self-explanatory.

E: Variable name for fix setforce does not exist

Self-explanatory.

E: Variable for fix setforce is invalid style

Only equal-style variables can be used.

E: Cannot use non-zero forces in an energy minimization

Fix setforce cannot be used in this manner.  Use fix addforce
instead.

*/
