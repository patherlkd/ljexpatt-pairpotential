/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.


Luke kristopher Davis @ucl    luke.davis.16@ucl.ac.uk
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/exp/att,PairLEA)

#else

#ifndef LMP_PAIR_LEA_H
#define LMP_PAIR_LEA_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLEA : public Pair {
 public:
  PairLEA(class LAMMPS *);
  virtual ~PairLEA();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);

 protected:
  double cut_global,cut_lj, min_global;
  double **cut,**cut_exp, **cut_exp2;
  double **epsilon,**epsilon_exp,**sigma;
  double **lambda;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  double *cut_respa;
  double **exp_sl,**exp_rcl;
  unsigned int lj_switch, extracted;
  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

E: cut_one > cut_two

This means the inner repulsive cutoff is larger than the attraction cutoff. This cannot be for obvious reasons. The attraction happens beyond the 'excluded volume' region. 

*/
