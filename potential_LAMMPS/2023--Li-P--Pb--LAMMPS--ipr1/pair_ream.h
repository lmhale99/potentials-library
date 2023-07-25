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

#ifdef PAIR_CLASS

PairStyle(ream,PairREAM)

#else

#ifndef LMP_PAIR_REAM_H
#define LMP_PAIR_REAM_H

#include <stdio.h>
#include "pair.h"

namespace LAMMPS_NS {


class PairREAM : public Pair {
 public:
  friend class FixSemiGrandCanonicalMC;   // Alex Stukowski option

  // public variables so USER-ATC package can access them

  double cutmax;

  // potentials as array data

  int nrho,nP,nr;
  int nfrho,nrhor,nz2r,nRp;
  double **frho,**rhor,**z2r,**Rp,**gr;
  int *type2frho,**type2rhor,**type2z2r;

  // potentials in spline form used for force computation

  double dr,rdr,drho,dP,rdrho,rdP,rhomax,Pmax;
  double ***rhor_spline,***frho_spline,***z2r_spline;
  double ***Rp_spline,***gr_spline;

  PairREAM(class LAMMPS *);
  virtual ~PairREAM();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void swap_eam(double *, double **);

 protected:
  int comm_tag;

  int nmax;                   // allocated size of per-atom arrays
  double cutforcesq;

  // per-atom arrays

  double *rho,*prho,*envf,*fp;

  // potentials as file data

  int *map;                   // which element each atom type maps to

  struct FsX {
    char **elements;
    int nelements,nrho,nP,nr;
    double drho,dP,dr,cut;
    double *mass;
    double **frho,***rhor,***z2r;
    double ***Rp,***gr;
  };
  FsX *fsx;

  virtual void allocate();
  virtual void array2spline();
  void interpolate(int, double, double *, double **);
  void grab(FILE *, int, double *);

  virtual void read_file(char *);
  virtual void file2array();

  virtual double compute_env_interaction(int, int);
  virtual double func_emod(int, int, double);
  virtual double func_emodp(int, int, double);
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

E: Cannot open REAM potential file %s

The specified REAM potential file cannot be opened.  Check that the
path and name are correct.

*/
