/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Anthony Tom (University of Tennessee - Knoxville)
   Code based off that pair_eam.cpp file in LAMMPS made by
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(eam,PairEAM)
// clang-format on
#else

#ifndef SPK_PAIR_EAM_H
#define SPK_PAIR_EAM_H

#include "pair.h"

namespace SPPARKS_NS{
    
class PairEAM : public Pair {
 public:
  
  // public variables

  double cutmax;

  // potentials as array data

  int nrho, nr;
  int nfrho, nrhor, nz2r;
  double **frho, **rhor, **z2r;
  int *type2frho, **type2rhor, **type2z2r;

  // potentials in spline form used for force computation

  double dr, rdr, drho, rdrho, rhomax, rhomin;
  double ***rhor_spline, ***frho_spline, ***z2r_spline;

  PairEAM(class SPPARKS *);
  ~PairEAM() override;
  double energy(int, int, int *, double **, int *) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);
  void *extract_peratom(const char *, int &);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void swap_eam(double *, double **);

  class CommOffLattice *comm;

 protected:
  double nmax; // allocated size of per-atom arrays
  double cutforcesq;
  double **scale;
  bigint embedstep; // timestep, the embedding term was computed

  virtual void ev_setup(int, int, int alloc = 1);
  void ev_unset();

  // per-atom arrays

  double *rho, *fp;
  int *numforce;

  // potentials as file data

  struct Funcfl {
    char *file;
    int nrho, nr;
    double drho, dr, cut, mass;
    double *frho, *rhor, *zr;
  };
  Funcfl *funcfl;
  int nfuncfl;

  struct Setfl {
    char **elements;
    int nelements, nrho, nr;
    double drho, dr, cut;
    double *mass;
    double **frho, **rhor, ***z2r;
  };
  Setfl * setfl;

  struct Fs
  {
    char **elements;
    int nelements, nrho, nr;
    double drho, dr, cut;
    double *mass;
    double **frho, ***rhor, ***z2r;
  };
  Fs * fs;

  virtual void allocate();
  virtual void array2spline();
  void interpolate(int, double, double *, double **);

  virtual void read_file(char *);
  virtual void file2array();
};

}  // namespace SPPARKS_NS

#endif
#endif