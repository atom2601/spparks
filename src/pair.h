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

#ifndef SPK_PAIR_H
#define SPK_PAIR_H

#include "pointers.h"

namespace SPPARKS_NS {

class Pair : protected Pointers {
 public:
  int ntypes;
  double cutoff;

  Pair(class SPPARKS *);
  virtual ~Pair() {}
  void init();

  virtual void settings(int, char **) = 0;
  virtual void coeff(int, char **) = 0;
  virtual void init_style() {}
  virtual double init_one(int, int) {return 0.0;}
  virtual double energy(int, int, int *, double **, int *) = 0;

   // int comm_forward;          // size of forward communication (0 if none)
   // int comm_reverse;          // size of reverse communication (0 if none)

  int restartinfo;                     // 1 if pair style writes restart info
  int one_coeff;                       // 1 if allows only one coeff * * call
  int manybody_flag;                   // 1 if a manybody potential
  int unit_convert_flag;               // value != 0 indicates support for unit conversion

//   int evflag;
//   int eflag_either, eflag_global, eflag_atom;
//   int vflag_either, vflag_global, vflag_atom, cvflag_atom;

 protected:
  int allocated;                       // 0/1 = whether arrays are allocated
  int **setflag;
  double **cutsq;
  int mix_flag;
  int vflag_fdotr;

  double mix_energy(double, double, double, double);
  double mix_distance(double, double);

  // for mapping of elements to atom types and parameters
  // mostly used for manybody potentials (here the eam)
  // "elements" here will be the different number of sites
  int nelements;   // # of unique elements
  char **elements; // names of unique elements
  int *map;
  void map_element2type(int, char **, bool update_setflag = true); 

   // various lammps functions not used here
   //     void ev_init(int eflag, int vflag, int alloc = 1)
   //   {
   //     if (eflag || vflag)
   //       ev_setup(eflag, vflag, alloc);
   //     else
   //       ev_unset();
   //   }
   //   virtual void ev_setup(int, int, int alloc = 1);
   //   void ev_unset();
};

}

#endif

/* ERROR/WARNING messages:

E: All pair coeffs are not set

Self-explanatory.

*/
