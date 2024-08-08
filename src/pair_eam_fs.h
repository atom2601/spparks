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
   Code based off that pair_eam_fs.cpp file in LAMMPS 
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(eam/fs,PairEAMFS);
// clang-format on
#else

#ifndef SPK_PAIR_EAM_FS_H
#define SPK_PAIR_EAM_FS_H

#include "pair_eam.h"

namespace SPPARKS_NS {

// need virtual public b/c of how eam/fs/opt inherits from it

class PairEAMFS : virtual public PairEAM {
 public:
  PairEAMFS(class SPPARKS *);

  void coeff(int, char **) override;
 
 protected:
  void read_file(char *) override;
  void file2array() override;
  int he_flag;
};

} // namespace SPPARKS_NS

#endif
#endif