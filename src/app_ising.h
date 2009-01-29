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

#ifndef APP_ISING_H
#define APP_ISING_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppIsing : public AppLattice {
 public:
  AppIsing(class SPPARKS *, int, char **);
  ~AppIsing();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int rejectstyle;
  int *sites;

  void site_event_rejection_single(int, class RandomPark *);
  void site_event_rejection_double(int, class RandomPark *);
};

}

#endif
