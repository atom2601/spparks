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

#include "math.h"
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_3d.h"
#include "comm_lattice3d.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPotts3d::AppPotts3d(SPPARKS *spk, int& narg, char **arg) : 
  AppLattice3d(spk,narg,arg)
{

  spinfile = NULL;

  // parse arguments

  if (narg < 6) error->all("Illegal app_style command");

  nx_global = atoi(arg[1]);
  ny_global = atoi(arg[2]);
  nz_global = atoi(arg[3]);
  nspins = atoi(arg[4]);
  int seed = atoi(arg[5]);
  random = new RandomPark(seed);
  init_style = RANDOM;

  // parse optional arguments

  int iarg = 6;
  int jarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"random") == 0) {
      init_style = RANDOM;
      iarg ++;
    } else if (strcmp(arg[iarg],"read") == 0) {
      init_style = READ;
      iarg ++;
      if (iarg >= narg) error->all("Illegal app_style command");
      int n = strlen(arg[iarg]) + 1;
      spinfile = new char[n];
      spinfile = strcpy(spinfile,arg[iarg]);
      iarg ++;
    } else {
      // leave other arguments for child app
      arg[jarg] = arg[iarg];
      iarg++;
      jarg++;
    }
  }
  narg = jarg;
}

/* ---------------------------------------------------------------------- */

AppPotts3d::~AppPotts3d()
{
  delete [] spinfile;
}

