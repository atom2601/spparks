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
#include "random_park.h"

using namespace SPPARKS_NS;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836

/* ---------------------------------------------------------------------- */

RandomPark::RandomPark(int seed_init=0)
{
  seed = seed_init;
}

/* ---------------------------------------------------------------------- 
    reset the seed (if given) and warm up the sequence
 ------------------------------------------------------------------------*/

void RandomPark::init(int seed_init=0)
{
  if (seed_init != 0) {
    seed = seed_init;
  }
  for (int ii = 0; ii < 100; ii++) uniform();
}

/* ----------------------------------------------------------------------
   uniform RN 
------------------------------------------------------------------------- */

double RandomPark::uniform()
{
  int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  double ans = AM*seed;
  return ans;
}

/* ----------------------------------------------------------------------
   integer RN between 1 and N inclusive
------------------------------------------------------------------------- */

int RandomPark::irandom(int n)
{
  int i = (int) (uniform()*n) + 1;
  if (i > n) i = n;
  return i;
}
