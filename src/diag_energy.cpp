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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"
#include "diag_energy.h"
#include "app_lattice.h"
#include "comm_lattice.h"

using namespace SPPARKS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

DiagEnergy::DiagEnergy(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass != App::LATTICE)
    error->all("Diag style incompatible with app style");
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::init()
{
  applattice = (AppLattice *) app;
  nlocal = applattice->nlocal;
  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::compute()
{
  applattice->comm->all();

  double etmp = 0.0;
  for (int i = 0; i < nlocal; i++) etmp += applattice->site_energy(i);
  MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats(char *strtmp)
{
  sprintf(strtmp,"%10g",energy);
}

/* ---------------------------------------------------------------------- */

void DiagEnergy::stats_header(char *strtmp)
{
  sprintf(strtmp,"%10s","Energy");
}
