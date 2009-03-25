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

// APP as a library
// C or Fortran style interface
// new APP-specific functions can be added

#include "mpi.h"
#include "library.h"
#include "app.h"
#include "app_lattice.h"
#include "comm_lattice.h"
#include "input.h"

using namespace SPPARKS_NS;

/* ----------------------------------------------------------------------
   create an instance of SPPARKS and return pointer to it
------------------------------------------------------------------------- */

void spparks_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  SPPARKS *spk = new SPPARKS(argc,argv,communicator);
  *ptr = (void *) spk;
}

/* ----------------------------------------------------------------------
   destruct an instance of SPPARKS
------------------------------------------------------------------------- */

void spparks_close(void *ptr)
{
  SPPARKS *spk = (SPPARKS *) ptr;
  delete spk;
}

/* ----------------------------------------------------------------------
   process an input script in filename str
------------------------------------------------------------------------- */

void spparks_file(void *ptr, char *str)
{
  SPPARKS *spk = (SPPARKS *) ptr;
  spk->input->file(str);
}

/* ----------------------------------------------------------------------
   process a single input command in str
------------------------------------------------------------------------- */

char *spparks_command(void *ptr, char *str)
{
  SPPARKS *spk = (SPPARKS *) ptr;
  return spk->input->one(str);
}

/* ----------------------------------------------------------------------
   add SPPARKS-specific library functions
   all must receive SPPARKS pointer as argument
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   extract a pointer to an internal SPPARKS value or data structure
   name = desired quantity, e.g. lattice or nlocal
   returns a void pointer which the caller can cast to the desired data type
   returns a NULL if app does not recognize the name
------------------------------------------------------------------------- */

void *spparks_extract(void *ptr, char *name)
{
  SPPARKS *spk = (SPPARKS *) ptr;
  return spk->app->extract(name);
}

/* ----------------------------------------------------------------------
   return total energy of system
------------------------------------------------------------------------- */

double spparks_energy(void *ptr)
{
  SPPARKS *spk = (SPPARKS *) ptr;
  if (spk->app->appclass != App::LATTICE) return 0.0;
  AppLattice *applattice = (AppLattice *) spk->app;
  int nlocal = applattice->nlocal;

  applattice->comm->all();

  double etmp = 0.0;
  for (int i = 0; i < nlocal; i++) etmp += applattice->site_energy(i);
  double energy;
  MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,spk->world);
  return energy;
}
