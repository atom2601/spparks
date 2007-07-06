/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "universe.h"
#include "memory.h"

using namespace SPPARKS;

/* ----------------------------------------------------------------------
   create & initialize the universe of processors in communicator
------------------------------------------------------------------------- */

Universe::Universe(SPK *spk, MPI_Comm communicator) : SysPtr(spk)
{
  version = (char *)"1 Oct 2006";

  uworld = communicator;
  MPI_Comm_rank(uworld,&me);
  MPI_Comm_size(uworld,&nprocs);

  uscreen = stdout;
  ulogfile = NULL;

  nworlds = 0;
  procs_per_world = NULL;
  root_proc = NULL;
}

/* ---------------------------------------------------------------------- */

Universe::~Universe()
{
  memory->sfree(procs_per_world);
  memory->sfree(root_proc);
}

/* ----------------------------------------------------------------------
   add 1 or more worlds to universe
   str == NULL -> add 1 world with all procs in universe
   str = NxM -> add N worlds, each with M procs
   str = P -> add 1 world with P procs
------------------------------------------------------------------------- */

void Universe::add_world(char *str)
{
  int n,nper;
  char *ptr;

  if (str == NULL) {
    n = 1;
    nper = nprocs;
  } else if ((ptr = strchr(str,'x')) != NULL) {
    *ptr = '\0';
    n = atoi(str);
    nper = atoi(ptr+1);
  } else {
    n = 1;
    nper = atoi(str);
  }

  procs_per_world = 
    (int *) memory->srealloc(procs_per_world,(nworlds+n)*sizeof(int),
				  "universe:procs_per_world");
  root_proc = 
    (int *) memory->srealloc(root_proc,(nworlds+n)*sizeof(int),
				  "universe:root_proc");

  for (int i = 0; i < n; i++) {
    procs_per_world[nworlds] = nper;
    if (nworlds == 0) root_proc[nworlds] = 0;
    else
      root_proc[nworlds] = root_proc[nworlds-1] + procs_per_world[nworlds-1];
    if (me >= root_proc[nworlds]) iworld = nworlds;
    nworlds++;
  }
}

/* ----------------------------------------------------------------------
   check if total procs in all worlds = procs in universe
------------------------------------------------------------------------- */

int Universe::consistent()
{
  int n = 0;
  for (int i = 0; i < nworlds; i++) n += procs_per_world[i];
  if (n == nprocs) return 1;
  else return 0;
}
