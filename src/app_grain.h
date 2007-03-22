/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef APP_GRAIN_H
#define APP_GRAIN_H

#include "stdio.h"
#include "app.h"

namespace SPPARKS {

class AppGrain : public App {
 public:
  AppGrain(class SPK *, int, char **);
  virtual ~AppGrain();
  virtual void init();
  virtual void input(char *, int, char **);
  virtual void run(int, char **);
  // 2D Hamiltonian functions
  int energy_site(const int, const int, const int);
  void update_mask(char**, const int, const int);
  // 3D Hamiltonian functions
  int energy_site(const int, const int, const int, const int);
  void update_mask(char***, const int, const int, const int);

 protected:
  int me,nprocs;
  int ntimestep;
  int dimension;
  int nx_global,ny_global,nz_global;  // size of global lattice [0,Nglobal-1]
                                      // global lattice has no ghosts
  int nx_local,ny_local,nz_local;     // size of local lattice [0,Nlocal+1]
                                      // local lattice includes ghosts
  int nx_offset,ny_offset,nz_offset;  // global indices (0:Nglobal-1) of
                                      //   my lower-left owned cell (1,1)
                                      // So, for single proc, offsets are zero
  int ***lattice;                     // local lattice with ghost cells
                                      //   allocated to Nlocal+2 in each dim
                                      //   owned cells indexed 1 to Nlocal
                                      //   ghost cells = 0 and Nlocal+1
  int **lat_2d, ***lat_3d;            // 2D and 3D pointers to lattice

  int nx_procs,ny_procs,nz_procs;     // # of procs in each dim
  int procwest,proceast;               // neighbor procs
  int procsouth,procnorth;
  int procdown,procup;    
  int nspins;                          // # of possible spins
  int nsweep;                          // # of sweeps to perform
  int seed;                            // random number generator seed
  int nstats,stats_next;
  int ndump,dump_next;
  FILE *fp;
  double temperature;
  int *dumpbuf;
  int maxdumpbuf;

  class RandomPark *random;

  void iterate();
  void stats();
  void dump_header();
  void dump();
  void dump_detailed(char*);

  void procs2lattice_2d();
  void procs2lattice_3d();
  void set_stats(int, char **);
  void set_dump(int, char **);
  void set_temperature(int, char **);
};

}

#endif
