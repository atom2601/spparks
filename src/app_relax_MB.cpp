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
#include "string.h"
#include "stdlib.h"
#include "app_relax_MB.h"
#include "potential.h"
#include "pair.h"
#include "random_park.h"
#include "error.h"

#include <iostream>

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppRelaxMB::AppRelaxMB(SPPARKS *spk, int narg, char **arg) : 
  AppOffLattice(spk,narg,arg)
{
  ninteger = 1;
  ndouble = 0;
  allow_kmc = 1;
  allow_rejection = 1;
  dt_sweep = 1.0;

  create_arrays();

  // parse arguments

  if (narg != 2) error->all(FLERR,"Illegal app_style command");

  delta = atof(arg[1]);
  deltasq = delta*delta;
}

/* ---------------------------------------------------------------------- */

AppRelaxMB::~AppRelaxMB() {}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppRelaxMB::grow_app()
{
  type = iarray[0];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppRelaxMB::init_app()
{
  potential->init();
  pair = potential->pair;
  if (pair == NULL) error->all(FLERR,"App relax requires a pair potential");

  delpropensity = pair->cutoff;
  delevent = delta;

  int ntypes = potential->pair->ntypes;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (type[i] < 1 || type[i] > ntypes) flag = 1;
  int flagall = 0;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   compute energy of site
   form a neighbor list first
------------------------------------------------------------------------- */

double AppRelaxMB::site_energy(int i)
{
  neighbor(i,pair->cutoff);
  return site_energy_neighbor(i);
}

/* ----------------------------------------------------------------------
   compute energy of site
   assume a neighbor list already exists
------------------------------------------------------------------------- */

double AppRelaxMB::site_energy_neighbor(int i)
{
  return pair->energy(i,numneigh,neighs,xyz,type);
}

/* ----------------------------------------------------------------------
   perform a particle event with Metropolis algorithm
   event = translation move of up to delta
------------------------------------------------------------------------- */

void AppRelaxMB::site_event_rejection(int i, RandomPark *random)
{
  double xold[3];
  double dx, dy, dz;
  double rc = pair->cutoff;

  // build i's neighbor list
  neighbor(i, rc);

  // save i's neighbor list
  int numneigh_i = numneigh;
  int *neighs_i = new int[numneigh_i];
  for (int jj = 0; jj < numneigh_i; jj++) neighs_i[jj] = neighs[jj];

  // initial energy — use site_energy to match diagnostic exactly
  double einitial = site_energy(i);
  for (int jj = 0; jj < numneigh_i; jj++) {
    int j = neighs_i[jj];
    if (j >= nlocal) continue;
    einitial += site_energy(j);
  }

  // save old position and bin
  xold[0] = xyz[i][0]; xold[1] = xyz[i][1]; xold[2] = xyz[i][2];
  int binold = bin[i];

  // generate random displacement
  double rsq = 1.0e20;
  while (rsq > deltasq) {
    dx = delta * (random->uniform() - 0.5);
    dy = delta * (random->uniform() - 0.5);
    if (dimension == 3) dz = delta * (random->uniform() - 0.5);
    else dz = 0.0;
    rsq = dx*dx + dy*dy + dz*dz;
  }

  // apply move and PBC
  xyz[i][0] += dx; xyz[i][1] += dy; xyz[i][2] += dz;
  if (xyz[i][0] < subxlo) xyz[i][0] += xprd;
  if (xyz[i][0] >= subxhi) xyz[i][0] -= xprd;
  if (xyz[i][1] < subylo) xyz[i][1] += yprd;
  if (xyz[i][1] >= subyhi) xyz[i][1] -= yprd;
  if (xyz[i][2] < subzlo) xyz[i][2] += zprd;
  if (xyz[i][2] >= subzhi) xyz[i][2] -= zprd;

  // update bin before neighbor search
  bin[i] = site2bin(i);

  // final energy — use site_energy to match diagnostic exactly
  double efinal = site_energy(i);
  for (int jj = 0; jj < numneigh_i; jj++) {
    int j = neighs_i[jj];
    if (j >= nlocal) continue;
    efinal += site_energy(j);
  }

  delete[] neighs_i;

  // accept or reject
  int success = 0;
  if (efinal <= einitial) {
    success = 1;
  } else if (temperature == 0.0) {
    xyz[i][0] = xold[0]; xyz[i][1] = xold[1]; xyz[i][2] = xold[2];
    bin[i] = binold;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    xyz[i][0] = xold[0]; xyz[i][1] = xold[1]; xyz[i][2] = xold[2];
    bin[i] = binold;
  } else success = 1;

  if (success) {
    bin[i] = binold;
    move(i);
    naccept++;
  }
}



// void AppRelaxMB::site_event_rejection(int i, RandomPark *random)
// {
//   double xold[3];
//   double dx,dy,dz;

//   double rc = pair->cutoff;
//   neighbor(i,rc+delta);

//   // neighbor list of i before the move
//   int numneigh_i = numneigh;
//   int *neighs_i = new int[numneigh_i];
//   for (int jj = 0; jj < numneigh_i; jj++) neighs_i[jj] = neighs[jj];

//   //  initial energy for atom i and all its neighbors 
//   double einitial = site_energy_neighbor(i);
//   for (int jj = 0; jj < numneigh_i; jj++){
//     int j = neighs_i[jj];
//     if (j >= nlocal) continue; // skip ghost atoms
//     neighbor(j, rc+delta);
//     // einitial += site_energy_neighbor(j);
//     einitial += pair->energy_neighbor_contribution(j, numneigh, neighs, xyz, type, i);
//   }

//   xold[0] = xyz[i][0];
//   xold[1] = xyz[i][1];
//   xold[2] = xyz[i][2];

//   double rsq = 1.0e20;
//   while (rsq > deltasq) {
//     dx = delta * (random->uniform() - 0.5);
//     dy = delta * (random->uniform() - 0.5);
//     if (dimension == 3) dz = delta * (random->uniform() - 0.5);
//     else dz = 0.0;
//     rsq = dx*dx + dy*dy + dz*dz;
//   }

//   xyz[i][0] += dx;
//   xyz[i][1] += dy;
//   xyz[i][2] += dz;

//   if (i == 0) {
//     printf("move: dx=%.6f dy=%.6f dz=%.6f  |dr|=%.6f\n",
//            dx, dy, dz, sqrt(dx*dx+dy*dy+dz*dz));
//     fflush(stdout);
//   }

//   // handle periodic bounary conditions
//   if (xyz[i][0] < subxlo) xyz[i][0] += xprd;
//   if (xyz[i][0] >= subxhi) xyz[i][0] -= xprd;
//   if (xyz[i][1] < subylo) xyz[i][1] += yprd;
//   if (xyz[i][1] >= subyhi) xyz[i][1] -= yprd;
//   if (xyz[i][2] < subzlo) xyz[i][2] += zprd;
//   if (xyz[i][2] >= subzhi) xyz[i][2] -= zprd;

//   // final energy of i and all its neighbors

//   neighbor(i, rc+delta);
//   double efinal = site_energy_neighbor(i);
//   for (int jj = 0; jj < numneigh_i; jj++){
//     int j = neighs_i[jj];
//     if (j >= nlocal) continue; // skip ghost atoms
//     neighbor(j, rc+delta);
//     efinal += pair->energy_neighbor_contribution(j, numneigh, neighs, xyz, type, i);
//   }

//   delete[] neighs_i;

//   // accept or reject via Boltzmann criterion

//   // app_relax_MB.cpp - temporary debug
//   if (i == 0) {
//       printf("dE = %f, einitial = %f, efinal = %f\n", 
//             efinal-einitial, einitial, efinal);
//       fflush(stdout);
//   }

//   int success = 0;

//   if (efinal <= einitial) {
//     if (i == 0) printf("  -> accepted by downhill branch\n");
//     success = 1;
//   } else if (temperature == 0.0) {
//     xyz[i][0] = xold[0];
//     xyz[i][1] = xold[1];
//     xyz[i][2] = xold[2];
//   } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
//     if (i == 0) printf("  -> rejected by Boltzmann\n");
//     xyz[i][0] = xold[0];
//     xyz[i][1] = xold[1];
//     xyz[i][2] = xold[2];
//   } else {
//     success = 1;
//     printf("UPHILL ACCEPTED: atom=%d dE=%f einitial=%f efinal=%f\n",
//            i, efinal-einitial, einitial, efinal);
//   }

//   if (success) {
//     move(i);
//     naccept++;
//   }
// }
