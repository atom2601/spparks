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
#include "app_potts.h"
#include "comm_lattice.h"
#include "solve.h"
#include "random_park.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

#include <map>

using namespace SPPARKS_NS;

enum{SPIN,NEIGH,NEIGHONLY};

/* ---------------------------------------------------------------------- */

AppPotts::AppPotts(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  // parse arguments

  if (narg < 4) error->all("Illegal app_style command");

  nspins = atoi(arg[1]);
  if (strcmp(arg[2],"spin") == 0) rejectstyle = SPIN;
  else if (strcmp(arg[2],"neigh") == 0) rejectstyle = NEIGH;
  else if (strcmp(arg[2],"neighonly") == 0) rejectstyle = NEIGHONLY;
  else error->all("Illegal app_style command");
  int seed = atoi(arg[3]);
  RandomPark *random = new RandomPark(seed);

  if (rejectstyle == SPIN) dt_sweep = 1.0/nspins;
  else if (rejectstyle == NEIGH) dt_sweep = 1.0/maxneigh;
  else if (rejectstyle == NEIGHONLY) dt_sweep = 1.0/maxneigh;

  options(narg-4,&arg[4]);

  // define lattice and partition it across processors
  
  create_lattice();
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];

  // initialize my portion of lattice
  // each site = one of nspins
  // loop over global list so assignment is independent of # of procs
  // use map to see if I own global site

  if (infile) read_file();

  else {
    std::map<int,int> hash;
    for (int i = 0; i < nlocal; i++)
      hash.insert(std::pair<int,int> (id[i],i));
    std::map<int,int>::iterator loc;
    
    int isite;
    for (int iglobal = 1; iglobal <= nglobal; iglobal++) {
      isite = random->irandom(nspins);
      loc = hash.find(iglobal);
      if (loc != hash.end()) lattice[loc->second] = isite;
    }
  }

  delete random;
}

/* ---------------------------------------------------------------------- */

AppPotts::~AppPotts()
{
  delete [] sites;
  delete [] unique;
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppPotts::site_energy(int i)
{
  int isite = lattice[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != lattice[neighbor[i][j]]) eng++;
  return (double) eng;
}

/* ----------------------------------------------------------------------
   perform a site event with null bin rejection
------------------------------------------------------------------------- */

void AppPotts::site_event_rejection(int i, RandomPark *random)
{
  if (rejectstyle == SPIN) site_event_rejection_spins(i,random);
  else if (rejectstyle == NEIGH) site_event_rejection_neighbors(i,random);
  else site_event_rejection_neighbors_only(i,random);
}

/* ----------------------------------------------------------------------
   perform a site event with null bin rejection
   flip to random spin from 1 to nspins
------------------------------------------------------------------------- */

void AppPotts::site_event_rejection_spins(int i, RandomPark *random)
{
  int oldstate = lattice[i];
  double einitial = site_energy(i);

  // event = random spin from 1 to nspins, including self

  int iran = (int) (nspins*random->uniform()) + 1;
  if (iran > nspins) iran = nspins;
  lattice[i] = iran;
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i] = oldstate;
  }

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5*numneigh[i]) mask[i] = 1;
    if (lattice[i] != oldstate)
      for (int j = 0; j < numneigh[i]; j++)
	mask[neighbor[i][j]] = 0;
  }
}

/* ----------------------------------------------------------------------
   perform a site event with null bin rejection
   flip to random neighbor spin with null bin extending to maxneigh
------------------------------------------------------------------------- */

void AppPotts::site_event_rejection_neighbors(int i, RandomPark *random)
{
  int oldstate = lattice[i];
  double einitial = site_energy(i);

  // events = spin flips to neighboring site different than self

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == lattice[i]) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  int iran = (int) (maxneigh*random->uniform());
  if (iran >= nevent) return;
  lattice[i] = unique[iran];
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i] = oldstate;
  }

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5*numneigh[i]) mask[i] = 1;
    if (lattice[i] != oldstate)
      for (int j = 0; j < numneigh[i]; j++)
	mask[neighbor[i][j]] = 0;
  }
}

/* ----------------------------------------------------------------------
   perform a site event with no null bin rejection
   flip to random neighbor spin
   technically this is an incorrect rejection-KMC algorithm
------------------------------------------------------------------------- */

void AppPotts::site_event_rejection_neighbors_only(int i, RandomPark *random)
{
  int oldstate = lattice[i];
  double einitial = site_energy(i);

  // events = spin flips to neighboring site different than self

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == lattice[i]) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  if (nevent == 0) return;
  int iran = (int) (nevent*random->uniform());
  if (iran >= nevent) iran = nevent-1;
  lattice[i] = unique[iran];
  double efinal = site_energy(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    lattice[i] = oldstate;
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    lattice[i] = oldstate;
  }

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5*numneigh[i]) mask[i] = 1;
    if (lattice[i] != oldstate)
      for (int j = 0; j < numneigh[i]; j++)
	mask[neighbor[i][j]] = 0;
  }
}

/* ----------------------------------------------------------------------
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppPotts::site_propensity(int i)
{
  // events = spin flips to neighboring site different than self
  // disallows wild flips = flips to value different than all neighs

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == lattice[i]) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  // for each flip:
  // compute energy difference between initial and final state
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor

  int oldstate = lattice[i];
  double einitial = site_energy(i);
  double efinal;
  double prob = 0.0;

  for (m = 0; m < nevent; m++) {
    lattice[i] = unique[m];
    efinal = site_energy(i);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
  }

  lattice[i] = oldstate;
  return prob;
}

/* ----------------------------------------------------------------------
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppPotts::site_event(int i, RandomPark *random)
{
  int j,m,value;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double efinal;

  int oldstate = lattice[i];
  double einitial = site_energy(i);
  double prob = 0.0;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = lattice[neighbor[i][j]];
    if (value == oldstate) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;

    lattice[i] = value;
    efinal = site_energy(i);
    if (efinal <= einitial) prob += 1.0;
    else if (temperature > 0.0) prob += exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) break;
  }

  // compute propensity changes for self and neighbor sites
  // ignore update of neighbor sites with isite < 0

  int nsites = 0;
  int isite = i2site[i];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(i);

  for (j = 0; j < numneigh[i]; j++) {
    m = neighbor[i][j];
    isite = i2site[m];
    if (isite < 0) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
  }

  solve->update(nsites,sites,propensity);
}
