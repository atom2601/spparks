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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair.h"
#include "error.h"
#include "utils.h"

using namespace SPPARKS_NS;

enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};
enum{R,RSQ,BMP};

/* ---------------------------------------------------------------------- */

Pair::Pair(SPPARKS *spk) : 
  Pointers(spk), map(nullptr)
{
  mix_flag = GEOMETRIC;
  allocated = 0;

  one_coeff = 0;
  manybody_flag = 0;
  restartinfo = 1;
  unit_convert_flag = utils::NOCONVERT;

  // comm_forward = comm_reverse = 0;
}

// Pair::~Pair()
// {
//   delete[] map;
// }

/* ---------------------------------------------------------------------- */

void Pair::init()
{
  int i,j;

  if (!allocated) error->all(FLERR,"All pair coeffs are not set");

  // I,I coeffs must be set
  // init_one() will check if I,J is set explicitly or inferred by mixing

  for (i = 1; i <= ntypes; i++)
    if (setflag[i][i] == 0) error->all(FLERR,"All pair coeffs are not set");

  // style-specific initialization

  init_style();

  // call init_one() for each I,J
  // set cutsq for each I,J, used to neighbor
  // cutforce = max of all I,J cutoffs

  double cut;
  cutoff = 0.0;
  for (i = 1; i <= ntypes; i++)
    for (j = i; j <= ntypes; j++) {
      cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
      cutoff = MAX(cutoff,cut);
    }
}

/* ----------------------------------------------------------------------
   mixing of pair potential prefactors (epsilon)
------------------------------------------------------------------------- */

double Pair::mix_energy(double eps1, double eps2, double sig1, double sig2)
{
  double value;
  if (mix_flag == GEOMETRIC)
    value = sqrt(eps1*eps2);
  else if (mix_flag == ARITHMETIC)
    value = sqrt(eps1*eps2);
  else if (mix_flag == SIXTHPOWER)
    value = 2.0 * sqrt(eps1*eps2) *
      pow(sig1,3.0) * pow(sig2,3.0) / (pow(sig1,6.0) * pow(sig2,6.0));
  return value;
}

/* ----------------------------------------------------------------------
   mixing of pair potential distances (sigma, cutoff)
------------------------------------------------------------------------- */

double Pair::mix_distance(double sig1, double sig2)
{
  double value;
  if (mix_flag == GEOMETRIC)
    value = sqrt(sig1*sig2);
  else if (mix_flag == ARITHMETIC)
    value = 0.5 * (sig1+sig2);
  else if (mix_flag == SIXTHPOWER)
    value = pow((0.5 * (pow(sig1,6.0) + pow(sig2,6.0))),1.0/6.0);
  return value;
}

/* -------------------------------------------------------------------
   build element to atom type mapping for manybody potentials
   also clear and reset setflag[][] array and check missing entries
---------------------------------------------------------------------- */

void Pair::map_element2type(int narg, char **arg, bool update_setflag)
{
  int i,j;
  // const int ntypes = atom->ntypes;

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if "NULL"
  // nelements = # of unique elements
  // elements = list of element names

  if (narg != ntypes)
    error->all(FLERR, "Number of element to type mappings does not match number of atom types");

  if (elements) {
    for (i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;
  }
  elements = new char*[ntypes];
  for (i = 0; i < ntypes; i++) elements[i] = nullptr;

  nelements = 0;
  map[0] = -1;
  for (i = 1; i <= narg; i++) {
    std::string entry = arg[i-1];
    if (entry == "NULL") {
      map[i] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (entry == elements[j]) break;
    map[i] = j;
    if (j == nelements) {
      elements[j] = utils::strdup(entry);
      nelements++;
    }
  }
  
  // if requested, clear setflag[i][j] and set it for type pairs
  // where both are mapped to elements in map.

  if (update_setflag) {

    int count = 0;
    for (i = 1; i <= ntypes; i++) {
      for (j = i; j <= ntypes; j++) {
        setflag[i][j] = 0;
        if ((map[i] >= 0) && (map[j] >= 0)) {
          setflag[i][j] = 1;
          count++;
        }
      }
    }

    if (count == 0) error->all(FLERR, "incorrect args for pair coefficients");
  }
}

// misc, half-compete functions from lammps

// void Pair::ev_setup(int eflag, int vflag, int alloc)
// {
//   int i,n;

//   eflag_either = eflag;
//   eflag_global = eflag & ENERGY_GLOBAL;
//   eflag_atom = eflag & ENERGY_ATOM;

//   vflag_global = vflag & VIRTUAL_PAIR;

//   evflag = eflag_either || vflag_either;



// }

// void Pair::ev_unset()
// {
//   evflag = 0;

//   eflag_either = 0;
//   eflag_global = 0;
//   eflag_atom = 0;

//   vflag_either = 0;
//   vflag_global = 0;
//   vflag_atom = 0;
//   cvflag_atom = 0;
//   vflag_fdotr = 0;
// }