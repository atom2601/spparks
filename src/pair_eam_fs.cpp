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

/* ----------------------------------------------------------------------
   Contributing author: Anthony Tom (University of Tennessee - Knoxville)
   Code based off that pair_eam_fs.cpp file in LAMMPS made by
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "app_off_lattice.h"
#include "comm_off_lattice.h"
#include "pair_eam_fs.h"
#include "potential.h"
#include "potential_file_reader.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

#include <cstring>
#include <iostream>

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

PairEAMFS::PairEAMFS(SPPARKS *spk) : PairEAM(spk)
{
    MPI_Comm_rank(world,&me);
    one_coeff = 1;
    manybody_flag = 1;
    he_flag = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
   read EAM Finnis-Sinclair file
------------------------------------------------------------------------- */

void PairEAMFS::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  // if (narg != 3 + potential->pair->ntypes) {
  //   std::cout << potential->pair->ntypes <<" , "<< narg << std::endl;
  //   error->all(FLERR, "Number of element to type mappings does not match number of atom types"); }

  // read EAM Finnis-Sinclair file

  if (fs) {
    for (i = 0; i < fs->nelements; i++) delete[] fs->elements[i];
    delete[] fs->elements;
    memory->destroy(fs->mass);
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    delete fs;
  }
  fs = new Fs();
  read_file(arg[2]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if "NULL"

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i], "NULL") == 0 ) {
        map[i - 2] = -1;
        continue;
    }
    for (j = 0; j < fs->nelements; j ++)
      if (strcmp(arg[i], fs->elements[j]) == 0) break;
    if (j < fs->nelements)
      map[i - 2] = j;
    else
      error->all(FLERR, "No matching element in EAM potential file");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = potential->pair->ntypes;
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++) setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
        if (map[i] >= 0 && map[j] >= 0) {
            setflag[i][j] = 1;
            // if (i ==j) atom->set_mass(FLERR, i, fs->mass[map[i]]); // not using mass in spparks
            count++;
        }
    scale[i][j] = 1.0;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficenits");
}

/* ----------------------------------------------------------------------
   read a multi-element DYNAMO setfl file
------------------------------------------------------------------------- */

void PairEAMFS::read_file(char *filename)
{
  Fs *file = fs;

  // read potential file
  if (me == 0) {
    PotentialFileReader reader(spk, filename, he_flag ? "eam/he" : "eam/fs", unit_convert_flag);

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY, unit_convert);
    try {
        reader.skip_line();
        reader.skip_line();
        reader.skip_line();

        // extract element names from nelements line
        ValueTokenizer values = reader.next_values(1);
        file->nelements = values.next_int();

        if ((int) values.count() != file->nelements + 1)
          error->one(FLERR, "Incorrect element names in EAM potential file");
        
        file->elements = new char *[file->nelements];
        for (int i = 0; i < file->nelements; i++)
          file->elements[i] = utils::strdup(values.next_string());
        
        if (he_flag)
          values = reader.next_values(6);
        else
          values = reader.next_values(5);
        file->nrho = values.next_int();
        file->drho = values.next_double();
        file->nr = values.next_int();
        file->dr = values.next_double();
        file->cut = values.next_double();
        if (he_flag) rhomax = values.next_double();

        if ((file->nrho <= 0) || (file->nr <= 0) || (file->dr <= 0.0))
          error->one(FLERR, "Invalid EAM potential file");

        memory->create(file->mass, file->nelements, "pair:mass");
        memory->create(file->frho, file->nelements, file->nrho + 1, "pair:frho");
        memory->create(file->rhor, file->nelements, file->nelements, file->nr + 1, "pair:rhor");
        memory->create(file->z2r, file->nelements, file->nelements, file->nr + 1, "pair:z2r");

        for (int i = 0; i < file->nelements; i++) {
            values = reader.next_values(2);
            values.next_int(); // ignore
            file->mass[i] = values.next_double();
            // std::cout << file->mass[i] << std::endl;

            reader.next_dvector(&file->frho[i][1], file->nrho);
            if (unit_convert) {
              for (int j = 1; j <= file->nrho; ++j) file->frho[i][j]; // *= conversion_factor;
            }

            for (int j = 0; j < file->nelements; j++) {
              reader.next_dvector(&file->rhor[i][j][1], file->nr);
            }
        }

        for (int i = 0; i < file->nelements; i++) {
          for (int j = 0; j <= i; j++) {
            reader.next_dvector(&file->z2r[i][j][1], file->nr);
            if (unit_convert) {
              for (int k = 1; k <= file->nr; ++k) file->z2r[i][j][k] ;//*= conversion_factor;
            }
          }
        }
    } catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }

  // broadcast potential information

  // MPI STUFF ??????????

  // allocate memory on other procs
  if (me != 0) {
    file->elements = new char *[file->nelements];
    for (int i =0; i < file->nelements; i++) file->elements[i] = nullptr;
    memory->create(file->mass, file->nelements, "pair:mass");
    memory->create(file->frho, file->nelements, file->nrho + 1, "pair:frho");
    memory->create(file->rhor, file->nelements, file->nelements, file->nr + 1, "pair:rhor");
    memory->create(file->z2r, file->nelements, file->nelements, file->nr + 1, "pair:z2r");
  }

  // broadcase file->elements string array
//     for (int i = 0; i < file->nelements; i++) {
//     int n;
//     if (comm->me == 0) n = strlen(file->elements[i]) + 1;
//     M
//   }
}

/* ----------------------------------------------------------------------
   copy read-in setfl potential to standard array format
------------------------------------------------------------------------- */

void PairEAMFS::file2array()
{
  int i, j, m, n;
  int ntypes = potential->pair->ntypes;

  // set function params directly from fs file

  nrho = fs->nrho;
  nr = fs->nr;
  drho = fs->drho;
  dr = fs->dr;
  if (he_flag)
    rhomin = rhomax - (nrho - 1) * drho;
  else
    rhomax = (nrho - 1) * drho;

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of fs elements + 1 for zero array

  nfrho = fs->nelements + 1;
  memory->destroy(frho);
  memory->create(frho, nfrho, nrho + 1, "pair:frho");

  // copy each elements frho to global frho

  for (i = 0; i < fs->nelements; i++)
    for (m = 1; m <= nrho; m++) frho[i][m] = fs->frho[i][m];

  // add extra frho of zeroes for non-EAM types to point (pair hybrid)
  // this is necessary b/c fp is still cpmuted for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho - 1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to 
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0)
      type2frho[i] = map[i];
    else
      type2frho[i] = nfrho -1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = square of # of fs elements

  nrhor = fs->nelements * fs->nelements;
  memory->destroy(rhor);
  memory->create(rhor, nrhor, nr + 1, "pair:rhor");

  // copy each element pair rhor to global rhor

  n = 0;
  for (i = 0; i < fs->nelements; i++)
    for (j = 0; j < fs->nelements; j++) {
      for (m = 1; m <= nr; m++) rhor[n][m] = fs->rhor[i][j][m];
      n++;
    }

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to 
  // for fs files, there is a full NxN set of rhor arrays
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j ++) type2rhor[i][j] = map[i] * fs->nelements + map[j];
  
  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of fs elements

  nz2r = fs->nelements * (fs->nelements + 1) /2;
  memory->destroy(z2r);
  memory->create(z2r, nz2r, nr + 1, "pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < fs->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = fs->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow, icol;
  for (i = 1; i <= ntypes; i++) {
    for (j=1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}