/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Anthony Tom (University of Tennessee - Knoxville)
   Code based off that potential_file_reader.cpp file in LAMMPS made by
   Contributing authors: SRichard Berger (Temple U)
------------------------------------------------------------------------- */

#ifndef SPK_POTENTIAL_FILE_READER_H
#define SPK_POTENTIAL_FILE_READER_H

#include "pointers.h"
#include "tokenizer.h"

namespace SPPARKS_NS {
class TextFileReader;

class PotentialFileReader : protected Pointers {
 protected:
  TextFileReader *reader;
  std::string filename;
  std::string filetype;
  int unit_convert;

  TextFileReader *open_potential(const std::string &path);

 public:
  PotentialFileReader(class SPPARKS *spk, const std::string &filename,
                      const std::string &potential_name, const int auto_convert = 0);
  PotentialFileReader(class SPPARKS *spk, const std::string &filename,
                      const std::string &potential_name, const std::string &name_suffix,
                      const int auto_convert = 0);
  ~PotentialFileReader() override;
  
  class CommOffLattice *comm;

  void ignore_comments(bool value);

  void set_bufsize(int);
  void rewind();
  void skip_line();
  char *next_line(int nparms = 0);
  void next_dvector(double *list, int n);
  ValueTokenizer next_values(int nparms,
                             const std::string &seperators = TOKENIZER_DEFAULT_SEPARATORS);
  
  // convenience functions
  double next_double();
  int next_int();
  tagint next_tagint();
  bigint next_bigint();
  std::string next_string();

  // unit conversion info
  int get_unit_convert() const { return unit_convert; }
};

} // namespace SPPARKS_NS

#endif