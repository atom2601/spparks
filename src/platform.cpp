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
   Code based off the platform.cpp file in LAMMPS 
   Only copies over code that was necessary, tried to keep everything in
   the same order as that path file
------------------------------------------------------------------------- */

#include "platform.h"

#include "text_file_reader.h"
#include "utils.h"

#include <fmt/format.h>
#include <deque>
#include <exception>
// #include <mpi.h>

using namespace SPPARKS_NS;

/* ----------------------------------------------------------------------
   split a "path" environment variable into a list
------------------------------------------------------------------------- */
std::vector<std::string> platform::list_pathenv(const std::string &var)
{
   std::vector<std::string> dirs;
   const char *ptr = getenv(var.c_str());
   if (ptr == nullptr) return dirs;

   std::string pathvar = ptr;
   std::size_t first = 0, next;
   while (true) {
      next = pathvar.find_first_of(pathvarsep, first);
      if (next == std::string::npos) {
         dirs.push_back(pathvar.substr(first));
         break;
      } else {
         dirs.push_back(pathvar.substr(first, next - first));
         first = next + 1;
      }
   }
   return dirs;
}

/* ----------------------------------------------------------------------
   strip off leading part of path, return just the filename
------------------------------------------------------------------------- */
std::string platform::path_basename(const std::string &path)
{
   size_t start = path.find_last_of(platform::filepathsep);

   if (start == std::string::npos) {
      start = 0;
   } else {
      start += 1;
   }

   return path.substr(start);
}

/* ----------------------------------------------------------------------
   join two paths.
   if one of the two is an empty string just return the other unmodified
   if the first string ends in the separator or the second begins with one, trim them
------------------------------------------------------------------------- */
std::string platform::path_join(const std::string &a, const std::string &b)
{
   if (a.empty()) return b;
   if (b.empty()) return a;

   // remove trailing separator(s) in first part
   std::string joined = a;
   while (joined.find_last_of(platform::filepathsep) == joined.size() - 1) {
      for (const auto &s : platform::filepathsep)
         if (joined.back() == s) joined.pop_back();
   }

   // skip over leading separator(s) in second part
   std::size_t skip = 0;
   while (b.find_first_of(platform::filepathsep, skip) == skip) ++ skip;

   // combine and return
   joined += platform::filepathsep[0] + b.substr(skip);
   return joined;
}

/* ----------------------------------------------------------------------
   try to open file for reading to prove if it exists and is accessible
------------------------------------------------------------------------- */

bool platform::file_is_readable(const std::string &path)
{
   FILE *fp = fopen(path.c_str(), "r");
   if (fp) {
      fclose(fp);
      return true;
   }
   return false;
}