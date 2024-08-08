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
   Code based off the platform.h file in LAMMPS 
   Only copies over code that was necessary, tried to keep everything in
   the same order as that path file
------------------------------------------------------------------------- */

#ifndef SPK_PLATFORM_H
#define SPK_PLATFORM_H

/*! \file platform.h */

#include "spktype.h"
#include <cstdio>
#include <string>
#include <vector>

namespace SPPARKS_NS {
namespace platform{

  /*! Platform specific file path component separator
   *
   * This is a string with the character that separates directories and filename in paths on
   * a platform. If multiple are characters are provided, the first is the preferred one. */

#if !defined(_WIN32)
   constexpr char filepathsep[] = "/";
#else
   constexpr char filepathsep[] = "\\/";
#endif

     /*! Platform specific path environment variable component separator
   *
   * This is the character that separates entries in "PATH"-style environment variables. */

#if !defined(_WIN32)
   constexpr char pathvarsep = ':';
#else
   constexpr char pathvarsep = ';';
#endif

  /*! Get list of entries in a path environment variable
   *
   * This provides a list of strings of the entries in an environment
   * variable that is containing a "path" like "PATH" or "LD_LIBRARY_PATH".
   *
   * \param   var  name of the environment variable
   * \return  vector with strings of all entries in that path variable */

  std::vector<std::string> list_pathenv(const std::string &var);

  /*! Strip off leading part of path, return just the filename
   *
   * \param path file path
   * \return file name */

  std::string path_basename(const std::string &path);

  /*! Join two pathname segments
   *
   * This uses the forward slash '/' character unless LAMMPS is compiled
   * for Windows where it uses the backward slash '\\'
   *
   * \param   a  first path
   * \param   b  second path
   * \return     combined path */

  std::string path_join(const std::string &a, const std::string &b);

   /*! Check if file exists and is readable
   *
   * \param path file path
   * \return true if file exists and is readable */

  bool file_is_readable(const std::string &path); 

} // namespace platform
} // namespace SPPARKS_NS

#endif