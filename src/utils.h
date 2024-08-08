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
   Code based off the utils.h file in LAMMPS 
   Only copies over code that was necessary, tried to keep everything in
   the same order as that utils file
------------------------------------------------------------------------- */

#ifndef SPK_UTILS_H
#define SPK_UTILS_H

/*! \file utils.h */

#include "format.h"
#include "spktype.h"

// #include <mpi.h>

#include <vector>

namespace SPPARKS_NS {

// forward declarations
class Error;
class SPPARKS;

namespace utils {

    class CommOffLattice *comm;

    /*! Match text against a simplified regex pattern
     * 
     * \param text the text to be matched against the pattern
     * \param pattern the search pattern, which may contain regexp markers
     * \return true if the pattern matches, false if not */

    bool strmatch(const std::string &text, const std::string &pattern);

    /* Internal function handling the arguement list for logmesg(). */

    void fmtargs_logmesg(SPPARKS *spk, fmt::string_view format, fmt::format_args args);

    /*! Send formatted message to screen and logfile, if available
     * 
     * This function simplifies the repetitive task of outputting some
     * message to both the screen and/or the log file. The template
     * wrapper with fmtlib format and arguement processing allows
     * this function to work similar to ``fmt::print()``.
     * 
     * \param spk pointer to SPPARKS class instance
     * \param format format string of message to be printed
     * \param args arguements to format string */

    template <typename... Args> void logmesg(SPPARKS *spk, const std::string &format, Args &&...args)
    {
        fmtargs_logmesg(spk, format, fmt::make_format_args(args...));
    }

    /*! \overload
     *
     * \param spk pointer to SPPARKS class instance
     * \param mesg string with message to be printed */

    void logmesg(SPPARKS *spk, const std::string &mesg);

    /*! Return a string representing the current system error status
     * 
     * This is a wrapper around calling strerror(errorno).
     * 
     *  \return error string*/

    std::string getsyserror();

    /*! Make C-style copy of string in new storage
    * 
    * This allocates a storage buffer and copies the C-style or 
    * C++ style string into it. The buffer is allocated with "new" 
    * and thus needs to be deallocated with "delete[]".
    * 
    * \param text string that should be copies
    * \return new buffer with copy of string*/
    char *strdup(const std::string &text);

/*! Check if a string will likely have UTF-8 encoded characters
   *
   * UTF-8 uses the 7-bit standard ASCII table for the first 127 characters and
   * all other characters are encoded as multiple bytes.  For the multi-byte
   * characters the first byte has either the highest two, three, or four bits
   * set followed by a zero bit and followed by one, two, or three more bytes,
   * respectively, where the highest bit is set and the second highest bit set
   * to 0.  The remaining bits combined are the character code, which is thus
   * limited to 21-bits.
   *
   * For the sake of efficiency this test only checks if a character in the string
   * has the highest bit set and thus is very likely an UTF-8 character.  It will
   * not be able to tell this this is a valid UTF-8 character or whether it is a
   * 2-byte, 3-byte, or 4-byte character.
   *
\verbatim embed:rst

*See also*
   :cpp:func:`utils::utf8_subst`

\endverbatim
   * \param line  string that should be checked
   * \return true if string contains UTF-8 encoded characters (bool) */

    inline bool has_utf8(const std::string &line)
    {
        for (auto c : line)
            if (c & 0x80U) return true;
        return false;
    }

  /*! Replace known UTF-8 characters with ASCII equivalents
   *
\verbatim embed:rst

*See also*
   :cpp:func:`utils::has_utf8`

\endverbatim
   * \param line  string that should be converted
   * \return new string with ascii replacements (string) */

    std::string utf8_subst(const std::string &line);

    /*! Count words in string with custon choice of seperating characters
     * 
     * \param text string that should be searched
     * \param separators string containing characters that will be treated as whitespace
     * \return number of words found */
    
    size_t count_words(const std::string &text, const std::string &separators);

    /*! Count words in string, ignore any whitespace matching " \t\r\nf"
     * 
     * \param text string that should be searched
     * \return number of words found */

    size_t count_words(const std::string &text);

    /*! Count words in C-string, ignore any whitespace matching " \t\r\n\f"
    * 
    * \param text string that should be searched
    * \return number of words found */

   size_t count_words(const char *text);

  /*! Take text and split into non-whitespace words.
   *
   * This can handle strings with single and double quotes, escaped quotes,
   * and escaped codes within quotes, but due to using an STL container and
   * STL strings is rather slow because of making copies. Designed for
   * parsing command lines and similar text and not for time critical
   * processing.  Use a tokenizer class if performance matters.
   *
\verbatim embed:rst

*See also*
   :cpp:class:`Tokenizer`, :cpp:class:`ValueTokenizer`

\endverbatim
   * \param text string that should be split
   * \return STL vector with the words */

    std::vector<std::string> split_words(const std::string &text);

    /*! Check if string can be converted to valid integer
     * 
     * \param str string that should be checked
     * \return true, if string contains valid a integer, false otherwise */

    bool is_integer(const std::string &str);

    /*! Check if string can be converted to valid floating-point number
     * 
     * \param str string that should be checked
     * \return true, if string contains valid number, false otherwise */

    bool is_double(const std::string &str);

    /*! Determine full path of potential file. If file is not found in current directory,
     * search directories listed in LAMMPS_POTENTIALS enviornment variable
     * Above text is taken from the LAMMPS utils file, this directory does not currently
     * exist in SPPARKS
     * 
     * \param path file path
     * \return full path to potential file */

    std::string get_potential_file_path(const std::string &path);

    /*! Read potential file and return DATE field if it is present
     * 
     * \param path file path
     * \param potential_name name of potential that is being read
     * \return DATE field if present */

    std::string get_potential_date(const std::string &path, const std::string &potential_name);

    /*! Read potential file and return UNITS field if it is present
     * 
     * \param path file path
     * \param potential_name name of potential that is being used
     * \return UNITS field if present */

    std::string get_potential_units(const std::string &path, const std::string &potential_name);

    enum { NOCONVERT = 0, METAL2REAL = 1, REAL2METAL = 1 << 1 };
    enum { UNKNOWN = 0, ENERGY };

    /*! Return bitmask of available conversion factors for a given property
     *
     * \param property property to be converted
     * \return bitmask indicating available conversions */

    int get_supported_conversions(const int property);

    /*! Return unit conversion factor for given property and selected from/to units
     *
     * \param property property to be converted
     * \param conversion constant indicating the conversion
     * \return conversion factor */

    double get_conversion_factor(const int property, const int conversion);

    /*! Open a potential file as specified by *name*
    *
    * If opening the file directly fails, the function will search for
    * it in the list of folder pointed to by the environment variable
    * ``LAMMPS_POTENTIALS`` (if it is set).
    *
    * If the potential file has a ``UNITS`` tag in the first line, the
    * tag's value is compared to the current unit style setting.
    * The behavior of the function then depends on the value of the
    * *auto_convert* parameter.  If it is a null pointer, then the unit
    * values must match or else the open will fail with an error.  Otherwise
    * the bitmask that *auto_convert* points to is used check for
    * compatibility with possible automatic conversions by the calling
    * function.  If compatible, the bitmask is set to the required
    * conversion or ``utils::NOCONVERT``.
    *
    * \param name          file- or pathname of the potential file
    * \param lmp           pointer to top-level LAMMPS class instance
    * \param auto_convert  pointer to unit conversion bitmask or ``nullptr``
    * \return              FILE pointer of the opened potential file or ``nullptr`` */

    FILE *open_potential(const std::string &name, SPPARKS *spk, int *auto_convert);

} // namespace utils
} // namespace SPPARKS_NS

#endif