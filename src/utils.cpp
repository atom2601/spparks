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
   Code based off the utils.cpp file in LAMMPS 
   Only copies over code that was necessary, tried to keep everything in
   the same order as that utils file
------------------------------------------------------------------------- */

#include "utils.h"

#include "error.h"
#include "text_file_reader.h"
#include "platform.h"
#include "format.h"
#include "comm_off_lattice.h"

#include <cctype>
#include <cerrno>
#include <cstring>
#include <ctime>

/*! \file utils.cpp*/

extern "C" {
/** Match text against a (simplified) regular expression
 * (regexp will be complied automatically). */
static int re_match(const char *text, const char *pattern);
}

using namespace SPPARKS_NS;

bool utils::strmatch(const std::string &text, const std::string &pattern)
{
    const int pos = re_match(text.c_str(), pattern.c_str());
    return (pos >= 0);
}

void utils::logmesg(SPPARKS *spk, const std::string &mesg)
{
    if (spk->screen) fputs(mesg.c_str(), spk->screen);
    if (spk->logfile) fputs(mesg.c_str(), spk->logfile);
}

void utils::fmtargs_logmesg(SPPARKS *spk, fmt::string_view format, fmt::format_args args)
{
    try {
        logmesg(spk, fmt::vformat(format, args));
    } catch (fmt::format_error &e) {
        logmesg(spk, std::string(e.what()) + "\n");
    }
}

std::string utils::getsyserror()
{
    return {strerror(errno)};
}

/* ----------------------------------------------------------------------
   Make copy of string in new storage. Works like the (non-portable)
   C-style strdup() but also accepts a C++ string as argument.
------------------------------------------------------------------------- */
char *utils::strdup(const std::string &text)
{
    auto tmp = new char[text.size() + 1];
    strcpy(tmp, text.c_str()); // NOLINT
    return tmp;
}

/* ----------------------------------------------------------------------
   Replace UTF-8 encoded chars with known ASCII equivalents
------------------------------------------------------------------------- */
std::string utils::utf8_subst(const std::string &line)
{
  const auto *const in = (const unsigned char *) line.c_str();
  const int len = line.size();
  std::string out;

  for (int i = 0; i < len; ++i) {

    // UTF-8 2-byte character
    if ((in[i] & 0xe0U) == 0xc0U) {
      if ((i + 1) < len) {
        // NON-BREAKING SPACE (U+00A0)
        if ((in[i] == 0xc2U) && (in[i + 1] == 0xa0U)) out += ' ', ++i;
        // MODIFIER LETTER PLUS SIGN (U+02D6)
        if ((in[i] == 0xcbU) && (in[i + 1] == 0x96U)) out += '+', ++i;
        // MODIFIER LETTER MINUS SIGN (U+02D7)
        if ((in[i] == 0xcbU) && (in[i + 1] == 0x97U)) out += '-', ++i;
      }
      // UTF-8 3-byte character
    } else if ((in[i] & 0xf0U) == 0xe0U) {
      if ((i + 2) < len) {
        // EN QUAD (U+2000)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x80U)) out += ' ', i += 2;
        // EM QUAD (U+2001)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x81U)) out += ' ', i += 2;
        // EN SPACE (U+2002)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x82U)) out += ' ', i += 2;
        // EM SPACE (U+2003)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x83U)) out += ' ', i += 2;
        // THREE-PER-EM SPACE (U+2004)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x84U)) out += ' ', i += 2;
        // FOUR-PER-EM SPACE (U+2005)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x85U)) out += ' ', i += 2;
        // SIX-PER-EM SPACE (U+2006)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x86U)) out += ' ', i += 2;
        // FIGURE SPACE (U+2007)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x87U)) out += ' ', i += 2;
        // PUNCTUATION SPACE (U+2008)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x88U)) out += ' ', i += 2;
        // THIN SPACE (U+2009)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x89U)) out += ' ', i += 2;
        // HAIR SPACE (U+200A)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x8aU)) out += ' ', i += 2;
        // ZERO WIDTH SPACE (U+200B)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x8bU)) out += ' ', i += 2;
        // LEFT SINGLE QUOTATION MARK (U+2018)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x98U)) out += '\'', i += 2;
        // RIGHT SINGLE QUOTATION MARK (U+2019)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x99U)) out += '\'', i += 2;
        // LEFT DOUBLE QUOTATION MARK (U+201C)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x9cU)) out += '"', i += 2;
        // RIGHT DOUBLE QUOTATION MARK (U+201D)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0x9dU)) out += '"', i += 2;
        // NARROW NO-BREAK SPACE (U+202F)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x80U) && (in[i + 2] == 0xafU)) out += ' ', i += 2;
        // WORD JOINER (U+2060)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa0U)) out += ' ', i += 2;
        // INVISIBLE SEPARATOR (U+2063)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa3U)) out += ' ', i += 2;
        // INVISIBLE PLUS (U+2064)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x81U) && (in[i + 2] == 0xa4U)) out += '+', i += 2;
        // MINUS SIGN (U+2212)
        if ((in[i] == 0xe2U) && (in[i + 1] == 0x88U) && (in[i + 2] == 0x92U)) out += '-', i += 2;
        // ZERO WIDTH NO-BREAK SPACE (U+FEFF)
        if ((in[i] == 0xefU) && (in[i + 1] == 0xbbU) && (in[i + 2] == 0xbfU)) out += ' ', i += 2;
      }
      // UTF-8 4-byte character
    } else if ((in[i] & 0xf8U) == 0xf0U) {
      if ((i + 3) < len) { ; }
    } else
      out += in[i];
  }
  return out;
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */
size_t utils::count_words(const char *text)
{
    size_t count = 0;
    const char *buf = text;
    char c = *buf;

    while (c) {
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') {
            c = *++buf;
            continue;
        };

        ++count;
        c = *++buf;

        while(c) {
            if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f') { break; }
            c = *++buf;
        }
    }

    return count;
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */
size_t utils::count_words(const std::string &text)
{
    return utils::count_words(text.c_str());
}

/* ----------------------------------------------------------------------
   return number of words
------------------------------------------------------------------------- */
size_t utils::count_words(const std::string &text, const std::string &separators)
{
    size_t count = 0;
    size_t start = text.find_first_not_of(separators);

    while (start != std::string::npos) {
        size_t end = text.find_first_of(separators, start);
        ++count;

        if (end == std::string::npos) {
            return count;
        } else {
            start = text.find_first_not_of(separators, end + 1);
        }
    } 
    
    return count;
}

/* ----------------------------------------------------------------------
   Convert string into words on whitespace while handling single and
   double quotes.
------------------------------------------------------------------------- */
std::vector<std::string> utils::split_words(const std::string &text)
{
    std::vector<std::string> list;
    const char *buf = text.c_str();
    std::size_t beg = 0;
    std::size_t len = 0;
    std::size_t add = 0;
    char c = *buf;

    while (c) {
        // leading whitespace
        if (c == ' ' || c == '\t' || c == '\r' || c == '\f') {
            c = *++buf;
            ++beg;
            continue;
        };
        len = 0;
    
    // handle escaped/quoted text
    quoted:

        // handle single quote
        if (c == '\'') {
            ++beg;
            add = 1;
            c = *++buf;
            while (((c != '\'') && (c != '\0')) || ((c == '\\') && (buf[1] == '\''))) {
                if ((c == '\\') && (buf[1] == '\'')) {
                    ++buf;
                    ++len;
                }
                c = *++ buf;
                ++len;
            }
            if (c != '\'') ++len;
            c = *++buf;

            // handle triple double quotation marks
        } else if ((c == '"') && (buf[1] == '"') && (buf[2] == '"') && (buf[3] != '"')) {
            len = 3;
            add = 1;
            buf += 3;
            c = *buf;

            // handle double quote
        } else if (c == '"') {
            ++beg;
            add = 1;
            c = *++ buf;
            while (((c != '"') && (c != '\0')) || ((c == '\\') && (buf[1] == '"'))) {
                if ((c = '\\') && (buf[1] == '"')) {
                    ++buf;
                    ++len;
                }
                c = *++ buf;
                ++len;
            }
            if (c != '"') ++len;
            c = *++buf;
        }

        // unquoted
        while (true) {
            if ((c == '\'') || (c == '"')) goto quoted;
            // skip escaped quote
            if ((c == '\\') && ((buf[1] == '\'') || (buf[1] == '"'))) {
                ++buf;
                ++len;
                c = *++buf;
                ++len;
            }
            if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n') || (c == '\f') || (c == '\0')) {
                list.push_back(text.substr(beg, len));
                beg += len + add;
                break;
            }
            c = *++buf;
            ++len;
        }
    }
    return list;
}

/* ----------------------------------------------------------------------
   Return whether string is a valid integer number
------------------------------------------------------------------------- */
bool utils::is_integer(const std::string &str)
{
    if (str.empty()) return false;
    
    return strmatch(str, "^[+-]?\\d+$");
}

bool utils::is_double(const std::string &str)
{
    if (str.empty()) return false;

    return strmatch(str, "^[+-]\\d+\\.?\\d*$") ||
        strmatch(str, "^[+-]?\\d+\\.?\\d*[eE][+-]?\\d+$") || 
        strmatch(str, "^[+-]?\\d*\\?\\d+$") ||
        strmatch(str, "^[+-]?\\d*\\.?\\d+[eE][+-]?\\d+$");
}

/* ----------------------------------------------------------------------
   try to find potential file as specified by name
   search current directory and the LAMMPS_POTENTIALS directory if
   specified
   second portion will not work for SPPARKS due to no potential library
------------------------------------------------------------------------- */
std::string utils::get_potential_file_path(const std::string &path)
{
    if (platform::file_is_readable(path)) {
        return path;
    } else {
        for (const auto &dir : platform::list_pathenv("LAMMPS_POTENTIALS")) {
            auto pot = platform::path_basename(path);
            auto filepath = platform::path_join(dir,pot);
            if (platform::file_is_readable(filepath)) return filepath;
        }
    }
    return "";
}

/* ----------------------------------------------------------------------
   read first line of potential file
   if it has a DATE field, return the following word
------------------------------------------------------------------------- */
std::string utils::get_potential_date(const std::string &path, const std::string &potential_name)
{
    TextFileReader reader(path, potential_name);
    reader.ignore_comments = false;

    char *line = reader.next_line();
    if (line == nullptr) return "";
    Tokenizer words(line);
    while (words.has_next()) {
        if (words.next() == "DATE:") {
            if (words.has_next()) return words.next();
        }
    }
    return "";
}

/* ----------------------------------------------------------------------
   read first line of potential file
   if it has UNITS field, return following word
------------------------------------------------------------------------- */
std::string utils::get_potential_units(const std::string &path, const std::string &potential_name)
{
    TextFileReader reader(path, potential_name);
    reader.ignore_comments = false;

    char *line = reader.next_line();
    if (line == nullptr) return "";
    Tokenizer words(line);
    while (words.has_next()) {
        if (words.next() == "UNITS:") {
            if (words.has_next()) return words.next();
        }
    }
    return "";
}

/* ----------------------------------------------------------------------
   return bitmask of supported conversions for a given property
------------------------------------------------------------------------- */
int utils::get_supported_conversions(const int property)
{
    if (property == ENERGY)
        return METAL2REAL | REAL2METAL;
    else
        return NOCONVERT;
}

/* ----------------------------------------------------------------------
   return conversion factor for a given property and conversion setting
   return 0.0 if unknown.
------------------------------------------------------------------------- */
double utils::get_conversion_factor(const int property, const int conversion)
{
    if (property == ENERGY) {
        if (conversion == NOCONVERT) {
            return 1.0;
        } else if (conversion == METAL2REAL) {
            return 23.060549;
        } else if (conversion == REAL2METAL) {
            return 1.0 / 23.060549;
        }
    }
    return 0.0;
}

/* ----------------------------------------------------------------------
   open a potential file as specified by name
   if fails, search in dir specified by env variable LAMMPS_POTENTIALS
------------------------------------------------------------------------- */
FILE *utils::open_potential(const std::string &name, SPPARKS *spk, int *auto_convert)
{
    auto error = spk->error;
    auto me = spk->comm->get_me();

    std::string filepath = get_potential_file_path(name);

    if (!filepath.empty()) {
        std::string unit_style = spk->update->unit_style;
        std::string date = get_potential_date(filepath, "potential");
        std::string units = get_potential_units(filepath, "potential");

        if (!date.empty() && (me == 0))
            logmesg(spk, "Reading potential file {} with DATE: {}\n", name, date);

        if (auto_convert == nullptr) {
            if (!units.empty() && (units != unit_style) && (me == 0)) {
                std::string formattedMesg = fmt::format("Potential file {} requires {} units but {} are in use", 
                                                        name, units, unit_style);
                const char* cstr = formattedMesg.c_str();
                error->one(FLERR, cstr);
                return nullptr;
            }
        } else {
            if (units.empty() || units == unit_style) {
                *auto_convert = NOCONVERT;
            } else {
                if ((units == "metal") && (unit_style == "real") && (*auto_convert & METAL2REAL)) {
                    *auto_convert = METAL2REAL;
                } else if ((units == "real") && (unit_style == "metal") && (*auto_convert & REAL2METAL)) {
                    *auto_convert = REAL2METAL;
                } else {
                    std::string formattedMesg = fmt::format("Potential file {} requires {} units but {} units are in use",
                                                            name, units, unit_style);
                    const char* cstr = formattedMesg.c_str();
                    error->one(FLERR, cstr);
                    return nullptr;
                }
            }
            if ((*auto_convert != NOCONVERT) && (me == 0)) {
                std::string formattedMesg = fmt::format("Converting potential file in {} units to {} units",
                                                        units, unit_style);
                const char* cstr = formattedMesg.c_str();
                error->warning(FLERR, cstr);
            }
        }
        return fopen(filepath.c_str(), "r");
    }
    return nullptr;
}