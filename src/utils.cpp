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

#include "app_off_lattice.h"
#include "comm_off_lattice.h"
#include "text_file_reader.h"
#include "platform.h"
#include "utils.h"
#include "error.h"

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
    // auto me = comm->get_me();

    std::string filepath = get_potential_file_path(name);

    if (!filepath.empty()) {
        // std::string unit_style = spk->update->unit_style;
        std::string date = get_potential_date(filepath, "potential");
        std::string units = get_potential_units(filepath, "potential");

        // if (!date.empty() && (me == 0))
        //     logmesg(spk, "Reading potential file {} with DATE: {}\n", name, date);

        // if (auto_convert == nullptr) {
        //     if (!units.empty() && (units != unit_style) && (me == 0)) {
        //         std::string formattedMesg = fmt::format("Potential file {} requires {} units but {} are in use", 
        //                                                 name, units, unit_style);
        //         const char* cstr = formattedMesg.c_str();
        //         error->one(FLERR, cstr);
        //         return nullptr;
        //     }
        // } else {
        //     if (units.empty() || units == unit_style) {
        //         *auto_convert = NOCONVERT;
        //     } else {
        //         if ((units == "metal") && (unit_style == "real") && (*auto_convert & METAL2REAL)) {
        //             *auto_convert = METAL2REAL;
        //         } else if ((units == "real") && (unit_style == "metal") && (*auto_convert & REAL2METAL)) {
        //             *auto_convert = REAL2METAL;
        //         } else {
        //             std::string formattedMesg = fmt::format("Potential file {} requires {} units but {} units are in use",
        //                                                     name, units, unit_style);
        //             const char* cstr = formattedMesg.c_str();
        //             error->one(FLERR, cstr);
        //             return nullptr;
        //         }
        //     }
        //     if ((*auto_convert != NOCONVERT) && (me == 0)) {
        //         std::string formattedMesg = fmt::format("Converting potential file in {} units to {} units",
        //                                                 units, unit_style);
        //         const char* cstr = formattedMesg.c_str();
        //         error->warning(FLERR, cstr);
        //     }
        // }
        return fopen(filepath.c_str(), "r");
    }
    return nullptr;
}

extern "C" {

/* Typedef'd pointer to get abstract datatype. */
typedef struct regex_t *re_t;
typedef struct regex_context_t *re_ctx_t;

/* Compile regex string pattern to a regex_t-array. */
static re_t re_compile(re_ctx_t context, const char *pattern);

/* Find matches of the compiled pattern inside text. */
static int re_matchp(const char *text, re_t pattern, int *matchlen);

/* Definitions: */

#define MAX_REGEXP_OBJECTS 256 /* Max number of regex symbols in expression. */
#define MAX_CHAR_CLASS_LEN 256 /* Max length of character-class buffer in.   */

enum {
  RX_UNUSED,
  RX_DOT,
  RX_BEGIN,
  RX_END,
  RX_QUESTIONMARK,
  RX_STAR,
  RX_PLUS,
  RX_CHAR,
  RX_CHAR_CLASS,
  RX_INV_CHAR_CLASS,
  RX_DIGIT,
  RX_NOT_DIGIT,
  RX_INTEGER,
  RX_NOT_INTEGER,
  RX_FLOAT,
  RX_NOT_FLOAT,
  RX_ALPHA,
  RX_NOT_ALPHA,
  RX_WHITESPACE,
  RX_NOT_WHITESPACE /*, BRANCH */
};

typedef struct regex_t {
  unsigned char type; /* CHAR, STAR, etc.                      */
  union {
    unsigned char ch;   /*      the character itself             */
    unsigned char *ccl; /*  OR  a pointer to characters in class */
  } u;
} regex_t;

typedef struct regex_context_t {
  /* MAX_REGEXP_OBJECTS is the max number of symbols in the expression.
       MAX_CHAR_CLASS_LEN determines the size of buffer for chars in all char-classes in the expression. */
  regex_t re_compiled[MAX_REGEXP_OBJECTS];
  unsigned char ccl_buf[MAX_CHAR_CLASS_LEN];
} regex_context_t;

int re_match(const char *text, const char *pattern)
{
  regex_context_t context;
  int dummy;
  return re_matchp(text, re_compile(&context, pattern), &dummy);
}

int re_find(const char *text, const char *pattern, int *matchlen)
{
  regex_context_t context;
  return re_matchp(text, re_compile(&context, pattern), matchlen);
}

/* Private function declarations: */
static int matchpattern(regex_t *pattern, const char *text, int *matchlen);
static int matchcharclass(char c, const char *str);
static int matchstar(regex_t p, regex_t *pattern, const char *text, int *matchlen);
static int matchplus(regex_t p, regex_t *pattern, const char *text, int *matchlen);
static int matchone(regex_t p, char c);
static int matchdigit(char c);
static int matchint(char c);
static int matchfloat(char c);
static int matchalpha(char c);
static int matchwhitespace(char c);
static int matchmetachar(char c, const char *str);
static int matchrange(char c, const char *str);
static int matchdot(char c);
static int ismetachar(char c);

/* Semi-public functions: */
int re_matchp(const char *text, re_t pattern, int *matchlen)
{
  *matchlen = 0;
  if (pattern != nullptr) {
    if (pattern[0].type == RX_BEGIN) {
      return ((matchpattern(&pattern[1], text, matchlen)) ? 0 : -1);
    } else {
      int idx = -1;

      do {
        idx += 1;

        if (matchpattern(pattern, text, matchlen)) {
          if (text[0] == '\0') return -1;

          return idx;
        }
      } while (*text++ != '\0');
    }
  }
  return -1;
}

re_t re_compile(re_ctx_t context, const char *pattern)
{
  regex_t *const re_compiled = context->re_compiled;
  unsigned char *const ccl_buf = context->ccl_buf;
  int ccl_bufidx = 1;

  char c;    /* current char in pattern   */
  int i = 0; /* index into pattern        */
  int j = 0; /* index into re_compiled    */

  while (pattern[i] != '\0' && (j + 1 < MAX_REGEXP_OBJECTS)) {
    c = pattern[i];

    switch (c) {
        /* Meta-characters: */
      case '^': {
        re_compiled[j].type = RX_BEGIN;
      } break;
      case '$': {
        re_compiled[j].type = RX_END;
      } break;
      case '.': {
        re_compiled[j].type = RX_DOT;
      } break;
      case '*': {
        re_compiled[j].type = RX_STAR;
      } break;
      case '+': {
        re_compiled[j].type = RX_PLUS;
      } break;
      case '?': {
        re_compiled[j].type = RX_QUESTIONMARK;
      } break;

        /* Escaped character-classes (\s \w ...): */
      case '\\': {
        if (pattern[i + 1] != '\0') {
          /* Skip the escape-char '\\' */
          i += 1;
          /* ... and check the next */
          switch (pattern[i]) {
              /* Meta-character: */
            case 'd': {
              re_compiled[j].type = RX_DIGIT;
            } break;
            case 'D': {
              re_compiled[j].type = RX_NOT_DIGIT;
            } break;
            case 'i': {
              re_compiled[j].type = RX_INTEGER;
            } break;
            case 'I': {
              re_compiled[j].type = RX_NOT_INTEGER;
            } break;
            case 'f': {
              re_compiled[j].type = RX_FLOAT;
            } break;
            case 'F': {
              re_compiled[j].type = RX_NOT_FLOAT;
            } break;
            case 'w': {
              re_compiled[j].type = RX_ALPHA;
            } break;
            case 'W': {
              re_compiled[j].type = RX_NOT_ALPHA;
            } break;
            case 's': {
              re_compiled[j].type = RX_WHITESPACE;
            } break;
            case 'S': {
              re_compiled[j].type = RX_NOT_WHITESPACE;
            } break;

              /* Escaped character, e.g. '.' or '$' */
            default: {
              re_compiled[j].type = RX_CHAR;
              re_compiled[j].u.ch = pattern[i];
            } break;
          }
        }
        /* '\\' as last char in pattern -> invalid regular expression. */
      } break;

        /* Character class: */
      case '[': {
        /* Remember where the char-buffer starts. */
        int buf_begin = ccl_bufidx;

        /* Look-ahead to determine if negated */
        if (pattern[i + 1] == '^') {
          re_compiled[j].type = RX_INV_CHAR_CLASS;
          i += 1;                  /* Increment i to avoid including '^' in the char-buffer */
          if (pattern[i + 1] == 0) /* incomplete pattern, missing non-zero char after '^' */
          {
            return nullptr;
          }
        } else {
          re_compiled[j].type = RX_CHAR_CLASS;
        }

        /* Copy characters inside [..] to buffer */
        while ((pattern[++i] != ']') && (pattern[i] != '\0')) {
          /* Missing ] */
          if (pattern[i] == '\\') {
            if (ccl_bufidx >= MAX_CHAR_CLASS_LEN - 1) { return nullptr; }
            if (pattern[i + 1] == 0) /* incomplete pattern, missing non-zero char after '\\' */
            {
              return nullptr;
            }
            ccl_buf[ccl_bufidx++] = pattern[i++];
          } else if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
            return nullptr;
          }
          ccl_buf[ccl_bufidx++] = pattern[i];
        }
        if (ccl_bufidx >= MAX_CHAR_CLASS_LEN) {
          /* Catches cases such as [00000000000000000000000000000000000000][ */
          return nullptr;
        }
        /* Null-terminate string end */
        ccl_buf[ccl_bufidx++] = 0;
        re_compiled[j].u.ccl = &ccl_buf[buf_begin];
      } break;

        /* Other characters: */
      default: {
        re_compiled[j].type = RX_CHAR;
        re_compiled[j].u.ch = c;
      } break;
    }
    /* no buffer-out-of-bounds access on invalid patterns -
     * see https://github.com/kokke/tiny-regex-c/commit/1a279e04014b70b0695fba559a7c05d55e6ee90b */
    if (pattern[i] == 0) { return nullptr; }

    i += 1;
    j += 1;
  }
  /* 'RX_UNUSED' is a sentinel used to indicate end-of-pattern */
  re_compiled[j].type = RX_UNUSED;

  return (re_t) re_compiled;
}

/* Private functions: */
static int matchdigit(char c)
{
  return isdigit(c);
}

static int matchint(char c)
{
  return (matchdigit(c) || (c == '-') || (c == '+'));
}

static int matchfloat(char c)
{
  return (matchint(c) || (c == '.') || (c == 'e') || (c == 'E'));
}

static int matchalpha(char c)
{
  return isalpha(c);
}

static int matchwhitespace(char c)
{
  return isspace(c);
}

static int matchalphanum(char c)
{
  return ((c == '_') || matchalpha(c) || matchdigit(c));
}

static int matchrange(char c, const char *str)
{
  return ((c != '-') && (str[0] != '\0') && (str[0] != '-') && (str[1] == '-') &&
          (str[1] != '\0') && (str[2] != '\0') && ((c >= str[0]) && (c <= str[2])));
}

static int matchdot(char c)
{
#if defined(RE_DOT_MATCHES_NEWLINE) && (RE_DOT_MATCHES_NEWLINE == 1)
  (void) c;
  return 1;
#else
  return c != '\n' && c != '\r';
#endif
}

static int ismetachar(char c)
{
  return ((c == 's') || (c == 'S') || (c == 'w') || (c == 'W') || (c == 'd') || (c == 'D'));
}

static int matchmetachar(char c, const char *str)
{
  switch (str[0]) {
    case 'd':
      return matchdigit(c);
    case 'D':
      return !matchdigit(c);
    case 'i':
      return matchint(c);
    case 'I':
      return !matchint(c);
    case 'f':
      return matchfloat(c);
    case 'F':
      return !matchfloat(c);
    case 'w':
      return matchalphanum(c);
    case 'W':
      return !matchalphanum(c);
    case 's':
      return matchwhitespace(c);
    case 'S':
      return !matchwhitespace(c);
    default:
      return (c == str[0]);
  }
}

static int matchcharclass(char c, const char *str)
{
  do {
    if (matchrange(c, str)) {
      return 1;
    } else if (str[0] == '\\') {
      /* Escape-char: increment str-ptr and match on next char */
      str += 1;
      if (matchmetachar(c, str)) {
        return 1;
      } else if ((c == str[0]) && !ismetachar(c)) {
        return 1;
      }
    } else if (c == str[0]) {
      if (c == '-') {
        return ((str[-1] == '\0') || (str[1] == '\0'));
      } else {
        return 1;
      }
    }
  } while (*str++ != '\0');

  return 0;
}

static int matchone(regex_t p, char c)
{
  switch (p.type) {
    case RX_DOT:
      return matchdot(c);
    case RX_CHAR_CLASS:
      return matchcharclass(c, (const char *) p.u.ccl);
    case RX_INV_CHAR_CLASS:
      return !matchcharclass(c, (const char *) p.u.ccl);
    case RX_DIGIT:
      return matchdigit(c);
    case RX_NOT_DIGIT:
      return !matchdigit(c);
    case RX_INTEGER:
      return matchint(c);
    case RX_NOT_INTEGER:
      return !matchint(c);
    case RX_FLOAT:
      return matchfloat(c);
    case RX_NOT_FLOAT:
      return !matchfloat(c);
    case RX_ALPHA:
      return matchalphanum(c);
    case RX_NOT_ALPHA:
      return !matchalphanum(c);
    case RX_WHITESPACE:
      return matchwhitespace(c);
    case RX_NOT_WHITESPACE:
      return !matchwhitespace(c);
    default:
      return (p.u.ch == c);
  }
}

static int matchstar(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  int prelen = *matchlen;
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text >= prepos) {
    if (matchpattern(pattern, text--, matchlen)) return 1;
    (*matchlen)--;
  }

  *matchlen = prelen;
  return 0;
}

static int matchplus(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  const char *prepos = text;
  while ((text[0] != '\0') && matchone(p, *text)) {
    text++;
    (*matchlen)++;
  }
  while (text > prepos) {
    if (matchpattern(pattern, text--, matchlen)) return 1;
    (*matchlen)--;
  }
  return 0;
}

static int matchquestion(regex_t p, regex_t *pattern, const char *text, int *matchlen)
{
  if (p.type == RX_UNUSED) return 1;
  if (matchpattern(pattern, text, matchlen)) return 1;
  if (*text && matchone(p, *text++)) {
    if (matchpattern(pattern, text, matchlen)) {
      (*matchlen)++;
      return 1;
    }
  }
  return 0;
}

/* Iterative matching */
static int matchpattern(regex_t *pattern, const char *text, int *matchlen)
{
  int pre = *matchlen;
  do {
    if ((pattern[0].type == RX_UNUSED) || (pattern[1].type == RX_QUESTIONMARK)) {
      return matchquestion(pattern[0], &pattern[2], text, matchlen);
    } else if (pattern[1].type == RX_STAR) {
      return matchstar(pattern[0], &pattern[2], text, matchlen);
    } else if (pattern[1].type == RX_PLUS) {
      return matchplus(pattern[0], &pattern[2], text, matchlen);
    } else if ((pattern[0].type == RX_END) && pattern[1].type == RX_UNUSED) {
      return (text[0] == '\0');
    }
    (*matchlen)++;
  } while ((text[0] != '\0') && matchone(*pattern++, *text++));

  *matchlen = pre;
  return 0;
}
}