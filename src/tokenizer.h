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
   Code based off that tokenizer.h file in LAMMPS made by
   Contributing authors: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#ifndef SPK_TOKENIZER_H
#define SPK_TOKENIZER_H

#include "spktype.h"

#include <exception>
#include <string>
#include <vector>

namespace SPPARKS_NS {

#define TOKENIZER_DEFAULT_SEPARATORS " \t\r\n\f"

class Tokenizer {
    std::string text;
    std::string separators;
    size_t start;
    size_t ntokens;

    public:
        Tokenizer(std::string str, std::string separators = TOKENIZER_DEFAULT_SEPARATORS);
        Tokenizer(Tokenizer &&);
        Tokenizer(const Tokenizer &);
        Tokenizer &operator=(const Tokenizer &);
        Tokenizer &operator=(Tokenizer &&);
        void swap(Tokenizer &);

        void reset();
        void skip(int n = 1);
        bool has_next() const;
        bool contains(const std::string &str) const;
        std::string next();

        size_t count();
        std::vector<std::string> as_vector();
};

/** General Tokenizer exception class */

class TokenizerException : public std::exception {
    std::string message;

    public:
        // remove unused default constructor
        TokenizerException() = delete;

        /** Thrown during retrieving or skipping tokens
         * 
         *  \param msg String with error message
         *  \param token String of the token/work that casued the error */
        explicit TokenizerException(const std::string &msg, const std::string &token);

        /** Retrieve message describing the thrown exception
         * \return string with error message  */
        const char *what() const noexcept override { return message.c_str(); }
};

/** Exception thrown by ValueTokenizer with trying to convert an invalid integer string */

class InvalidIntegerException : public TokenizerException {

    public:
    /** Thrown during converting string to integer number
     * 
     * \param token String of the token/word that caused the error */
    explicit InvalidIntegerException(const std::string &token) :
    TokenizerException("Not a valid integer number", token)
    {
    }
};

/** Exception thrown by ValueTokenizer when trying to convert an floating point string */

class InvalidFloatException : public TokenizerException {
    public:
    /** Thrown during converting string to floating point number
     * 
     * \public token String of the token/word that caused the error */
    explicit InvalidFloatException(const std::string &token) :
        TokenizerException("Not a valid floating-point number", token)
    {
    }
};

class ValueTokenizer {
    Tokenizer tokens;

    public:
        ValueTokenizer(const std::string &str,
                       const std::string &separators = TOKENIZER_DEFAULT_SEPARATORS);
        ValueTokenizer(const ValueTokenizer &) = default;
        ValueTokenizer(ValueTokenizer &&);
        ValueTokenizer &operator=(const ValueTokenizer &);
        ValueTokenizer &operator=(ValueTokenizer &&);
        void swap(ValueTokenizer &);

        std::string next_string();
        tagint next_tagint();
        bigint next_bigint();
        int next_int();
        double next_double();

        bool has_next() const;
        bool contains(const std::string &value) const;
        void skip(int ntokens = 1);

        size_t count();
};

} // namespace SPPARKS_NS

#endif