/* fget.h -- alternatives to fscanf that die on error. */
/* Last edited on 2023-01-28 01:48:36 by stolfi */

#ifndef fget_H
#define fget_H

/* This module provides alternatives to {fscanf} that will skip
  spaces, then parse one value of a given type, aborting the program
  in case of syntax error or unexpected end-of-file.
  
  Created on Dec/2002 by J. Stolfi, from Modula-3 version ca. 1995,
  inspired on Lex.m3 by L. Cardelli. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

void fget_skip_spaces(FILE *f); 
  /* Skips spaces (SPACE, NBSP, TAB, NUL) until the first non-space
    character or end-of-file.  Will NOT skip line-breaks or page-breaks. */

void fget_skip_formatting_chars(FILE *f); 
  /* Skips all blank `formatting' characters --- namely spaces (SPACE, NBSP,
    TAB, NUL), line-breaks (CR, LF), and page-breaks (FF, VT) ---
    until the first non-blank character, or end-of-file. */

void fget_match(FILE *f, char *t);
  /* Requires the string {t} to be the next thing on {f}, and
    consumes that string. */

bool_t fget_test_char(FILE *f, int c);
  /* Checks whether the next character is {c} (which may be anything, including
    space, line-break, page-break, or even EOF). If it is, consumes that
    character and returns TRUE. If it is something else (even if it is a
    space, line-break or page-break, or end-of-file), returns FALSE and
    leaves the character there.  */

void fget_skip_to_eol(FILE *f);
  /* Reads arbitrary characters from {f} until an end-of-line or end-of-file.
    The end-of-line character is consumed. */

/* The following procedures will skip spaces (SPACEs, NBSPs, TABs and NULs)
  before the desired input.  They will NOT skip line-breaks or page-breaks. */

char fget_char(FILE *f);
  /* Skips spaces, then reads a single character --- which must
    not be a space, line break, or page break. Fails
    if a line break, page break, or end-of-file occurs before the character. */

bool_t fget_bool(FILE *f);
  /* The procedure {fget_bool_t} will skip spaces, then read a single
    character: 'T','t','Y','y','1' for TRUE, 'F','f','N','n','0' for FALSE. Anything
    else is an error.  Note that, if the file contains "TRUE" or
    "FALSE", only the first letter will be consumed. */
    
char *fget_to_delim(FILE *f, char *dels);
  /* Skips spaces, then reads zero or more characters, until end-of-file
    or the first formatting character (space, line break, or page
    break), or one of the characters in {dels}, which is not consumed.
    If {dels} is {NULL}, it is ignored. The result is returned as a
    newly allocated, zero-terminated string (which may be empty). */

char *fget_string(FILE *f);
  /* Skips spaces. Fails if hits end_of_file, line break, or page break.
    Otherwise reads one or more characters, until end-of-file or the
    first formatting character (space, line break, or page break), which
    is not consumed. The result is returned as a newly allocated,
    zero-terminated string (which has at least one character). */
    
int fget_int(FILE *f);
int32_t fget_int32(FILE *f);
int64_t fget_int64(FILE *f);
  /* Skips spaces, reads the longest possible non-empty string that
    looks like a decimal integer (with optional '+' or '-' sign), and
    tries to convert it to a signed integer value of the specified
    size. A failure in any of these steps is a fatal error. */
    
unsigned int fget_uint(FILE *f, int base);
uint32_t fget_uint32(FILE *f, int base);
uint64_t fget_uint64(FILE *f, int base);
  /* Skips spaces, reads the longest possible non-empty string that
    looks like an unsigned integer, and tries to convert it to an 
    unsigned integer value of the specified size. A failure in any 
    of these steps is a fatal error.
    
    The digit string found is interpreted relative to the given
    {base}, which should be in {2..36}. If {base} is 10 or less, the
    function consumes any string of digits in '0'..'9'; if {base} is
    11 or more, it consumes any string of digits or letters 'a'..'z'
    (or 'A'..'Z') which represent values in 10..35. In either case,
    the procedure fails if any digit has a value {base} or more. */

uint64_t fget_bin64(FILE *f);
  /* Skips spaces, reads the longest possible non-empty string of
    digits, ans tries to parseit as the binary representation of an
    {uint64_t} value. Fails if there are no digits, if it runs into
    a digit greater than 1, or if the value is greater than {2^64-1}. */

double fget_double(FILE *f);
  /* Skips spaces, reads the longest possible non-empty string that
    looks like a floating-point number (including optional exponent),
    and tries to convert it to a {double} value. A failure in any of
    these steps is a fatal error. */

void fget_eol(FILE *f);
  /* Skips any spaces and requires the next non-space character to be a 
    newline, which it consumes.  It is a fatal error if the next non-space is not a newline,
    or end-of-file is found.  Equivalent to
    {fget_skip_spaces(f); fget_match(f, "\n")}. */

void fget_comment_or_eol(FILE *f, int cmtc);
  /* Similar to {fget_eol}, but allows arbitrary comment string,
    preceded by spaces and the character {cmtc}, before the newline.
    Namely, first does {fget_skip_spaces(rd)}. Then, if the next character is newline,
    consumes it; if it is {cmtc}, consumes it and any others up to and including the newline. 
    
    Otherwise it fails.  The character {cmtc} should not be a space or newline. */

bool_t fget_test_comment_or_eol(FILE *f, int cmtc);
  /* First does {fget_skip_spaces(rd)}. Then, if the next character
    is {cmtc}, consumes it and other characters up to the newline and returns {true}; if it is newline,
    consumes it and returns {TRUE}; otherwise, leaves that character there 
    and returns {FALSE}. The character {cmtc} should not be a space or newline. */

void fget_skip_and_match(FILE *f, char *t);
  /* Equivalent to {fget_skip_spaces(f); fget_match(f, t)}. */

bool_t fget_skip_and_test_char(FILE *f, int c);
  /* Equivalent to {fget_skip_spaces(f); fget_test_char(f, c)}. */

#endif
