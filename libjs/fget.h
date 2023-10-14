/* fget.h -- alternatives to fscanf that die on error. */
/* Last edited on 2023-10-13 13:48:50 by stolfi */

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

bool_t fget_test_eof(FILE *f);
  /* Returns {TRUE} iff there are no more characters to be read from {f}.
    Otherwise leaves {f} effectively unchanged.  */

bool_t fget_test_char(FILE *f, char c);
  /* If {f} is exhausted, returns {FALSE}.
    Otherwise, checks whether the next character is {c} (which may be anything,
    including space or formatting char). If it is, consumes that
    character and returns TRUE.  If it is some other character, (even if it is a
    space or formatting char), returns FALSE and leaves the character
    there. 
    
    Note that {c = -1 = '\377'} will match only the byte 0xff, not {EOF}. */

/* The following procedures will skip spaces (SPACEs, NBSPs, TABs and NULs)
  before the desired input.  They will NOT skip line-breaks or page-breaks. */

char fget_char(FILE *f);
  /* Skips spaces, then reads a single character --- which must
    not be a space or formatting char (such as newline). Fails
    if a formatting char or end-of-file occurs before the character.
    
    Note that a result of {-1 = '\377'} means the byte 0xff
    not {EOF}. */

bool_t fget_bool(FILE *f);
  /* The procedure {fget_bool_t} will skip spaces, then read a single
    character: 'T','t','Y','y','1' for TRUE, 'F','f','N','n','0' for FALSE. Anything
    else is an error.  Note that, if the file contains "TRUE" or
    "FALSE", only the first letter will be consumed. */
    
char *fget_to_delim(FILE *f, char del);
  /* Skips spaces, then takes zero or more characters, until end-of-file
    or the first formatting character (space, line break, or page
    break), or the character {del}, which is not consumed.
    
    If {del} is '\n', it just takes characters until a formatting char
    or end-of-file. The result is returned as a newly allocated,
    zero-terminated string (which may be empty). */

char *fget_to_delims(FILE *f, char *dels);
  /* Similar to {fget_to_delim}, but any of the characters in
    {dels} will act as a delimiter, and will not be consumed.
    If {dels} is {NULL}, empty, or "\n", it just takes characters
    unit a formatting char or end-of-file. */

char *fget_line(FILE *f);
  /* Collects all characters from the current point up to but 
    not including the end of the line or end of file into as a newly allocated,
    zero-terminated string, which is returned.  The end-of-file,
    if any, is not consumed. */

char *fget_string(FILE *f);
  /* Skips spaces. Fails if hits end_of_file, line break, or page break.
    Otherwise reads one or more characters, until end-of-file or the
    a formatting character (space, line break, or page break), which
    is not consumed.  The result is returned as a newly allocated,
    zero-terminated string (which has at least one character). */
    
int32_t fget_int32(FILE *f);
int64_t fget_int64(FILE *f);
  /* Skips spaces, reads the longest possible non-empty string that
    looks like a decimal integer (with optional '+' or '-' sign), and
    tries to convert it to a signed integer value of the specified
    size. A failure in any of these steps is a fatal error. */
    
uint32_t fget_uint32(FILE *f, uint32_t base);
uint64_t fget_uint64(FILE *f, uint32_t base);
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
    or end-of-file is found.  Equivalent to  {fget_skip_spaces(f); fget_match(f, "\n")}. */

void fget_skip_to_eol(FILE *f);
  /* Reads arbitrary characters from {f} until an end-of-line or end-of-file.
    The end-of-line character is consumed. */

void fget_comment_or_eol(FILE *f, char cmtc);
  /* Similar to {fget_eol}, but allows arbitrary comment string,
    preceded by spaces and the character {cmtc}, before the newline.
    
    Namely, first does {fget_skip_spaces(f)}. Then, if the next
    character is newline, consumes it; if it is {cmtc}, consumes it and
    any other characters up to and including the newline.
    
    Fails in all other cases, including if EOF occurs before the newline
    while skipping the comment. The character {cmtc} had better not be a
    space or {EOF}. If {cmtc} is '\n', the effect is the same as {fget_eol}. */

bool_t fget_test_comment_or_eol(FILE *f, char cmtc);
  /* Similar to {fget_comment_or_eol}, but returns {FALSE}
    if it finds anything (including end-of-file) before a newline
    or {cmtc}. In this case, that character is not consumed.
    Otherwise the procedure returns {TRUE}, after consuming all characteds
    up to and including the newline.
    
    The character {cmtc} should not be a space or {EOF}. */

void fget_skip_spaces_and_match(FILE *f, char *t);
  /* Equivalent to {fget_skip_spaces(f); fget_match(f, t)}. */

bool_t fget_skip_and_test_char(FILE *f, char c);
  /* Equivalent to {fget_skip_spaces(f); fget_test_char(f, c)}. */

/* HANDY TOOLS */

bool_t fget_is_formatting_char(char c);
  /* Returns TRUE iff {c} is a formatting character, namely a space
    (SPACE, TAB, NUL, NBSP) line break (CR, LF) or page break (FF, VT).
    Note that {c} cannot be {EOF}, since {(char)EOF} coincides with
    '\240', non-breaking space. */
    
bool_t fget_is_space(char c);
  /* Returns TRUE iff {c} is a space char (SPACE, TAB, NUL, NBSP). Note
    that {c} cannot be {EOF}, since {(char)EOF} coincides with '\240',
    non-breaking space.. */

#endif
