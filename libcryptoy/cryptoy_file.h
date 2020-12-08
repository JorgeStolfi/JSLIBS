/* File-based toy cryptography tools. */
/* Last edited on 2013-10-31 01:45:18 by stolfilocal */

#ifndef cryptoy_file_H
#define cryptoy_file_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#define cryptoy_file_H_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

int64_t cryptoy_file_xor(FILE *a, FILE *b, FILE *r);
  /* Writes to file {r} the bitwise XOR of the contents of files {a} and {b}.
    The files must be properly opened by the caller. Returns {-1} if read or write
    errors were detected, or file {a} is longer than file {b}. 
    sOtherwise returns the number of bytes written to file {r}. */

int64_t cryptoy_file_mix(FILE *a, FILE *b, FILE *x, FILE *r, FILE *s);
  /* Writes to files {r} and {s} the contents of files {a} and {b},
    respectively, after swapping the bits that correspond to the '1'
    bits of file {x}. The files must be properly opened by the caller.
    Returns {-1} if read or write errors were detected, or if files {a}
    and {b} have different lengths, or if they are longer than file {x}.
    Otherwise returns the number of bytes written to each of {r} and
    {s}. */

#define cryptoy_file_H_rights \
  "This file can be used and modified for any purpose as long as the source" \
  " is made available with any compiled code and the copyright and this" \
  " 'rights' note are preserved in all copies and derived versions."

#define cryptoy_file_H_warranty \
  "This file is provided 'as is'. The author and his employer are" \
  " not liable for any damages that may result from its use."

#endif
