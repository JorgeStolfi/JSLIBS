/* fget_geo.h -- extends fget.h for geometric args. */
/* Last edited on 2008-07-14 22:32:39 by stolfi */

#ifndef fget_geo_H
#define fget_geo_H

/* Copyright © 2008 Jorge Stolfi, Unicamp. See note at end of file. */

/* This interface provides convenient tools for parsing command
  line arguments whose values are real vectors. */

#include <stdio.h>

#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <r6.h>

r2_t fget_r2(FILE *rd);
r3_t fget_r3(FILE *rd);
r4_t fget_r4(FILE *rd);
r6_t fget_r6(FILE *rd);
  /* Reads from {rd} (with {fget_double}) the next {N} real numbers as
    the coordinates of a point; where {N} is 2,3,4, or 6. Does not
    skip over line breaks, so all numbers must be on the current line
    of {rd}. */

void fget_rn(FILE *rd, double p[], int n);
  /* Reads from {rd} (with {fget_double}) the next {n} real numbers,
    and stores them in {p[0.n-1]}. Does not skip over line breaks, so
    all numbers must be on the current line {rd}. */

r3_t fget_r3_dir(FILE *rd);
  /* Same as {fget_r3} but normalizes the result to unit length. */

/* Copyright © 2008 by Jorge Stolfi.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appears in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty of any kind.
*/

#endif
