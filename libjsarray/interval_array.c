/* See interval_array.h */
/* Last edited on 2010-04-23 10:29:56 by stolfi */ 

#include <interval_array.h>

#include <array.h>
#include <indexing.h>
#include <indexing_descr.h>

#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <affirm.h>
#include <bool.h>
#include <interval.h>
// #include <jsmath.h>
// 
// #include <limits.h>
// #include <assert.h>
// #include <string.h>
// #include <math.h>
#include <stdio.h>
#include <stdlib.h>

array_typeimpl(interval_array_t, interval_array, interval_t);

void interval_array_elem_write(FILE *wr, interval_t *valP) 
  { fprintf(wr, "[ %+23.15e _ %+23.15e ]", valP->end[0], valP->end[1]); }
  
void interval_array_elem_read(FILE *rd, interval_t *valP)
  { 
    fget_skip_spaces(rd);
    fget_match(rd, "[");
    valP->end[0] = fget_double(rd);
    fget_skip_spaces(rd);
    fget_match(rd, "_");
    valP->end[1] = fget_double(rd);
    fget_skip_spaces(rd);
    fget_match(rd, "]");
  }

array_io_impl(interval_array_t, interval_array, interval_t);

/* #define interval_array_elem_add(X,Y) ((X)+(Y)) */
/* #define interval_array_elem_mul(X,Y) ((X)*(Y)) */
/* #define interval_array_elem_zero (0.0) */
/* #define interval_array_elem_one (1.0) */

/* array_linalg_impl(interval_array_t, interval_array, interval_t); */
