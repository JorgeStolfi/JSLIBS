/* See float_array.h */
/* Last edited on 2016-04-01 01:39:18 by stolfilocal */ 

#include <float_array.h>

#include <array.h>
#include <ix.h>
#include <ix_descr.h>

#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <affirm.h>
#include <bool.h>
// #include <jsmath.h>
// 
// #include <limits.h>
// #include <assert.h>
// #include <string.h>
// #include <math.h>
#include <stdio.h>
#include <stdlib.h>

array_typeimpl(float_array_t, float_array, float);

void float_array_elem_write(FILE *wr, float *valP)
  { fprintf(wr, "%15.8e", *valP); }
  
void float_array_elem_read(FILE *rd, float *valP)
{ (*valP) = (float)fget_double(rd); }

array_io_impl(float_array_t, float_array, float);

/* #define float_array_elem_add(X,Y) ((X)+(Y)) */
/* #define float_array_elem_mul(X,Y) ((X)*(Y)) */
/* #define float_array_elem_zero (0.0) */
/* #define float_array_elem_one (1.0) */

/* array_linalg_impl(float_array_t, float_array, double); */
