/* See double_array.h */
/* Last edited on 2010-04-23 10:25:09 by stolfi */ 

#include <double_array.h>

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

array_typeimpl(double_array_t, double_array, double);

void double_array_elem_write(FILE *wr, double *valP)
  { fprintf(wr, "%+23.15e", *valP); }
  
void double_array_elem_read(FILE *rd, double *valP) 
  { (*valP) = fget_double(rd); }

array_io_impl(double_array_t, double_array, double);

/* #define double_array_elem_add(X,Y) ((X)+(Y)) */
/* #define double_array_elem_mul(X,Y) ((X)*(Y)) */
/* #define double_array_elem_zero (0.0) */
/* #define double_array_elem_one (1.0) */

/* array_linalg_impl(double_array_t, double_array, double); */
