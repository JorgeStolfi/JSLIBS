#ifndef i3_H
#define i3_H

/* i3.h --- operations on points and vectors of Z^3 */
/* Last edited on 2024-11-20 13:49:14 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <sign.h>

/* OVERFLOW

  Many of the operations below may overflow if the scalar arguments
  (including scale factors and vector coordinates) are too large in
  absolute value.
  
  The "Max {M}:" information in the comments below is the maximum {M} such
  that the ouput is guaranteed correct as long as all input scalars
  in the symmetric range {[-M..+M]}.
  
  The "Output range:" information is the smallest symmetric range that contains
  the output scalars when the input scalars and
  coordinates span the symmetric range {[-M..+M]}.
  
  Note that signed {int2_t} variables naturally have an *asymmetric* range
  {[-N..N-1]} where {N = 2^31}. */
  

typedef struct i3_t { int32_t c[3]; } i3_t;

void i3_zero (i3_t *r);
  /* Sets {r} to the zero vector. */
  
void i3_all (int32_t x, i3_t *r);
  /* Sets all coordinates of {r} to the value {x}. */

void i3_axis (uint32_t i, i3_t *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void i3_add (i3_t *a, i3_t *b, i3_t *r);
  /* Sets {r = a + b}. Input MMax: {2^30-1}. Output range: {-2*M..+2*M}. */

void i3_sub (i3_t *a, i3_t *b, i3_t *r);
  /* Sets {r = a - b}. Input MMax: {2^30-1}. Output range: {-2*M..+2*M}. */

void i3_neg (i3_t *a, i3_t *r);
  /* Sets {r} to {-a}. Input MMax: {2^31-1}. Output range: {-M..+M}. */

uint32_t i3_L_inf_norm (i3_t *a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate).
    Input MMax: {2^31-1}. Output range: {0..+M}. */

uint64_t i3_L_inf_dist (i3_t *a, i3_t *b);
  /* Returns the L-infinity distance between {a} and {b} (max absolute diff).
    Input MMax: {2^31-1}. Output range: {0..+2*M} */
  
uint64_t i3_norm_sqr (i3_t *a);
  /* Returns the square of the Euclidean norm of {a}.  
    Input MMax: {2^31-1}. Output range: {0..+3*M^2}. */

uint64_t i3_dist_sqr (i3_t *a, i3_t *b);
  /* Returns the square of the Euclidean distance between {a} and {b}.
    Input MMax: {2^31-1}. Output range: {0..+12*M^2}. */

int64_t i3_dot (i3_t *a, i3_t *b);
  /* Dot product of vectors {a} and {b}. 
    Input MMax: {2^31-1}. Output range: {-3*M^2 .. +3*M^2}. */

void i3_cross (i3_t *a, i3_t *b, i3_t *r);
  /* Sets {r} to the cross product of {a} and {b}. 
    Input MMax: {2^15-1 = 32767}. Output range: {-2*M^2 .. +2*M^2}.  */

int64_t i3_det (i3_t *a, i3_t *b, i3_t *c);
  /* Returns the determinant of the 3 x 3 matrix whose rows are {a,b,c}.
    Input MMax: {2^31-1}. Output range: {-4*M^3 .. +4*M^3}. */

bool_t i3_eq(i3_t *p, i3_t *q);
  /* True iff points {p} and {q} are identical. */

void i3_throw_cube (int32_t m, i3_t *r);
  /* Sets {r} to a uniformly random point of the 3-cube {[-|m| .. +|m|]^3}.
    Input MMax: {2^31-1}. Output range: {-|m| .. +|m|}. */

void i3_print (FILE *f, i3_t *a);
  /* Prints {a} on file {f}, with some default format. */

void i3_gen_print (FILE *f, i3_t *a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%d", "(", " ", and ")", respectively. */

/* DERIVED TYPES */

vec_typedef(i3_vec_t,i3_vec,i3_t);
  /* An {i3_vec_t} is a vector of {i3_t}s. */

typedef bool_t i3_pred_t(i3_t *a);
  /* Type of a function that returns a {bool_t} value from an {i3_t} value. */

typedef double i3_double_func_t(i3_t *a);
  /* Type of a function that returns a {double} value from an {i3_t} value. */

typedef i3_t i3_map_t(i3_t *a);
  /* Type of a function that returns a {double} value from an {i3_t} value. */

#endif
