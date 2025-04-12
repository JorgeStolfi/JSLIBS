#ifndef i2_H
#define i2_H

/* Operations on points and vectors of Z^2 */
/* Last edited on 2025-03-11 20:05:08 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <sign.h>

/* OVERFLOW

  Many of the operations below may overflow if the scalar arguments
  (including scale factors and vector coordinates) are too large in
  absolute value.
  
  The "Input MMax:" information in the comments below is the maximum {M} such
  that the ouput is guaranteed correct as long as all input scalars
  are in the symmetric range {[-M..+M]}.
  
  The "Output range:" information is the smallest symmetric range that contains
  the output scalars when the input scalars and
  coordinates span the symmetric range {[-M..+M]}.
  
  Note that signed {int2_t} variables naturally have an *asymmetric* range
  {[-N..N-1]} where {N = 2^31}. */

typedef struct i2_t { int32_t c[2]; } i2_t;

void i2_zero (i2_t *r);
  /* Sets {r} to the zero vector. */
  
void i2_all (int32_t x, i2_t *r);
  /* Sets all coordinates of {r} to the value {x}. */
  
void i2_axis (uint32_t i, i2_t *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void i2_add (i2_t *a, i2_t *b, i2_t *r);
  /* Sets {r = a + b}. Input MMax: {2^30-1}. Output range: {-2*M..+2*M}. */

void i2_sub (i2_t *a, i2_t *b, i2_t *r);
  /* Sets {r = a - b}. Input MMax: {2^30-1}. Output range: {-2*M..+2*M}. */

void i2_neg (i2_t *a, i2_t *r);
  /* Sets {r} to {-a}. Input MMax: {2^31-1}. Output range: {-M..+M}. */

uint32_t i2_L_inf_norm (i2_t *a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate).
    Input MMax: {2^31-1}. Output range: {0..+M}. */

uint64_t i2_L_inf_dist (i2_t *a, i2_t *b);
  /* Returns the L-infinity distance between {a} and {b} (max absolute diff).
    Input MMax: {2^31-1}. Output range: {0..+2*M}. */

uint64_t i2_norm_sqr (i2_t *a);
  /* Returns the square of the Euclidean norm of {a}. 
    Input MMax: {2^31-1}. Output range: {0..+2*M^2}.
    !!! Perhaps should be {uint64_t}? !!! */

uint64_t i2_dist_sqr (i2_t *a, i2_t *b);
  /* Returns the square of the Euclidean distance between {a} and {b}.
    Input MMax: {2^30-1}. Output range: {0..+8*M^2}.
    !!! Perhaps should be {uint64_t}? !!! */

int64_t i2_dot (i2_t *a, i2_t *b);
  /* Dot product of vectors {a} and {b}. 
    Input MMax: {2^30-1}. Output range: {-2*M^2 .. +2*M^2}. */

void i2_cross (i2_t *a, i2_t *r);
  /* Returns in {r} the vector {a} rotated 90 degrees counterclockwise. 
    Input MMax: {2^31-1}. Output range: {-M .. +M}. */

int64_t i2_det (i2_t *a, i2_t *b);
  /* Returns the determinant of the 2x2 matrix whose rows are {a,b}.  
    Input MMax: {2^30-1}. Output range: {-2*M^2 .. +2*M^2}. */

bool_t i2_eq(i2_t *p, i2_t *q);
  /* True iff points {p} and {q} are identical. */

sign_t i2_cyclic_order(i2_t *a, i2_t *b, i2_t *c);
  /* Returns {+1} if the directions of the vectors {a}, {b}, and {c}, in
    that order, turn counterclockwise around the origin. Returns {-1} if
    they turn clockwise. Returns 0 if any two have the same direction,
    or any one of them is zero. */

void i2_throw_cube (int32_t m, i2_t *r);
  /* Sets {r} to a uniformly random integer point of the 2-cube (square)
    {[-|m| .. +|m|]^2}.  The sign of {m} is ignored. 
    Input MMax: {2^31-1}. Output range: {-|m| .. +|m|}. */

void i2_print (FILE *f, i2_t *a);
  /* Prints {a} on file {f}, with some default format. */

void i2_gen_print (FILE *f, i2_t *a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%d", "(", " ", and ")", respectively. */

/* DERIVED TYPES */

vec_typedef(i2_vec_t,i2_vec,i2_t);
  /* An {i2_vec_t} is a vector of {i2_t}s. */

typedef bool_t i2_pred_t(i2_t *a);
  /* Type of a function that returns a {bool_t} value from an {i2_t} value. */

typedef double i2_double_func_t(i2_t *a);
  /* Type of a function that returns a {double} value from an {i2_t} value. */

typedef i2_t i2_map_t(i2_t *a);
  /* Type of a function that returns a {double} value from an {i2_t} value. */

#endif
