/* rfn.h --- operations on points and vectors of R^n (single precision version) */
/* Last edited on 2021-08-20 16:09:23 by stolfi */
/* 
  Based on VectorN.mg, created  95-02-27 by J. Stolfi.
*/

#ifndef rfn_H
#define rfn_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

void rfn_zero (int32_t n, float *r);
  /* Sets {r} to the zero vector. */
  
void rfn_all (int32_t n, float x, float *r);
  /* Sets all coordinates of {r} to the value {x}. */
  
void rfn_axis (int32_t n, int32_t i, float *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void rfn_copy (int32_t n, float *a, float *r);
  /* Copies {a} into {r}; i.e. sets {r[i] = a[i]} for {i} in {0..n-1}. */
  
void rfn_add (int32_t n, float *a, float *b, float *r);
  /* Sets {r = a + b}. */

void rfn_sub (int32_t n, float *a, float *b, float *r);
  /* Sets {r = a - b}. */

void rfn_neg (int32_t n, float *a, float *r);
  /* Sets {r} to {-a}. */

void rfn_scale (int32_t n, double s, float *a, float *r);
  /* Sets {r := s * a}. */

void rfn_shift (int32_t n, double s, float *a, float *r);
  /* Sets {r[i] := s + a[i]} for {i} in {0..n-1}. */

void rfn_mix (int32_t n, double s, float *a, double t, float *b, float *r);
  /* Sets {r := s * a + t * b}. */

void rfn_mix_in (int32_t n, double s, float *a, float *r);
  /* Sets {r := r + s * a}. */

void rfn_weigh (int32_t n, float *a, float *w, float *r);
  /* Sets {r[i] := a[i] * w[i]}. */

void rfn_unweigh (int32_t n, float *a, float *w, float *r);
  /* Sets {r[i] := a[i] / w[i]}. */

void rfn_rot_axis (int32_t n, float *a, int32_t i, int32_t j, double ang, float *r);
  /* Sets {r} to {a} after a rotation that moves axis {i} towards 
    axis {j} by {ang} radians, leaving all other coordinates unchanged. */

double rfn_sum (int32_t n, float *a);
  /* Returns the sum of all elements {a[0..n-1]}. */

double rfn_norm (int32_t n, float *a);
  /* Returns the Euclidean length of {a}. */

double rfn_norm_sqr (int32_t n, float *a);
  /* Returns the square of the Euclidean length of {a}. */

float rfn_L_inf_norm (int32_t n, float *a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double rfn_dist (int32_t n, float *a, float *b);
  /* Returns the Euclidean distance between {a} and {b}. */

double rfn_dist_sqr (int32_t n, float *a, float *b);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double rfn_L_inf_dist (int32_t n, float *a, float *b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */
 
double rfn_dir (int32_t n, float *a, float *r);
  /* Sets {r} to {a} normalized to unit Euclidean length; 
    returns the original length. */

float rfn_L_inf_dir (int32_t n, float *a, float *r);
  /* Sets {r} to the vector {a/rfn_L_inf_norm(a)}; 
    returns the original norm. */

double rfn_dot (int32_t n, float *a, float *b);
  /* Dot product of vectors {a} and {b}. */

double rfn_cos (int32_t n, float *a, float *b);
  /* Cosine of angle between vectors {a} and {b}. */

double rfn_sin (int32_t n, float *a, float *b);
  /* Absolute sine of angle between vectors {a} and {b}. */

double rfn_angle (int32_t n, float *a, float *b);
  /* Angle between vectors {a} and {b}, in radians. */

void rfn_cross (int32_t n, float **a, float *r);
  /* Sets {r} to the `cross product' of the {n-1} given {n}-vectors
    {a[0..n-2]}; namely, a vector perpendicular to {a[0..n-2]}, whose
    length is the {n}-dimensional measure of the parallelotope defined
    by those vectors. */
  
double rfn_det (int32_t n, float **a);
  /* Returns the determinant of the {n} x {n} matrix whose rows 
    are tne {n}-vectors {a[0..n-1]}. */

double rfn_decomp (int32_t n, float *a, float *u, float *para, float *perp);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = rfn_dot(n,a,u)/rfn_dot(n,u,u)}. Also returns {c}. */

double rfn_mirror (int32_t n, float *a, float *u, float *r);
  /* Stores into {r} the vector {a} mirrored in the direction
    of the unit vector {u}, namely {a - 2*dot(a,u)*u}.
    Also returns the value of {dot(a,u)}. */

void rfn_throw_cube (int32_t n, float *r);
  /* Sets {r} to a uniformly random point of the {n}-cube {[-1 _ +1]^n}. */
 
void rfn_throw_dir (int32_t n, float *r);
  /* Sets {r} to a random direction of {R^n}; that is, a 
    uniformly random point on {S^{n-1}}, the {(n-1)}-dimensional 
    unit sphere of {R^n}. */

void rfn_throw_ball (int32_t n, float *r);
  /* Sets {r} to a uniformly random point of the unit {n}-ball. */

void rfn_throw_normal (int32_t n, float *r);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

double rfn_abs_rel_diff(int32_t n, float *a, float *b, double abs_tol, double rel_tol);
  /* Computes the maximum difference between each pair {a[i],b[i]},
    divided by {abs_tol} or {rel_tol} times the largest of the two
    elements. See {abs_rel_diff} in {jsmath.h} for details. */

void rfn_print (FILE *f, int32_t n, float *a);
  /* Prints {a} on file {f}, with some default format. */

void rfn_gen_print
  ( FILE *f, int32_t n, float *a, 
    char *fmt, 
    char *lp, char *sep, char *rp
  );
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%16.8e", "(", " ", and ")", respectively. */

/* HEAP ALLOCATION */

float *rfn_alloc(int32_t n);
  /* Allocates {n} {float}s on the heap; bombs out if no mem. */

/* DERIVED TYPES */

typedef bool_t rfn_pred_t(int32_t n, double *a);
  /* Type of a function that returns a {bool_t} value from an {double} array. */

typedef double rfn_double_func_t(int32_t n, double *a);
  /* Type of a function that returns a {double} value from an {double} array. */

#endif
