/* rn.h --- operations on points and vectors of R^n */
/* Last edited on 2025-02-05 15:45:12 by stolfi */
/* 
  Based on VectorN.mg, created  95-02-27 by J. Stolfi.
*/

#ifndef rn_H
#define rn_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

void rn_zero (uint32_t n, double r[]);
  /* Sets {r} to the zero vector. */
  
void rn_all (uint32_t n, double x, double r[]);
  /* Sets all coordinates of {r} to the value {x}. */
  
void rn_axis (uint32_t n, uint32_t i, double r[]);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void rn_copy (uint32_t n, double a[], double r[]);
  /* Copies {a} into {r}; i.e. sets {r[i] = a[i]} for {i} in {0..n-1}. */
  
void rn_add (uint32_t n, double a[], double b[], double r[]);
  /* Sets {r = a + b}. */

void rn_sub (uint32_t n, double a[], double b[], double r[]);
  /* Sets {r = a - b}. */

void rn_neg (uint32_t n, double a[], double r[]);
  /* Sets {r} to {-a}. */

void rn_scale (uint32_t n, double s, double a[], double r[]);
  /* Sets {r := s * a}. */

void rn_shift (uint32_t n, double s, double a[], double r[]);
  /* Sets {r[i] := s + a[i]} for {i} in {0..n-1}. */

void rn_mix (uint32_t n, double s, double a[], double t, double b[], double r[]);
  /* Sets {r := s * a + t * b}. */

void rn_mix_in (uint32_t n, double s, double a[], double r[]);
  /* Sets {r := r + s * a}. */

void rn_weigh (uint32_t n, double a[], double w[], double r[]);
  /* Sets {r[i] := a[i] * w[i]}. */

void rn_unweigh (uint32_t n, double a[], double w[], double r[]);
  /* Sets {r[i] := a[i] / w[i]}.  */
  
void rn_rot_axis (uint32_t n, double a[], uint32_t i, uint32_t j, double ang, double r[]);
  /* Sets {r} to {a} after a rotation that moves axis {i} towards 
    axis {j} by {ang} radians, leaving all other coordinates unchanged. */

double rn_sum (uint32_t n, double a[]);
  /* Returns the sum of all elements {a[0..n-1]}. */

double rn_norm (uint32_t n, double a[]);
  /* Returns the Euclidean length of {a}. */

double rn_norm_sqr (uint32_t n, double a[]);
  /* Returns the square of the Euclidean length of {a}. */

double rn_L_inf_norm (uint32_t n, double a[]);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double rn_dist (uint32_t n, double a[], double b[]);
  /* Returns the Euclidean distance between {a} and {b}. */

double rn_dist_sqr (uint32_t n, double a[], double b[]);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double rn_L_inf_dist (uint32_t n, double a[], double b[]);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */

double rn_abs_rel_diff (uint32_t n, double a[], double b[], double abs_tol, double rel_tol);
  /* Computes the maximum difference between each pair {a[i],b[i]},
    divided by {abs_tol} or by {rel_tol} times the largest of the two
    elements. See {abs_rel_diff} in {jsmath.h} for details. */
 
double rn_dir (uint32_t n, double a[], double r[]);  
  /* Sets {r} to {a} normalized to unit Euclidean length.
    If the Euclidean length of {a} is zero, sets {r} to all {NAN}s.
    Returns the original Euclidean length of {a}. */

double rn_L_inf_dir (uint32_t n, double a[], double r[]);
  /* Sets {r} to the vector {a} divided by {rn_L_inf_norm(n,a)},
    the max absolute value of any coordinate.
    If that denominator is zero, sets {r} to all {NAN}s.
    Returns the original value of {rn_L_inf_norm(n,a)}.  */

double rn_dot (uint32_t n, double a[], double b[]);
  /* Dot product of vectors {a} and {b}. */

double rn_cos (uint32_t n, double a[], double b[]);
  /* Cosine of angle between vectors {a} and {b}. */

double rn_sin (uint32_t n, double a[], double b[]);
  /* Absolute sine of angle between vectors {a} and {b}. */

double rn_angle (uint32_t n, double a[], double b[]);
  /* Angle between vectors {a} and {b}, in radians. */

void rn_cross (uint32_t n, double *a[], double r[]);
  /* Sets {r} to the `cross product' of the {n-1} given {n}-vectors
    {a[0..n-2]}; namely, a vector perpendicular to {a[0..n-2]}, whose
    length is the {n}-dimensional measure of the parallelotope defined
    by those vectors. */
  
double rn_det (uint32_t n, double *a[]);
  /* Returns the determinant of the {n} x {n} matrix whose rows 
    are tne {n}-vectors {a[0..n-1]}. */

double rn_decomp (uint32_t n, double a[], double u[], double para[], double perp[]);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = rn_dot(n,a,u)/rn_dot(n,u,u)}. Also returns {c}. */

double rn_mirror (uint32_t n, double a[], double u[], double r[]);
  /* Stores into {r} the vector {a} mirrored in the direction
    of the unit vector {u}, namely {a - 2*dot(a,u)*u}.
    Also returns the value of {dot(a,u)}. */
    
uint32_t rn_remove_zeros(uint32_t n, double a[], double r[]);
  /* Finds the number {m} of non-zero elements of {a[0..n-1]}, copies
    them into {r[0..m-1]}, and returns {m}. Caller must make sure
    that the array {r} has at least {m} elements. */
    
uint32_t rn_insert_zeros(uint32_t n, double a[], double b[], double r[]);
  /* Copies each element of {a[0..n-1]} that is zero into the corrsponding
    element of {r[0..n-1]}, filling the other {m} elements of {r} with {b[0..m-1]}.
    in order. Returns the number {m}. Caller must make sure
    that the array {b} has at least {m} elements.  */ 

/* NON-UNIFORM DISTANCE FUNCTIONS */

double rn_rad_rel_max_diff (uint32_t n, double a[], double b[], double rad[]);
  /* Given two {n}-vectors {a,b}, returns the maximum
    differences between the coordinates, relative to the radius vector {rad[0..n-1]}.
    That is, returns the maximum of {fabs(a[i] - b[i])/rad[i]} for all {i} in
    {0..n-1}. 
    
    The radii {rad[0..n-1]} must be non-negative.  If any {rad[i]} is
    zero, the elements {a[i]} and {b[i]} must be equal, and are ignored.
    If all {rad}s are zero, or {n} is zero, the result is zero.
    
    If {a} or {b} is {NULL}, assumes a vector of {n} zeros.  If {rad}
    is {NULL}, assumes a vector of ones -- so that the result is just
    {rn_L_inf_dist(n,a,b)}. 
    
    Thus the axis-aligned, origin-centered box with with radius {rad[i]}
    along axis {i} consists of all vectors {v} of {\RR^n} such that
    {rn_rad_rel_max_diff(n,v,NULL,rad) <= 1}. */

double rn_rad_rel_dist_sqr (uint32_t n, double a[], double b[], double rad[]);
  /* Given two {n}-vectors {a,b}, returns the total squared coordinate
    differences between them, relative to the radius vector {rad[0..n-1]}. That
    is, returns the sum of {((a[i] - b[i])/rad[i])^2} for all {i} in
    {0..n-1}. 
    
    If {a} or {b} is {NULL}, assumes a vector of {n} zeros.  If {rad}
    is {NULL}, assumes a vector of ones -- so that the result is just
    {rn_dist_sqr(n,a,b)}. 
    
    Thus the axis-aligned, origin-centered ellipsoid with radius {rad[i]}
    along axis {i} consists of all vectors {a} such 
    that {rn_rad_rel_dist_sqr(n,a,NULL,rad) <= 1}. */

/* RANDOM VECTOR GENERATION */

void rn_throw_cube (uint32_t n, double r[]);
  /* Sets {r} to a uniformly random point of the {n}-cube {[-1 _ +1]^n}. */
 
void rn_throw_dir (uint32_t n, double r[]);
  /* Sets {r} to a random direction of {R^n}; that is, a 
    uniformly random point on {S^{n-1}}, the {(n-1)}-dimensional 
    unit sphere of {R^n}. */

void rn_throw_ball (uint32_t n, double r[]);
  /* Sets {r} to a uniformly random point of the unit {n}-ball. */

void rn_throw_normal (uint32_t n, double r[]);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

void rn_print (FILE *wr, uint32_t n, double a[]);
  /* Prints {a} on file {wr}, with some default format. The printout does
    NOT end with newline. */

void rn_gen_print
  ( FILE *wr, uint32_t n, double a[], 
    char *fmt, 
    char *lp, char *sep, char *rp
  );
  /* Prints {a} on file {wr}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%16.8e", "(", " ", and ")", respectively. */

void rn_rad_rel_print(FILE *wr, uint32_t n, double a[], double rad[]);
  /* Prints to {wr} three lines with the vector {a[0..n-1]}, the radii {rad[0..n-1]},
    and each {a[i]} divided by the corresponding {rad[i]}. */

/* HEAP ALLOCATION */

double *rn_alloc(uint32_t n);
  /* Allocates {n} {double}s on the heap; bombs out if no mem. */

/* DERIVED TYPES */

typedef bool_t rn_pred_t(uint32_t n, double a[]);
  /* Type of a function that returns a {bool_t} value from an {double} array. */

typedef double rn_double_func_t(uint32_t n, double a[]);
  /* Type of a function that returns a {double} value from an {double} array. */

#endif
