/* rf4.h --- operations on points and vectors of R^4 (single-precision version) */
/* Last edited on 2023-01-12 06:50:51 by stolfi */

#ifndef rf4_H
#define rf4_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <interval.h>

typedef struct rf4_t { float c[4]; } rf4_t;

#define INF INFINITY

rf4_t rf4_zero (void);
  /* Returns a zero vector. */
  
rf4_t rf4_all (float x);
  /* Returns a vector with all coordinates set to the value {x}. */

rf4_t rf4_axis (int32_t i);
  /* Returns the unit vector of coordinate axis {i}. */

rf4_t rf4_add (rf4_t* const a, rf4_t* const b);
  /* Returns the vector {a + b}. */

rf4_t rf4_sub (rf4_t* const a, rf4_t* const b);
  /* Returns the vector {a - b}. */

rf4_t rf4_neg (rf4_t* const a);
  /* Returns the vector {-a}. */

rf4_t rf4_scale (double s, rf4_t* const a);
  /* Returns the vector {s * a}. */

rf4_t rf4_mix (double s, rf4_t* const a, double t, rf4_t* const b);
  /* Returns the vector {s * a + t * b}. */

void rf4_mix_in (double s, rf4_t* const a, rf4_t* r);
  /* Sets {r := r + s * a}. */

rf4_t rf4_weigh (rf4_t* const a, rf4_t* const w);
  /* Returns the vector {r} with {r[i] = a[i] * w[i]}. */

rf4_t rf4_unweigh (rf4_t* const a, rf4_t* const w);
  /* Returns the vector {r} with {r[i] = a[i] / w[i]}. */

rf4_t rf4_rot_axis (rf4_t* const a, int32_t i, int32_t j, double ang);
  /* Returns {a} after a rotation that moves axis {i} towards 
    axis {j} by {ang} radians, leaving all other coordinates unchanged. */

double rf4_norm (rf4_t* const a);
  /* Returns the Euclidean length of {a}. */

double rf4_norm_sqr (rf4_t* const a);
  /* Returns the square of the Euclidean length of {a}. */

float rf4_L_inf_norm (rf4_t* const a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double rf4_dist (rf4_t* const a, rf4_t* const b);
  /* Returns the Euclidean distance between {a} and {b}. */

double rf4_dist_sqr (rf4_t* const a, rf4_t* const b);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double rf4_L_inf_dist (rf4_t* const a, rf4_t* const b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */

rf4_t rf4_dir (rf4_t* const a, double *normP); 
  /*  Returns {a} scaled to unit Euclidean length. 
    If {normP} is not {NULL}, stores the Euclidean norn of {a} in {*normP}. */
  
rf4_t rf4_L_inf_dir (rf4_t* const a, float *normP); 
  /* Returns {a} scaled to unit L-infinity norm,
    If {normP} is not {NULL}, stores the L-infinity norn of {a} in {*normP}. */

double rf4_dot (rf4_t* const a, rf4_t* const b);
  /* Dot product of vectors {a} and {b}. */

double rf4_cos (rf4_t* const a, rf4_t* const b);
  /* Cosine of angle between vectors {a} and {b}. */

double rf4_sin (rf4_t* const a, rf4_t* const b);
  /* Absolute sine of angle between vectors {a} and {b}. */

double rf4_angle (rf4_t* const a, rf4_t* const b);
  /* Angle between vectors {a} and {b}, in radians. */

rf4_t rf4_cross (rf4_t* const a, rf4_t* const b, rf4_t* const c);
  /* Returns the `cross product' of the vectors {a}, {b}, and {c}. */
  
double rf4_det (rf4_t* const a, rf4_t* const b, rf4_t* const c, rf4_t* const d);
  /* Returns the determinant of the 4 x 4 matrix whose rows
    are {a}, {b}, {c}, and {d}. */

double rf4_decomp (rf4_t* const a, rf4_t* const u, rf4_t* para, rf4_t* perp);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = rf4_dot(a,u)/rf4_dot(u,u)}. Also returns {c}. */

bool_t rf4_is_finite(rf4_t* const p);
  /* True iff all coordinates of {p} are finite (neither {±INF} nor {NAN}). */

bool_t rf4_eq(rf4_t* const p, rf4_t* const q);
  /* True iff points {p} and {q} are identical. */
  
rf4_t rf4_barycenter(int32_t np, rf4_t p[], double w[]);
  /*  Returns the barycenter of all points {p[0..np-1]} with weights {w[0..np-1]}.
    The weights must have positive sum.  Assumes equal weights if {w = NULL}. */

void rf4_bbox(int32_t np, rf4_t p[], interval_t B[], bool_t finite);
  /* Computes the X and Y ranges {B[0..3]} of the points 
    {p[0..np-1]}. If {finite} is true, ignores points 
    that have infinite or NAN coordinate(s). */

rf4_t rf4_throw_cube (void);
  /* Sets {r} to a uniformly random point of the 4-cube (square) {[-1 _ +1]^4}. */
  
rf4_t rf4_throw_dir (void);
  /* Sets {r} to a random direction of {R^4}; that is, a 
    uniformly random point on the unit circle {S^1}. */

rf4_t rf4_throw_ball (void);
  /* Sets {r} to a uniformly random point of the unit N-ball. */

rf4_t rf4_throw_normal (void);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

void rf4_print (FILE *f, rf4_t* const a);
  /* Prints {a} on file {f}, with some default format. */

void rf4_gen_print (FILE *f, rf4_t* const a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%13.6e", "(", " ", and ")", respectively. */

/* DERIVED TYPES */

vec_typedef(rf4_vec_t, rf4_vec, rf4_t);
  /* An {rf4_vec_t} is a vector of {rf4_t}s. */

typedef bool_t rf4_pred_t(rf4_t* const a);
  /* Type of a function that returns a {bool_t} value from an {rf4_t} value. */

typedef double rf4_double_func_t(rf4_t* const a);
  /* Type of a function that returns a {double} value from an {rf4_t} value. */

typedef rf4_t rf4_map_t(rf4_t* const a);
  /* Type of a function that returns an {rf4_t} value from an {rf4_t} value. */
 
#endif
