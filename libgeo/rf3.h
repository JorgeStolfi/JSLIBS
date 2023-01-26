/* rf3.h --- operations on points and vectors of R^3 (single-precision version) */
/* Last edited on 2023-01-12 06:50:42 by stolfi */

#ifndef rf3_H
#define rf3_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <interval.h>

typedef struct rf3_t { float c[3]; } rf3_t;

#define INF INFINITY

rf3_t rf3_zero (void);
  /* Returns a zero vector. */
  
rf3_t rf3_all (float x);
  /* Returns a vector with all coordinates set to the value {x}. */
  
rf3_t rf3_axis (int32_t i);
  /* Returns the unit vector of coordinate axis {i}. */

rf3_t rf3_add (rf3_t* const a, rf3_t* const b);
  /* Returns the vector {a + b}. */

rf3_t rf3_sub (rf3_t* const a, rf3_t* const b);
  /* Returns the vector {a - b}. */

rf3_t rf3_neg (rf3_t* const a);
  /* Returns the vector {-a}. */

rf3_t rf3_scale (double s, rf3_t* const a);
  /* Returns the vector {s * a}. */

rf3_t rf3_mix (double s, rf3_t* const a, double t, rf3_t* const b);
  /* Returns the vector {s * a + t * b}. */

void rf3_mix_in (double s, rf3_t* const a, rf3_t* const r);
  /* Sets {r := r + s * a}. */

rf3_t rf3_weigh (rf3_t* const a, rf3_t* const w);
  /* Returns the vector {r} with {r[i] = a[i] * w[i]}. */

rf3_t rf3_unweigh (rf3_t* const a, rf3_t* const w);
  /* Returns the vector {r} with {r[i] = a[i] / w[i]}. */

rf3_t rf3_rot_axis (rf3_t* const a, int32_t i, int32_t j, double ang);
  /* Returns {a} after a rotation that moves axis {i} towards 
    axis {j} by {ang} radians, leaving all other coordinates unchanged. */

rf3_t rf3_rot_gen (rf3_t* const a, rf3_t* const d, double ang);
  /* Returns {a} after a rotation by {ang} radians around the vector {r}
    in the sense of the right-hand rule. */

double rf3_norm (rf3_t* const a);
  /* Returns the Euclidean length of {a}. */

double rf3_norm_sqr (rf3_t* const a);
  /* Returns the square of the Euclidean length of {a}. */

float rf3_L_inf_norm (rf3_t* const a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double rf3_dist (rf3_t* const a, rf3_t* const b);
  /* Returns the Euclidean distance between {a} and {b}. */

double rf3_dist_sqr (rf3_t* const a, rf3_t* const b);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double rf3_L_inf_dist (rf3_t* const a, rf3_t* const b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */

rf3_t rf3_dir (rf3_t* const a, double *normP); 
  /*  Returns {a} scaled to unit Euclidean length. 
    If {normP} is not {NULL}, stores the Euclidean norn of {a} in {*normP}. */
  
rf3_t rf3_L_inf_dir (rf3_t* const a, float *normP); 
  /* Returns {a} scaled to unit L-infinity norm,
    If {normP} is not {NULL}, stores the L-infinity norn of {a} in {*normP}. */

double rf3_dot (rf3_t* const a, rf3_t* const b);
  /* Dot product of vectors {a} and {b}. */

double rf3_cos (rf3_t* const a, rf3_t* const b);
  /* Cosine of angle between vectors {a} and {b}. */

double rf3_sin (rf3_t* const a, rf3_t* const b);
  /* Absolute sine of angle between vectors {a} and {b}. */

double rf3_angle (rf3_t* const a, rf3_t* const b);
  /* Angle between vectors {a} and {b}, in radians. */

rf3_t rf3_cross (rf3_t* const a, rf3_t* const b);
  /* Returns the `cross product' of the vectors {a} and {b}. */
  
double rf3_det (rf3_t* const a, rf3_t* const b, rf3_t* const c);
  /* Returns the determinant of the 3 x 3 matrix whose rows
    are {a}, {b}, and {c}. */

double rf3_decomp (rf3_t* const a, rf3_t* const u, rf3_t* para, rf3_t* perp);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = rf3_dot(a,u)/rf3_dot(u,u)}. Also returns {c}. */

bool_t rf3_is_finite(rf3_t* const p);
  /* True iff all coordinates of {p} are finite (neither {±INF} nor {NAN}). */

bool_t rf3_eq(rf3_t* const p, rf3_t* const q);
  /* True iff points {p} and {q} are identical. */
  
rf3_t rf3_barycenter(int32_t np, rf3_t p[], double w[]);
  /*  Returns the barycenter of all points {p[0..np-1]} with weights {w[0..np-1]}.
    The weights must have positive sum.  Assumes equal weights if {w = NULL}. */

void rf3_bbox(int32_t np, rf3_t p[], interval_t B[], bool_t finite);
  /* Computes the X and Y ranges {B[0..2} of the points 
    {p[0..np-1]}. If {finite} is true, ignores points 
    that have infinite or NAN coordinate(s). */

rf3_t rf3_circumcenter(rf3_t* const a, rf3_t* const b, rf3_t* const c);
  /* The center of the circle passing through the three points {a,b,c}. */

int32_t rf3_orient(rf3_t* const a, rf3_t* const b, rf3_t* const c);
  /* The orientation of the triangle{a,b,c}: {+1} if CCW, {-1} if CW, 0 if flat.
    Note that the result is unreliable if {a,b,c} is nearly flat. */
    
bool_t rf3_incircle (rf3_t* const a, rf3_t* const b, rf3_t* const c, rf3_t* const d);
  /* TRUE iff the point {d} lies inside the circle {a,b,c}. */

rf3_t rf3_throw_cube (void);
  /* Sets {r} to a uniformly random point of the 3-cube (square) {[-1 _ +1]^3}. */
  
rf3_t rf3_throw_dir (void);
  /* Sets {r} to a random direction of {R^3}; that is, a 
    uniformly random point on the unit circle {S^1}. */

rf3_t rf3_throw_ball (void);
  /* Sets {r} to a uniformly random point of the unit N-ball. */

rf3_t rf3_throw_normal (void);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

void rf3_print (FILE *f, rf3_t* const a);
  /* Prints {a} on file {f}, with some default format. */

void rf3_gen_print (FILE *f, rf3_t* const a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%12.6e", "(", " ", and ")", respectively. */

/* DERIVED TYPES */

vec_typedef(rf3_vec_t, rf3_vec, rf3_t);
  /* An {rf3_vec_t} is a vector of {rf3_t}s. */

typedef bool_t rf3_pred_t(rf3_t* const a);
  /* Type of a function that returns a {bool_t} value from an {rf3_t} value. */

typedef double rf3_double_func_t(rf3_t* const a);
  /* Type of a function that returns a {double} value from an {rf3_t} value. */

typedef rf3_t rf3_map_t(rf3_t* const a);
  /* Type of a function that returns an {rf3_t} value from an {rf3_t} value. */
 
#endif
