/* rf2.h --- operations on points and vectors of R^2 (single-precision version) */
/* Last edited on 2021-08-17 08:51:28 by stolfi */

#ifndef rf2_H
#define rf2_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <interval.h>

typedef struct { float c[2]; } rf2_t;

#define INF INFINITY

rf2_t rf2_zero (void);
  /* Returns a zero vector. */
  
rf2_t rf2_all (float x);
  /* Returns a vector with all coordinates set to the value {x}. */
  
rf2_t rf2_axis (int32_t i);
  /* Returns the unit vector of coordinate axis {i}. */

rf2_t rf2_add (rf2_t *a, rf2_t *b);
  /* Returns the vector {a + b}. */

rf2_t rf2_sub (rf2_t *a, rf2_t *b);
  /* Returns the vector {a - b}. */

rf2_t rf2_neg (rf2_t *a);
  /* Returns the vector {-a}. */

rf2_t rf2_scale (double s, rf2_t *a);
  /* Returns the vector {s * a}. */

rf2_t rf2_mix (double s, rf2_t *a, double t, rf2_t *b);
  /* Returns the vector {s * a + t * b}. */

void rf2_mix_in (double s, rf2_t *a, rf2_t *r);
  /* Sets {r := r + s * a}. */

rf2_t rf2_weigh (rf2_t *a, rf2_t *w);
  /* Returns the vector {r} with {r[i] = a[i] * w[i]}. */

double rf2_norm (rf2_t *a);
  /* Returns the Euclidean length of {a}. */

double rf2_norm_sqr (rf2_t *a);
  /* Returns the square of the Euclidean length of {a}. */

float rf2_L_inf_norm (rf2_t *a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double rf2_dist (rf2_t *a, rf2_t *b);
  /* Returns the Euclidean distance between {a} and {b}. */

double rf2_dist_sqr (rf2_t *a, rf2_t *b);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double rf2_L_inf_dist (rf2_t *a, rf2_t *b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */

rf2_t rf2_dir (rf2_t *a, double *normP); 
  /*  Returns {a} scaled to unit Euclidean length. 
    If {normP} is not {NULL}, stores the Euclidean norn of {a} in {*normP}. */
  
rf2_t rf2_L_inf_dir (rf2_t *a, float *normP); 
  /* Returns {a} scaled to unit L-infinity norm,
    If {normP} is not {NULL}, stores the L-infinity norn of {a} in {*normP}. */

double rf2_dot (rf2_t *a, rf2_t *b);
  /* Dot product of vectors {a} and {b}. */

double rf2_cos (rf2_t *a, rf2_t *b);
  /* Cosine of angle between vectors {a} and {b}. */

double rf2_sin (rf2_t *a, rf2_t *b);
  /* Absolute sine of angle between vectors {a} and {b}. */

double rf2_angle (rf2_t *a, rf2_t *b);
  /* Angle between vectors {a} and {b}, in radians. */

rf2_t rf2_cross (rf2_t *a);
  /* Returns the `cross product' of the vector {a}, namely 
    {a} rotated 90 degrees counterclockwise. */
  
double rf2_det (rf2_t *a, rf2_t *b);
  /* Returns the determinant of the 2 x 2 matrix whose rows
    are {a} and {b}. */

double rf2_decomp (rf2_t *a, rf2_t *u, rf2_t *para, rf2_t *perp);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = rf2_dot(a,u)/rf2_dot(u,u)}. Also returns {c}. */

bool_t rf2_is_finite(rf2_t *p);
  /* True iff all coordinates of {p} are finite (neither {±INF} nor {NAN}). */

bool_t rf2_eq(rf2_t *p, rf2_t *q);
  /* True iff points {p} and {q} are identical. */
  
rf2_t rf2_barycenter(int32_t np, rf2_t p[], double w[]);
  /*  Returns the barycenter of all points {p[0..np-1]} with weights {w[0..np-1]}.
    The weights must have positive sum.  Assumes equal weights if {w = NULL}. */

void rf2_bbox(int32_t np, rf2_t p[], interval_t B[], bool_t finite);
  /* Computes the X and Y ranges {B[0],B[1]} of the points 
    {p.e[0..np-1]}. If {finite} is true, ignores points 
    that have infinite or NAN coordinate(s). */

rf2_t rf2_circumcenter(rf2_t *a, rf2_t *b, rf2_t *c);
  /* The center of the circle passing through the three points {a,b,c}. */

int32_t rf2_orient(rf2_t *a, rf2_t *b, rf2_t *c);
  /* The orientation of the triangle{a,b,c}: {+1} if CCW, {-1} if CW, 0 if flat.
    Note that the result is unreliable if {a,b,c} is nearly flat. */
    
bool_t rf2_incircle (rf2_t *a, rf2_t *b, rf2_t *c, rf2_t *d);
  /* TRUE iff the point {d} lies inside the circle {a,b,c}. */

rf2_t rf2_throw_cube (void);
  /* Sets {r} to a uniformly random point of the 2-cube (square) {[-1 _ +1]^2}. */
  
rf2_t rf2_throw_dir (void);
  /* Sets {r} to a random direction of {R^2}; that is, a 
    uniformly random point on the unit circle {S^1}. */

rf2_t rf2_throw_ball (void);
  /* Sets {r} to a uniformly random point of the unit N-ball. */

rf2_t rf2_throw_normal (void);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

void rf2_print (FILE *f, rf2_t *a);
  /* Prints {a} on file {f}, with some default format. */

void rf2_gen_print (FILE *f, rf2_t *a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%13.6e", "(", " ", and ")", respectively. */

/* DERIVED TYPES */

vec_typedef(rf2_vec_t, rf2_vec, rf2_t);
  /* An {rf2_vec_t} is a vector of {rf2_t}s. */

typedef bool_t rf2_pred_t(rf2_t *a);
  /* Type of a function that returns a {bool_t} value from an {rf2_t} value. */

typedef double rf2_double_func_t(rf2_t *a);
  /* Type of a function that returns a {double} value from an {rf2_t} value. */

typedef rf2_t rf2_map_t(rf2_t *a);
  /* Type of a function that returns an {rf2_t} value from an {rf2_t} value. */
 
#endif
