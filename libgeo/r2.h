/* r2.h --- operations on points and vectors of R^2 */
/* Last edited on 2021-06-09 20:43:33 by jstolfi */

#ifndef r2_H
#define r2_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <interval.h>

typedef struct { double c[2]; } r2_t;

#define INF INFINITY

void r2_zero (r2_t *r);
  /* Sets {r} to the zero vector. */
  
void r2_all (double x, r2_t *r);
  /* Sets all coordinates of {r} to the value {x}. */
  
void r2_axis (int32_t i, r2_t *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void r2_add (r2_t *a, r2_t *b, r2_t *r);
  /* Sets {r = a + b}. */

void r2_sub (r2_t *a, r2_t *b, r2_t *r);
  /* Sets {r = a - b}. */

void r2_neg (r2_t *a, r2_t *r);
  /* Sets {r} to {-a}. */

void r2_scale (double s, r2_t *a, r2_t *r);
  /* Sets {r := s * a}. */

void r2_mix (double s, r2_t *a, double t, r2_t *b, r2_t *r);
  /* Sets {r := s * a + t * b}. */

void r2_mix_in (double s, r2_t *a, r2_t *r);
  /* Sets {r := r + s * a}. */

void r2_weigh (r2_t *a, r2_t *w, r2_t *r);
  /* Sets {r[i] := a[i] * w[i]}. */

double r2_norm (r2_t *a);
  /* Returns the Euclidean length of {a}. */

double r2_norm_sqr (r2_t *a);
  /* Returns the square of the Euclidean length of {a}. */

double r2_L_inf_norm (r2_t *a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double r2_dist (r2_t *a, r2_t *b);
  /* Returns the Euclidean distance between {a} and {b}. */

double r2_dist_sqr (r2_t *a, r2_t *b);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double r2_L_inf_dist (r2_t *a, r2_t *b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */

double r2_dir (r2_t *a, r2_t *r); 
  /* Sets {r} to {a} normalized to unit Euclidean length; 
    returns the original length. */
  
double r2_L_inf_dir (r2_t *a, r2_t *r); 
  /* Sets {r} to the vector {a/r2_L_inf_norm(a)}; 
    returns the original norm. */

double r2_dot (r2_t *a, r2_t *b);
  /* Dot product of vectors {a} and {b}. */

double r2_cos (r2_t *a, r2_t *b);
  /* Cosine of angle between vectors {a} and {b}. */

double r2_sin (r2_t *a, r2_t *b);
  /* Absolute sine of angle between vectors {a} and {b}. */

double r2_angle (r2_t *a, r2_t *b);
  /* Angle between vectors {a} and {b}, in radians. */

void r2_cross (r2_t *a, r2_t *r);
  /* Sets {r} to the `cross product' of the vector {a}, namely 
    {a} rotated 90 degrees counterclockwise. */
  
double r2_det (r2_t *a, r2_t *b);
  /* Returns the determinant of the 2 x 2 matrix whose rows
    are {a} and {b}. */

double r2_decomp (r2_t *a, r2_t *u, r2_t *para, r2_t *perp);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = r2_dot(a,u)/r2_dot(u,u)}. Also returns {c}. */

int32_t r2_is_finite(r2_t *p);
  /* True iff both coordinates of {p} are finite (neither {±INF} nor {NAN}). */

bool_t r2_eq(r2_t *p, r2_t *q);
  /* True iff points {p} and {q} are identical. */
  
void r2_barycenter(int32_t np, r2_t p[], double w[], r2_t *barP);
  /* The barycenter of all points {p[0..np-1]} with weights {w[0..np-1]}.
    The weights must have positive sum.
    Assumes equal weights if {w = NULL}. */

void r2_bbox(int32_t np, r2_t p[], interval_t B[], bool_t finite);
  /* Computes the X and Y ranges {B[0],B[1]} of the points 
    {p.e[0..np-1]}. If {finite} is true, ignores points 
    that have infinite or NAN coordinate(s). */

r2_t r2_circumcenter(r2_t *a, r2_t *b, r2_t *c);
  /* The center of the circle passing through the three points {a,b,c}. */

int32_t r2_orient(r2_t *a, r2_t *b, r2_t *c);
  /* The orientation of the triangle{a,b,c}: {+1} if CCW, {-1} if CW, 0 if flat.
    Note that the result is unreliable if {a,b,c} is nearly flat. */
    
bool_t r2_incircle (r2_t *a, r2_t *b, r2_t *c, r2_t *d);
  /* TRUE iff the point {d} lies inside the circle {a,b,c}. */

void r2_throw_cube (r2_t *r);
  /* Sets {r} to a uniformly random point of the 2-cube (square) {[-1 _ +1]^2}. */
  
void r2_throw_dir (r2_t *r);
  /* Sets {r} to a random direction of {R^2}; that is, a 
    uniformly random point on the unit circle {S^1}. */

void r2_throw_ball (r2_t *r);
  /* Sets {r} to a uniformly random point of the unit N-ball. */

void r2_throw_normal (r2_t *r);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

void r2_print (FILE *f, r2_t *a);
  /* Prints {a} on file {f}, with some default format. */

void r2_gen_print (FILE *f, r2_t *a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%16.8e", "(", " ", and ")", respectively. */

/* DERIVED TYPES */

vec_typedef(r2_vec_t, r2_vec, r2_t);
  /* An {r2_vec_t} is a vector of {r2_t}s. */

typedef bool_t r2_pred_t(r2_t *a);
  /* Type of a function that returns a {bool_t} value from an {r2_t} value. */

typedef double r2_double_func_t(r2_t *a);
  /* Type of a function that returns a {double} value from an {r2_t} value. */

typedef r2_t r2_map_t(r2_t *a);
  /* Type of a function that returns a {double} value from an {r2_t} value. */
 
#endif
