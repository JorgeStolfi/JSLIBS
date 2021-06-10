/* r4.h --- operations on points and vectors of R^4 */
/* Last edited on 2021-06-09 20:44:03 by jstolfi */

#ifndef r4_H
#define r4_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>

typedef struct { double c[4]; } r4_t;

#define INF INFINITY

void r4_zero (r4_t *r);
  /* Sets {r} to the zero vector. */
  
void r4_all (double x, r4_t *r);
  /* Sets all coordinates of {r} to the value {x}. */
  
void r4_axis (int32_t i, r4_t *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void r4_add (r4_t *a, r4_t *b, r4_t *r);
  /* Sets {r = a + b}. */

void r4_sub (r4_t *a, r4_t *b, r4_t *r);
  /* Sets {r = a - b}. */

void r4_neg (r4_t *a, r4_t *r);
  /* Sets {r} to {-a}. */

void r4_scale (double s, r4_t *a, r4_t *r);
  /* Sets {r := s * a}. */

void r4_mix (double s, r4_t *a, double t, r4_t *b, r4_t *r);
  /* Sets {r := s * a + t * b}. */

void r4_mix_in (double s, r4_t *a, r4_t *r);
  /* Sets {r := r + s * a}. */

void r4_weigh (r4_t *a, r4_t *w, r4_t *r);
  /* Sets {r[i] := a[i] * w[i]}. */

double r4_norm (r4_t *a);
  /* Returns the Euclidean length of {a}. */

double r4_norm_sqr (r4_t *a);
  /* Returns the square of the Euclidean length of {a}. */

double r4_L_inf_norm (r4_t *a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double r4_dist (r4_t *a, r4_t *b);
  /* Returns the Euclidean distance between {a} and {b}. */

double r4_dist_sqr (r4_t *a, r4_t *b);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double r4_L_inf_dist (r4_t *a, r4_t *b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */

double r4_dir (r4_t *a, r4_t *r); 
  /* Sets {r} to {a} normalized to unit Euclidean length; 
    returns the original length. */
  
double r4_L_inf_dir (r4_t *a, r4_t *r); 
  /* Sets {r} to the vector {a/r4_L_inf_norm(a)}; 
    returns the original norm. */
  
double r4_dot (r4_t *a, r4_t *b);
  /* Dot product of vectors {a} and {b}. */

double r4_cos (r4_t *a, r4_t *b);
  /* Cosine of angle between vectors {a} and {b}. */

double r4_sin (r4_t *a, r4_t *b);
  /* Absolute sine of angle between vectors {a} and {b}. */

double r4_angle (r4_t *a, r4_t *b);
  /* Angle between vectors {a} and {b}, in radians. */

void r4_cross (r4_t *a, r4_t *b, r4_t *c, r4_t *r);
  /* Sets {r} to the cross product of {a}, {b}, and {c}. */
  
double r4_det (r4_t *a, r4_t *b, r4_t *c, r4_t *d);
  /* Returns the determinant of the 4 x 4 matrix whose rows 
    are {a,b,c,d}. */

double r4_decomp (r4_t *a, r4_t *u, r4_t *para, r4_t *perp);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = r4_dot(a,u)/r4_dot(u,u)}. Also returns {c}. */

int32_t r4_is_finite(r4_t *p);
  /* True iff all coordinates of {p} are finite. */

bool_t r4_eq(r4_t *p, r4_t *q);
  /* True iff points {p} and {q} are identical. */

void r4_throw_cube (r4_t *r);
  /* Sets {r} to a uniformly random point of the 4-cube {[-1 _ +1]^4}. */
  
void r4_throw_dir (r4_t *r);
  /* Sets {r} to a random direction of {R^4}; that is, a 
    uniformly random point on the unit 3-sphere {S^3}. */

void r4_throw_ball (r4_t *r);
  /* Sets {r} to a uniformly random point of the unit N-ball. */

void r4_throw_normal (r4_t *r);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

void r4_print (FILE *f, r4_t *a);
  /* Prints {a} on file {f}, with some default format. */

void r4_gen_print (FILE *f, r4_t *a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}.
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%16.8e", "(", " ", and ")", respectively. */

/* DERIVED TYPES */

vec_typedef(r4_vec_t,r4_vec,r4_t);
  /* An {r4_vec_t} is a vector of {r4_t}s. */

typedef bool_t r4_pred_t(r4_t *a);
  /* Type of a function that returns a {bool_t} value from an {r4_t} value. */

typedef double r4_double_func_t(r4_t *a);
  /* Type of a function that returns a {double} value from an {r4_t} value. */

typedef r4_t r4_map_t(r4_t *a);
  /* Type of a function that returns a {double} value from an {r4_t} value. */

#endif
