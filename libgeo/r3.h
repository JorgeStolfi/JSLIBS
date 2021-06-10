/* r3.h --- operations on points and vectors of R^3 */
/* Last edited on 2021-06-09 20:43:49 by jstolfi */

#ifndef r3_H
#define r3_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>

typedef struct { double c[3]; } r3_t;

#define INF INFINITY

void r3_zero (r3_t *r);
  /* Sets {r} to the zero vector. */
  
void r3_all (double x, r3_t *r);
  /* Sets all coordinates of {r} to the value {x}. */

void r3_axis (int32_t i, r3_t *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void r3_add (r3_t *a, r3_t *b, r3_t *r);
  /* Sets {r = a + b}. */

void r3_sub (r3_t *a, r3_t *b, r3_t *r);
  /* Sets {r = a - b}. */

void r3_neg (r3_t *a, r3_t *r);
  /* Sets {r} to {-a}. */

void r3_scale (double s, r3_t *a, r3_t *r);
  /* Sets {r := s * a}. */

void r3_mix (double s, r3_t *a, double t, r3_t *b, r3_t *r);
  /* Sets {r := s * a + t * b}. */

void r3_mix_in (double s, r3_t *a, r3_t *r);
  /* Sets {r := r + s * a}. */

void r3_weigh (r3_t *a, r3_t *w, r3_t *r);
  /* Sets {r[i] := a[i] * w[i]}. */

double r3_norm (r3_t *a);
  /* Returns the Euclidean length of {a}. */

double r3_norm_sqr (r3_t *a);
  /* Returns the square of the Euclidean length of {a}. */

double r3_L_inf_norm (r3_t *a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double r3_dist (r3_t *a, r3_t *b);
  /* Returns the Euclidean distance between {a} and {b}. */

double r3_dist_sqr (r3_t *a, r3_t *b);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double r3_L_inf_dist (r3_t *a, r3_t *b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */
  
double r3_dir (r3_t *a, r3_t *r); 
  /* Sets {r} to {a} normalized to unit Euclidean length; 
    returns the original length. */

double r3_L_inf_dir (r3_t *a, r3_t *r); 
  /* Sets {r} to the vector {a/r3_L_inf_norm(a)}; 
    returns the original norm. */

double r3_dot (r3_t *a, r3_t *b);
  /* Dot product of vectors {a} and {b}. */

double r3_cos (r3_t *a, r3_t *b);
  /* Cosine of angle between vectors {a} and {b}. */

double r3_sin (r3_t *a, r3_t *b);
  /* Absolute sine of angle between vectors {a} and {b}. */

double r3_angle (r3_t *a, r3_t *b);
  /* Angle between vectors {a} and {b}, in radians. */

void r3_cross (r3_t *a, r3_t *b, r3_t *r);
  /* Sets {r} to the cross product of {a} and {b}. */

double r3_det (r3_t *a, r3_t *b, r3_t *c);
  /* Returns the determinant of the 3 x 3 matrix whose rows 
    are {a,b,c}. */

double r3_decomp (r3_t *a, r3_t *u, r3_t *para, r3_t *perp);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = r3_dot(a,u)/r3_dot(u,u)}. Also returns {c}. */

int32_t r3_is_finite(r3_t *p);
  /* True iff all coordinates of {p} are finite. */

bool_t r3_eq(r3_t *p, r3_t *q);
  /* True iff points {p} and {q} are identical. */

void r3_throw_cube (r3_t *r);
  /* Sets {r} to a uniformly random point of the 3-cube {[-1 _ +1]^3}. */
  
void r3_throw_dir (r3_t *r);
  /* Sets {r} to a random direction of {R^3}; that is, a 
    uniformly random point on the unit sphere {S^2}. */

void r3_throw_ball (r3_t *r);
  /* Sets {r} to a uniformly random point of the unit 3-ball. */

void r3_throw_normal (r3_t *r);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

void r3_print (FILE *f, r3_t *a);
  /* Prints {a} on file {f}, with some default format. */

void r3_gen_print (FILE *f, r3_t *a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%16.8e", "(", " ", and ")", respectively. */

/* DERIVED TYPES */

vec_typedef(r3_vec_t,r3_vec,r3_t);
  /* An {r3_vec_t} is a vector of {r3_t}s. */

typedef bool_t r3_pred_t(r3_t *a);
  /* Type of a function that returns a {bool_t} value from an {r3_t} value. */

typedef double r3_double_func_t(r3_t *a);
  /* Type of a function that returns a {double} value from an {r3_t} value. */

typedef r3_t r3_map_t(r3_t *a);
  /* Type of a function that returns a {double} value from an {r3_t} value. */

#endif
