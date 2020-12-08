/* r6.h --- operations on points and vectors of R^6 */
/* Last edited on 2014-03-24 23:33:39 by stolfilocal */

#ifndef r6_H
#define r6_H

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#include <vec.h>

typedef struct { double c[6]; } r6_t;

#define INF INFINITY

void r6_zero (r6_t *r);
  /* Sets {r} to the zero vector. */
  
void r6_all (double x, r6_t *r);
  /* Sets all coordinates of {r} to the value {x}. */
  
void r6_axis (int i, r6_t *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void r6_add (r6_t *a, r6_t *b, r6_t *r);
  /* Sets {r = a + b}. */

void r6_sub (r6_t *a, r6_t *b, r6_t *r);
  /* Sets {r = a - b}. */

void r6_neg (r6_t *a, r6_t *r);
  /* Sets {r} to {-a}. */

void r6_scale (double s, r6_t *a, r6_t *r);
  /* Sets {r := s * a}. */

void r6_mix (double s, r6_t *a, double t, r6_t *b, r6_t *r);
  /* Sets {r := s * a + t * b}. */

void r6_mix_in (double s, r6_t *a, r6_t *r);
  /* Sets {r := r + s * a}. */

void r6_weigh (r6_t *a, r6_t *w, r6_t *r);
  /* Sets {r[i] := a[i] * w[i]}. */

double r6_norm (r6_t *a);
  /* Returns the Euclidean length of {a}. */

double r6_norm_sqr (r6_t *a);
  /* Returns the square of the Euclidean length of {a}. */

double r6_L_inf_norm (r6_t *a);
  /* Returns the L-infinity norm of {a} (max absolute coordinate). */

double r6_dist (r6_t *a, r6_t *b);
  /* Returns the Euclidean distance between {a} and {b}. */

double r6_dist_sqr (r6_t *a, r6_t *b);
  /* Returns the square of the Euclidean distance between {a} and {b}. */

double r6_L_inf_dist (r6_t *a, r6_t *b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */

double r6_dir (r6_t *a, r6_t *r); 
  /* Sets {r} to {a} normalized to unit Euclidean length; 
    returns the original length. */
  
double r6_L_inf_dir (r6_t *a, r6_t *r); 
  /* Sets {r} to the vector {a/r6_L_inf_norm(a)}; 
    returns the original norm. */
  
double r6_dot (r6_t *a, r6_t *b);
  /* Dot product of vectors {a} and {b}. */

double r6_cos (r6_t *a, r6_t *b);
  /* Cosine of angle between vectors {a} and {b}. */

double r6_sin (r6_t *a, r6_t *b);
  /* Absolute sine of angle between vectors {a} and {b}. */

double r6_angle (r6_t *a, r6_t *b);
  /* Angle between vectors {a} and {b}, in radians. */

void r6_cross (r6_t *a, r6_t *b, r6_t *c, r6_t *d, r6_t *e, r6_t *r);
  /* Sets {r} to the cross product of {a,b,c,d,e}. */
  
double r6_det (r6_t *a, r6_t *b, r6_t *c, r6_t *d, r6_t *e, r6_t *f);
  /* Returns the determinant of the 6 x 6 matrix whose rows 
    are {a,b,c,d,e,f}. */

double r6_decomp (r6_t *a, r6_t *u, r6_t *para, r6_t *perp);
  /* Sets {para} and {perp} (when not NULL) to the components of 
    {a} that are parallel and perpendicular to to {u}, respectively.
    Namely, {para = c*u} and {perp = a - c*u}, where 
    {c = r6_dot(a,u)/r6_dot(u,u)}. Also returns {c}. */

int r6_is_finite(r6_t *p);
  /* True iff all coordinates of {p} are finite. */

int r6_eq(r6_t *p, r6_t *q);
  /* True iff points {p} and {q} are identical. */

void r6_throw_cube (r6_t *r);
  /* Sets {r} to a uniformly random point of the 6-cube {[-1 _ +1]^6}. */
  
void r6_throw_dir (r6_t *r);
  /* Sets {r} to a random direction of {R^6}; that is, a 
    uniformly random point on the unit 5-sphere {S^5}. */

void r6_throw_ball (r6_t *r);
  /* Sets {r} to a uniformly random point of the unit N-ball. */

void r6_throw_normal (r6_t *r);
  /* Sets each coordinate {r[i]} to an independent Gaussian random
    number with zero mean and unit standard deviation. */

void r6_print (FILE *f, r6_t *a);
  /* Prints {a} on file {f}, with some default format. */

void r6_gen_print (FILE *f, r6_t *a, char *fmt, char *lp, char *sep, char *rp);
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}.
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%16.8e", "(", " ", and ")", respectively. */

vec_typedef(r6_vec_t,r6_vec,r6_t);
  /* An {r6_vec_t} is a vector of {r6_t}s. */

#endif
