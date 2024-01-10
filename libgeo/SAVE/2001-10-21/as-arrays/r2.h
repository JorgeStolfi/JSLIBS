/* r2.h --- operations on points and vectors of R^2 */
/* Last edited on 2001-10-21 21:17:28 by stolfi */

#ifndef r2_H
#define r2_H

#include <stdio.h>

typedef double r2_t [2];

void r2_zero (r2_t res);
  /* Sets "res" to the zero vector. */
  
void r2_axis (int i, r2_t res);
  /* Sets "res" to the "i"th vector of the canonical basis. */

void r2_add (r2_t a, r2_t b, r2_t res);
  /* Sets "res = a + b". */

void r2_sub (r2_t a, r2_t b, r2_t res);
  /* Sets "res = a - b". */

void r2_scale (double s, r2_t a, r2_t res);
  /* Sets "res := s * a". */

void r2_mix_in (double s, r2_t a, r2_t res);
  /* Sets "res := res + s * a". */

double r2_dist (r2_t a, r2_t b);
  /* Returns Euclidean distance between "a" and "b" */

double r2_orthize (r2_t a, r2_t u, r2_t res);
  /* Sets "res" to the component of "a" that is orthogonal to "u". */
  /* Returns "c" such that "a == res + c * u". */

double r2_normalize_inf (r2_t a);
  /* Normalizes "a" to unit L_infinity norm; returns original length. */

double r2_normalize (r2_t a); 
  /* Normalizes "a" to unit Euclidean length; returns original length. */
  
double r2_dot (r2_t a, r2_t b);
  /* Dot product of vectors "a" and "b" */

double r2_cross (r2_t a, r2_t b);
  /* Returns the cross product of "a" and "b". */

void r2_print (FILE *f, r2_t a);
  /* Prints x on file f. */

void r2_throw_cube (r2_t res);
  /* Sets "res" to a uniformly random point of the [-1 .. +1]^2 cube. */
  
void r2_throw_ball (r2_t res);
  /* Sets "res" to a uniformly random point of the unit ball. */

#endif
