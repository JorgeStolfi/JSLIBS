/* r4.h --- operations on points and vectors of R^4 */
/* Last edited on 2001-09-30 22:12:27 by stolfi */

#ifndef r4_H
#define r4_H

#include <stdio.h>

typedef double r4_t[4];

void r4_zero (r4_t res);
  /* Sets "res" to the zero vector. */
  
void r4_axis (int i, r4_t res);
  /* Sets "res" to the "i"th vector of the canonical basis. */

void r4_add (r4_t a, r4_t b, r4_t res);
  /* Sets "res = a + b". */

void r4_sub (r4_t a, r4_t b, r4_t res);
  /* Sets "res = a - b". */

void r4_scale (double s, r4_t a, r4_t res);
  /* Sets "res := s * a". */

double r4_dist (r4_t a, r4_t b);
  /* Returns Euclidean distance between "a" and "b" */

double r4_orthize (r4_t a, r4_t u, r4_t res);
  /* Sets "res" to the component of "a" that is orthogonal to "u". */
  /* Returns "c" such that "a == res + c * u". */

double r4_normalize_inf (r4_t a);
  /* Normalizes "a" to unit L_infinity norm; returns original length. */

double r4_normalize (r4_t a); 
  /* Normalizes "a" to unit Euclidean length; returns original length. */
  
double r4_dot (r4_t a, r4_t b);
  /* Dot product of vectors "a" and "b" */

void r4_cross (r4_t a, r4_t b, r4_t c, r4_t res);
  /* Sets "res" to the cross product of "a" and "b". */

void r4_print (FILE *f, r4_t a);
  /* Prints x on file f. */

void r4_throw_cube (r4_t res);
  /* Sets "res" to a uniformly random point of the [-1 .. +1]^N cube. */
  
void r4_throw_ball (r4_t res);
  /* Sets "res" to a uniformly random point of the unit ball. */

#endif
