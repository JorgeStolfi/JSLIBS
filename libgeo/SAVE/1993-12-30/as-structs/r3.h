/* r3.h --- operations on points and vectors of R^3 */

#ifndef R3_H
#define R3_H

#include <stdio.h>

typedef struct { double c[3]; } r3_t;

void r3_zero (r3_t *res);
  /* Sets "res" to the zero vector. */
  
void r3_axis (int i, r3_t *res);
  /* Sets "res" to the "i"th vector of the canonical basis. */

void r3_add (r3_t *a, r3_t *b, r3_t *res);
  /* Sets "res = a + b". */

void r3_sub (r3_t *a, r3_t *b, r3_t *res);
  /* Sets "res = a - b". */

void r3_scale (double s, r3_t *a, r3_t *res);
  /* Sets "res := s * a". */

void r3_mix_in (double s, r3_t *a, r3_t *res);
  /* Sets "res := res + s * a". */

double r3_dist (r3_t *a, r3_t *b);
  /* Returns Euclidean distance between "a" and "b" */

double r3_orthize (r3_t *a, r3_t *u, r3_t *res);
  /* Sets "res" to the component of "a" that is orthogonal to "u". */
  /* Returns "c" such that "a == res + c * u". */

double r3_normalize_inf (r3_t *a);
  /* Normalizes "a" to unit L_infinity norm; returns original length. */

double r3_normalize (r3_t *a); 
  /* Normalizes "a" to unit Euclidean length; returns original length. */
  
double r3_dot (r3_t *a, r3_t *b);
  /* Dot product of vectors "a" and "b" */

void r3_cross (r3_t *a, r3_t *b, r3_t *res);
  /* Sets "res" to the cross product of "a" and "b". */

void r3_print (FILE *f, r3_t *a);
  /* Prints x on file f. */

void r3_throw_cube (r3_t *res);
  /* Sets "res" to a uniformly random point of the [-1 .. +1]^N cube. */
  
void r3_throw_ball (r3_t *res);
  /* Sets "res" to a uniformly random point of the unit ball. */

#endif
