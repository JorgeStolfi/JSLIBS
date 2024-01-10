/* r4.h --- operations on points and vectors of R^4 */

#ifndef R4_H
#define R4_H

#include <stdio.h>

typedef double r4_t[4];

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

#endif
