/* in.h --- operations on points and vectors of Z^n */
/*
  Last edited on 2014-03-24 23:33:38 by stolfilocal
  Based on VectorN.mg, created  95-02-27 by J. Stolfi.
  Last edited by stolfi 
*/

#ifndef in_H
#define in_H

#include <stdio.h>
#include <stdint.h>

void in_zero (int n, int32_t *r);
  /* Sets {r} to the zero vector. */
  
void in_all (int n, int32_t x, int32_t *r);
  /* Sets all coordinates of {r} to the value {x}. */
  
void in_axis (int n, int i, int32_t *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void in_copy (int n, int32_t *a, int32_t *r);
  /* Copies {a} into {r}; i.e. sets {r[i] = a[i]} for {i} in {0..n-1}. */
  
void in_add (int n, int32_t *a, int32_t *b, int32_t *r);
  /* Sets {r = a + b}. */

void in_sub (int n, int32_t *a, int32_t *b, int32_t *r);
  /* Sets {r = a - b}. */

void in_neg (int n, int32_t *a, int32_t *r);
  /* Sets {r} to {-a}. */

void in_scale (int n, int32_t s, int32_t *a, int32_t *r);
  /* Sets {r := s * a}. */

void in_shift (int n, int32_t s, int32_t *a, int32_t *r);
  /* Sets {r[i] := s + a[i]} for {i} in {0..n-1}. */

void in_weigh (int n, int32_t *a, int32_t *w, int32_t *r);
  /* Sets {r[i] := a[i] * w[i]}. */

int64_t in_sum (int n, int32_t *a);
  /* Returns the sum of all elements {a[0..n-1]}. */

int64_t in_L_inf_dist (int n, int32_t *a, int32_t *b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */
 
int64_t in_dot (int n, int32_t *a, int32_t *b);
  /* Dot product of vectors {a} and {b}.  May overflow. */

void in_throw_cube (int n, int32_t *r, int32_t a, int32_t b);
  /* Sets {r} to a uniformly random point of the {n}-cube {[a .. b]^n}. */

void in_print (FILE *f, int n, int32_t *a);
  /* Prints {a} on file {f}, with some default format. */

void in_gen_print
  ( FILE *f, int n, int32_t *a, 
    char *fmt, 
    char *lp, char *sep, char *rp
  );
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%10d", "(", " ", and ")", respectively. */

/* HEAP ALLOCATION */

int32_t *in_alloc(int n);
  /* Allocates {n} {int32_t}s on the heap; bombs out if no mem. */

#endif
