/* in.h --- operations on points and vectors of Z^n */
/* Last edited on 2024-12-05 10:27:40 by stolfi */
/* Based on VectorN.mg, created  95-02-27 by J. Stolfi. */

#ifndef in_H
#define in_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

void in_zero (uint32_t n, int32_t *r);
  /* Sets {r} to the zero vector. */
  
void in_all (uint32_t n, int32_t x, int32_t *r);
  /* Sets all coordinates of {r} to the value {x}. */
  
void in_axis (uint32_t n, uint32_t i, int32_t *r);
  /* Sets {r} to the {i}th vector of the canonical basis. */

void in_copy (uint32_t n, int32_t *a, int32_t *r);
  /* Copies {a} into {r}; i.e. sets {r[i] = a[i]} for {i} in {0..n-1}. */
  
void in_add (uint32_t n, int32_t *a, int32_t *b, int32_t *r);
  /* Sets {r = a + b}. */

void in_sub (uint32_t n, int32_t *a, int32_t *b, int32_t *r);
  /* Sets {r = a - b}. */

void in_neg (uint32_t n, int32_t *a, int32_t *r);
  /* Sets {r} to {-a}. */

void in_scale (uint32_t n, int32_t s, int32_t *a, int32_t *r);
  /* Sets {r := s * a}. */

void in_shift (uint32_t n, int32_t s, int32_t *a, int32_t *r);
  /* Sets {r[i] := s + a[i]} for {i} in {0..n-1}. */

void in_weigh (uint32_t n, int32_t *a, int32_t *w, int32_t *r);
  /* Sets {r[i] := a[i] * w[i]}. */

int64_t in_sum (uint32_t n, int32_t *a);
  /* Returns the sum of all elements {a[0..n-1]}. */

uint64_t in_L_inf_dist (uint32_t n, int32_t *a, int32_t *b);
  /* Returns the L-infinity distance between {a} and {b} 
    (max absolute diff). */
 
int64_t in_dot (uint32_t n, int32_t *a, int32_t *b);
  /* Dot product of vectors {a} and {b}.  May overflow. */

void in_throw_cube (uint32_t n, int32_t *r, int32_t a, int32_t b);
  /* Sets {r} to a uniformly random point of the {n}-cube {[a .. b]^n}
    (or{[b .. a]^n}, if a > b}). */

void in_print (FILE *f, uint32_t n, int32_t *a);
  /* Prints {a} on file {f}, with some default format. */

void in_gen_print
  ( FILE *f, uint32_t n, int32_t *a, 
    char *fmt, 
    char *lp, char *sep, char *rp
  );
  /* Prints {a} on file {f}, formatting each coordinate with {fmt}. 
    The strings {lp}, {sep}, and {rp} are printed respectively before,
    between, and after all the coordinates of {a}.  When NULL, they default 
    to "%10d", "(", " ", and ")", respectively. */

/* HEAP ALLOCATION */

int32_t *in_alloc(uint32_t n);
  /* Allocates a vector of {n} {int32_t}s on the heap;
    bombs out if no mem. */

/* DERIVED TYPES */

typedef bool_t in_pred_t(uint32_t n, int32_t *a);
  /* Type of a function that returns a {bool_t} value from an {int32_t} vector. */

typedef double in_double_func_t(uint32_t n, int32_t *a);
  /* Type of a function that returns a {double} value from an {int32_t} vector. */

#endif
