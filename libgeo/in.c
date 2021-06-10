/*
  Last edited on 2018-03-04 22:54:31 by stolfilocal
  Based on VectorN.mg, created  95-02-27 by J. Stolfi.
  Last edited by stolfi 
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include <in.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_elim.h>

void in_zero (int32_t n, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = 0; }
  }

void in_all (int32_t n, int32_t x, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = x; }
  }

void in_axis (int32_t n, int32_t i, int32_t *r)
  { int32_t j;
    affirm((i >= 0) && (i < n), "in_axis: bad index");
    for (j = 0; j < n; j++) { r[j] = 0; }
    r[i] = 1;
  }

void in_copy (int32_t n, int32_t *a, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = a[i]; }
  }

void in_add (int32_t n, int32_t *a, int32_t *b, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = a[i] + b[i]; }
  }

void in_sub (int32_t n, int32_t *a, int32_t *b, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = a[i] - b[i]; }
  }

void in_neg (int32_t n, int32_t *a, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = - a[i]; }
  }

void in_scale (int32_t n, int32_t s, int32_t *a, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = s * a[i]; }
  }

void in_shift (int32_t n, int32_t s, int32_t *a, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = s + a[i]; }
  }

void in_weigh (int32_t n, int32_t *a, int32_t *w, int32_t *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = a[i] * w[i]; }
  }

int64_t in_sum (int32_t n, int32_t *a)
  { int32_t i;
    int64_t sum = 0;
    for (i = 0; i < n; i++) { sum += (int64_t)(a[i]); }
    return sum;
  }

int64_t in_L_inf_dist (int32_t n, int32_t *a, int32_t *b)
  { int64_t mag = 0;
    int32_t i;
    for (i = 0; i < n; i++) 
      { int64_t ai = a[i];
        int64_t bi = b[i];
        int64_t mi = llabs(ai - bi); if (mi > mag) { mag = mi; } }
    return mag;
  }

int64_t in_dot (int32_t n, int32_t *a, int32_t *b)
  { int64_t sum = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) 
      { int64_t ai = a[i];
        int64_t bi = b[i];
        sum += ai*bi;
      }
    return sum;
  }

void in_throw_cube (int32_t n, int32_t *r, int32_t a, int32_t b)
  { int32_t i;
    for (i = 0; i < n; i++)
      { r[i] = int32_abrandom(a, b); }
  }

void in_print (FILE *f, int32_t n, int32_t *a)
  { in_gen_print(f, n, a, NULL, NULL, NULL, NULL); }

void in_gen_print 
  ( FILE *f, int32_t n, int32_t *a, 
    char *fmt, 
    char *lp, char *sep, char *rp
  )
  { int32_t i;
    if (fmt == NULL) { fmt = "%10d"; }
    if (lp == NULL) { lp = "("; }
    if (sep == NULL) { sep = " "; }
    if (rp == NULL) { rp = ")"; }
    fputs(lp, f);
    for (i = 0; i < n; i++)
      { if (i > 0) { fputs(sep, f); }
        fprintf(f, fmt, a[i]);
      }
    fputs(rp, f);
  }

int32_t *in_alloc(int32_t n)
  { void *p = malloc(n*sizeof(int32_t));
    affirm(p != NULL, "not enough memory");
    return (int32_t *)p;
  }
