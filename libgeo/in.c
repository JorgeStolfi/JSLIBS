/*
  Last edited on 2018-03-04 22:54:31 by stolfilocal
  Based on VectorN.mg, created  95-02-27 by J. Stolfi.
  Last edited by stolfi 
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <in.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_elim.h>

void in_zero (int n, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = 0; }
  }

void in_all (int n, int32_t x, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = x; }
  }

void in_axis (int n, int i, int32_t *r)
  { int j;
    affirm((i >= 0) && (i < n), "in_axis: bad index");
    for (j = 0; j < n; j++) { r[j] = 0; }
    r[i] = 1;
  }

void in_copy (int n, int32_t *a, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = a[i]; }
  }

void in_add (int n, int32_t *a, int32_t *b, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = a[i] + b[i]; }
  }

void in_sub (int n, int32_t *a, int32_t *b, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = a[i] - b[i]; }
  }

void in_neg (int n, int32_t *a, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = - a[i]; }
  }

void in_scale (int n, int32_t s, int32_t *a, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = s * a[i]; }
  }

void in_shift (int n, int32_t s, int32_t *a, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = s + a[i]; }
  }

void in_weigh (int n, int32_t *a, int32_t *w, int32_t *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = a[i] * w[i]; }
  }

int64_t in_sum (int n, int32_t *a)
  { int i;
    int64_t sum = 0;
    for (i = 0; i < n; i++) { sum += (int64_t)(a[i]); }
    return sum;
  }

int64_t in_L_inf_dist (int n, int32_t *a, int32_t *b)
  { int64_t mag = 0;
    int i;
    for (i = 0; i < n; i++) 
      { int64_t ai = a[i];
        int64_t bi = b[i];
        int64_t mi = llabs(ai - bi); if (mi > mag) { mag = mi; } }
    return mag;
  }

int64_t in_dot (int n, int32_t *a, int32_t *b)
  { int64_t sum = 0.0;
    int i;
    for (i = 0; i < n; i++) 
      { int64_t ai = a[i];
        int64_t bi = b[i];
        sum += ai*bi;
      }
    return sum;
  }

void in_throw_cube (int n, int32_t *r, int32_t a, int32_t b)
  { int i;
    for (i = 0; i < n; i++)
      { r[i] = int32_abrandom(a, b); }
  }

void in_print (FILE *f, int n, int32_t *a)
  { in_gen_print(f, n, a, NULL, NULL, NULL, NULL); }

void in_gen_print 
  ( FILE *f, int n, int32_t *a, 
    char *fmt, 
    char *lp, char *sep, char *rp
  )
  { int i;
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

int32_t *in_alloc(int n)
  { void *p = malloc(n*sizeof(int32_t));
    affirm(p != NULL, "not enough memory");
    return (int32_t *)p;
  }
