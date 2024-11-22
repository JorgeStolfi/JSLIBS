/* Last edited on 2024-11-20 15:42:31 by stolfi */

/* Based on VectorN.mg, created  95-02-27 by J. Stolfi. */

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

void in_zero (uint32_t n, int32_t *r)
  { for (int32_t i = 0; i < n; i++) { r[i] = 0; }
  }

void in_all (uint32_t n, int32_t x, int32_t *r)
  { for (int32_t i = 0; i < n; i++) 
      { r[i] = x; }
  }

void in_axis (uint32_t n, uint32_t i, int32_t *r)
  { affirm((i >= 0) && (i < n), "in_axis: bad index");
    for (int32_t j = 0; j < n; j++) { r[j] = 0; }
    r[i] = 1;
  }

void in_copy (uint32_t n, int32_t *a, int32_t *r)
  { for (int32_t i = 0; i < n; i++)
      { r[i] = a[i]; }
  }

void in_add (uint32_t n, int32_t *a, int32_t *b, int32_t *r)
  { for (int32_t i = 0; i < n; i++)
      { r[i] = a[i] + b[i]; }
  }

void in_sub (uint32_t n, int32_t *a, int32_t *b, int32_t *r)
  { for (int32_t i = 0; i < n; i++) 
      { r[i] = a[i] - b[i]; }
  }

void in_neg (uint32_t n, int32_t *a, int32_t *r)
  { for (int32_t i = 0; i < n; i++) 
      { r[i] = - a[i]; }
  }

void in_scale (uint32_t n, int32_t s, int32_t *a, int32_t *r)
  { for (int32_t i = 0; i < n; i++) 
      { r[i] = s * a[i]; }
  }

void in_shift (uint32_t n, int32_t s, int32_t *a, int32_t *r)
  { for (int32_t i = 0; i < n; i++) 
      { r[i] = s + a[i]; }
  }

void in_weigh (uint32_t n, int32_t *a, int32_t *w, int32_t *r)
  { for (int32_t i = 0; i < n; i++) 
      { r[i] = a[i] * w[i]; }
  }

int64_t in_sum (uint32_t n, int32_t *a)
  { int64_t sum = 0;
    for (int32_t i = 0; i < n; i++) { sum += (int64_t)(a[i]); }
    return sum;
  }

uint64_t in_L_inf_dist (uint32_t n, int32_t *a, int32_t *b)
  { uint64_t mag = 0;
    for (int32_t i = 0; i < n; i++) 
      { int64_t ai = (int64_t)a[i];
        int64_t bi = (int64_t)b[i];
        uint64_t mi = (uint64_t)llabs(ai - bi);
        if (mi > mag) { mag = mi; }
      }
    return mag;
  }

int64_t in_dot (uint32_t n, int32_t *a, int32_t *b)
  { int64_t sum = 0.0;
    for (int32_t i = 0; i < n; i++) 
      { int64_t ai = (int64_t)a[i];
        int64_t bi = (int64_t)b[i];
        sum += ai*bi;
      }
    return sum;
  }

void in_throw_cube (uint32_t n, int32_t *r, int32_t a, int32_t b)
  { if (a > b) { int32_t tmp = a; a = b; b = tmp; }
    for (int32_t i = 0; i < n; i++)
      { r[i] = int32_abrandom(a, b); }
  }

void in_print (FILE *f, uint32_t n, int32_t *a)
  { in_gen_print(f, n, a, NULL, NULL, NULL, NULL); }

void in_gen_print 
  ( FILE *f, uint32_t n, int32_t *a, 
    char *fmt, 
    char *lp, char *sep, char *rp
  )
  { if (fmt == NULL) { fmt = "%10d"; }
    if (lp == NULL) { lp = "("; }
    if (sep == NULL) { sep = " "; }
    if (rp == NULL) { rp = ")"; }
    fputs(lp, f);
    for (int32_t i = 0; i < n; i++)
      { if (i > 0) { fputs(sep, f); }
        fprintf(f, fmt, a[i]);
      }
    fputs(rp, f);
  }

int32_t *in_alloc(uint32_t n)
  { int32_t *p = talloc(n, int32_t);
    return p;
  }
