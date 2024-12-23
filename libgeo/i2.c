/* See i2.h */
/* Last edited on 2024-11-20 14:10:02 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <i2.h>

#include <jsrandom.h>
#include <affirm.h>
#include <vec.h>

#define N 2

void i2_zero (i2_t *r)
  { r->c[0] = 0;
    r->c[1] = 0;
  }

void i2_all (int32_t x, i2_t *r)
  { r->c[0] = x;
    r->c[1] = x;
  }

void i2_axis (uint32_t i, i2_t *r)
  { demand(i < N, "i2_axis: bad index");
    r->c[0] = 0;
    r->c[1] = 0;
    r->c[i] = 1;
  }

void i2_add (i2_t *a, i2_t *b, i2_t *r)
  { r->c[0] = a->c[0] + b->c[0];
    r->c[1] = a->c[1] + b->c[1];
  }

void i2_sub (i2_t *a, i2_t *b, i2_t *r)
  { r->c[0] = a->c[0] - b->c[0];
    r->c[1] = a->c[1] - b->c[1];
  }

void i2_neg (i2_t *a, i2_t *r)
  { r->c[0] = - a->c[0];
    r->c[1] = - a->c[1];
  }

uint32_t i2_L_inf_norm (i2_t *a)
  { uint32_t d = 0;
    uint32_t a0 = (uint32_t)labs(a->c[0]);
    uint32_t a1 = (uint32_t)labs(a->c[1]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    return d;
  }

uint64_t i2_L_inf_dist (i2_t *a, i2_t *b)
  { uint64_t d = 0;
    uint64_t d0 = (uint64_t)llabs((int64_t)a->c[0] - (int64_t)b->c[0]);
    uint64_t d1 = (uint64_t)llabs((int64_t)a->c[1] - (int64_t)b->c[1]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    return (d);
  }

uint64_t i2_norm_sqr (i2_t *a)
  { int64_t a0 = (int64_t)a->c[0];
    int64_t a1 = (int64_t)a->c[1];
    return (uint64_t)(a0*a0 + a1*a1);
  }

uint64_t i2_dist_sqr (i2_t *a, i2_t *b)
  { int64_t d0 = (int64_t)a->c[0] - (int64_t)b->c[0];
    int64_t d1 = (int64_t)a->c[1] - (int64_t)b->c[1];
    return (uint64_t)(d0*d0 + d1*d1);
  }

int64_t i2_dot (i2_t *a, i2_t *b)
  { int64_t a0 = (int64_t)a->c[0];
    int64_t a1 = (int64_t)a->c[1];

    int64_t b0 = (int64_t)b->c[0];
    int64_t b1 = (int64_t)b->c[1];
    
    return a0*b0 + a1*(int64_t)b1;
  }

void i2_cross (i2_t *a, i2_t *r)
  { int32_t a0 = a->c[0];
    int32_t a1 = a->c[1];
    r->c[0] = - a1;
    r->c[1] =   a0;
  }

int64_t i2_det (i2_t *a, i2_t *b)
  { int64_t a0 = (int64_t)a->c[0];
    int64_t a1 = (int64_t)a->c[1];

    int64_t b0 = (int64_t)b->c[0];
    int64_t b1 = (int64_t)b->c[1];

    return a0*b1 - a1*b0; 
  }

bool_t i2_eq(i2_t *p, i2_t *q)
  { if (p->c[0] != q->c[0]) return FALSE;
    if (p->c[1] != q->c[1]) return FALSE;
    return TRUE;
  }

void i2_throw_cube (int32_t m, i2_t *r)
  { m = abs(m);
    r->c[0] = int32_abrandom(-m, m);
    r->c[1] = int32_abrandom(-m, m);
  }

void i2_print (FILE *f, i2_t *a)
  { i2_gen_print(f, a, NULL, NULL, NULL, NULL); }

void i2_gen_print (FILE *f, i2_t *a, char *fmt, char *lp, char *sep, char *rp)
  { if (fmt == NULL) { fmt = "%d"; }
    if (lp == NULL) { lp = "("; }
    if (sep == NULL) { sep = " "; }
    if (rp == NULL) { rp = ")"; }
    fputs(lp, f);
    for (uint32_t i = 0;  i < N; i++)
      { if (i > 0) { fputs(sep, f); }
        fprintf(f, fmt, a->c[i]);
      }
    fputs(rp, f);
  }

vec_typeimpl(i2_vec_t,i2_vec,i2_t);
