/* See i3.h */
/* Last edited on 2018-03-04 22:54:27 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <vec.h>

#include <i3.h>

#define N 3

void i3_zero (i3_t *r)
  { r->c[0] = 0.0;
    r->c[1] = 0.0;
    r->c[2] = 0.0;
  }

void i3_all (int32_t x, i3_t *r)
  { r->c[0] = x;
    r->c[1] = x;
    r->c[2] = x;
  }

void i3_axis (int i, i3_t *r)
  { affirm((i >= 0) && (i < N), "i3_axis: bad index");
    r->c[0] = 0.0;
    r->c[1] = 0.0;
    r->c[2] = 0.0;

    r->c[i] = 1.0;
  }

void i3_add (i3_t *a, i3_t *b, i3_t *r)
  { r->c[0] = a->c[0] + b->c[0];
    r->c[1] = a->c[1] + b->c[1];
    r->c[2] = a->c[2] + b->c[2];
  }

void i3_sub (i3_t *a, i3_t *b, i3_t *r)
  { r->c[0] = a->c[0] - b->c[0];
    r->c[1] = a->c[1] - b->c[1];
    r->c[2] = a->c[2] - b->c[2];
  }

void i3_neg (i3_t *a, i3_t *r)
  { r->c[0] = - a->c[0];
    r->c[1] = - a->c[1];
    r->c[2] = - a->c[2];
  }

int32_t i3_L_inf_norm (i3_t *a)
  { int32_t d = 0;
    int32_t a0 = abs(a->c[0]);
    int32_t a1 = abs(a->c[1]);
    int32_t a2 = abs(a->c[2]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    if (a2 > d) d = a2;
    return (d);
  }

int32_t i3_L_inf_dist (i3_t *a, i3_t *b)
  { int32_t d = 0;
    int32_t d0 = abs(a->c[0] - b->c[0]);
    int32_t d1 = abs(a->c[1] - b->c[1]);
    int32_t d2 = abs(a->c[2] - b->c[2]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    if (d2 > d) d = d2;
    return (d);
  }

int64_t i3_norm_sqr (i3_t *a)
  { int32_t a0 = a->c[0];
    int32_t a1 = a->c[1];
    int32_t a2 = a->c[2];
    return a0*(int64_t)a0 + a1*(int64_t)a1 + a2*(int64_t)a2;
  }

int64_t i3_dist_sqr (i3_t *a, i3_t *b)
  { int64_t d0 = (a->c[0] - (int64_t)b->c[0]);
    int64_t d1 = (a->c[1] - (int64_t)b->c[1]);
    int64_t d2 = (a->c[2] - (int64_t)b->c[2]);
    return d0*d0 + d1*d1 + d2*d2;
  }

int64_t i3_dot (i3_t *a, i3_t *b)
  { int32_t a0 = a->c[0];
    int32_t a1 = a->c[1];
    int32_t a2 = a->c[2];
    
    int32_t b0 = b->c[0];
    int32_t b1 = b->c[1];
    int32_t b2 = b->c[2];
    
    return a0*(int64_t)b0 + a1*(int64_t)b1 + a2*(int64_t)b2;
  }

void i3_cross (i3_t *a, i3_t *b, i3_t *r)
  { int32_t a0 = a->c[0];
    int32_t a1 = a->c[1];
    int32_t a2 = a->c[2];

    int32_t b0 = b->c[0];
    int32_t b1 = b->c[1];
    int32_t b2 = b->c[2];

    r->c[0] = (int32_t)(a1*(int64_t)b2 - a2*(int64_t)b1);
    r->c[1] = (int32_t)(a2*(int64_t)b0 - a0*(int64_t)b2);
    r->c[2] = (int32_t)(a0*(int64_t)b1 - a1*(int64_t)b0);
  }

int64_t i3_det (i3_t *a, i3_t *b, i3_t *c)
  { int32_t a0 = a->c[0];
    int32_t a1 = a->c[1];
    int32_t a2 = a->c[2];

    int32_t b0 = b->c[0];
    int32_t b1 = b->c[1];
    int32_t b2 = b->c[2];

    int32_t c0 = c->c[0];
    int32_t c1 = c->c[1];
    int32_t c2 = c->c[2];

    int64_t ab0 = a1*(int64_t)b2 - a2*(int64_t)b1;
    int64_t ab1 = a2*(int64_t)b0 - a0*(int64_t)b2;
    int64_t ab2 = a0*(int64_t)b1 - a1*(int64_t)b0;
    
    return ab0*c0 + ab1*c1 + ab2*c2;
  }

int i3_eq(i3_t *p, i3_t *q)
  { if (p->c[0] != q->c[0]) return FALSE;
    if (p->c[1] != q->c[1]) return FALSE;
    if (p->c[2] != q->c[2]) return FALSE;
    return TRUE;
  }

void i3_throw_cube (int32_t m, i3_t *r)
  { r->c[0] = int32_abrandom(-m, m);
    r->c[1] = int32_abrandom(-m, m);
    r->c[2] = int32_abrandom(-m, m);
  }

void i3_print (FILE *f, i3_t *a)
  { i3_gen_print(f, a, NULL, NULL, NULL, NULL); }

void i3_gen_print (FILE *f, i3_t *a, char *fmt, char *lp, char *sep, char *rp)
  { int i;
    if (fmt == NULL) { fmt = "%d"; }
    if (lp == NULL) { lp = "("; }
    if (sep == NULL) { sep = " "; }
    if (rp == NULL) { rp = ")"; }
    fputs(lp, f);
    for (i = 0; i < N; i++)
      { if (i > 0) { fputs(sep, f); }
        fprintf(f, fmt, a->c[i]);
      }
    fputs(rp, f);
  }

vec_typeimpl(i3_vec_t,i3_vec,i3_t);
