/* See r6.h. */
/* Last edited on 2021-08-20 16:10:13 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <r6.h>

#include <jsrandom.h>
#include <affirm.h>
#include <rn.h>
#include <vec.h>

#define N 6

void r6_zero (r6_t *r)
  { r->c[0] = 0.0;
    r->c[1] = 0.0;
    r->c[2] = 0.0;
    r->c[3] = 0.0;
    r->c[4] = 0.0;
    r->c[5] = 0.0;
  }

void r6_all (double x, r6_t *r)
  { r->c[0] = x;
    r->c[1] = x;
    r->c[2] = x;
    r->c[3] = x;
    r->c[4] = x;
    r->c[5] = x;
  }

void r6_axis (int32_t i, r6_t *r)
  { affirm((i >= 0) && (i < N), "r6_axis: bad index");
    r->c[0] = 0.0;
    r->c[1] = 0.0;
    r->c[2] = 0.0;
    r->c[3] = 0.0;
    r->c[4] = 0.0;
    r->c[5] = 0.0;

    r->c[i] = 1.0;
  }

void r6_add (r6_t *a, r6_t *b, r6_t *r)
  { r->c[0] = a->c[0] + b->c[0];
    r->c[1] = a->c[1] + b->c[1];
    r->c[2] = a->c[2] + b->c[2];
    r->c[3] = a->c[3] + b->c[3];
    r->c[4] = a->c[4] + b->c[4];
    r->c[5] = a->c[5] + b->c[5];
  }

void r6_sub (r6_t *a, r6_t *b, r6_t *r)
  { r->c[0] = a->c[0] - b->c[0];
    r->c[1] = a->c[1] - b->c[1];
    r->c[2] = a->c[2] - b->c[2];
    r->c[3] = a->c[3] - b->c[3];
    r->c[4] = a->c[4] - b->c[4];
    r->c[5] = a->c[5] - b->c[5];
  }

void r6_neg (r6_t *a, r6_t *r)
  { r->c[0] = - a->c[0];
    r->c[1] = - a->c[1];
    r->c[2] = - a->c[2];
    r->c[3] = - a->c[3];
    r->c[4] = - a->c[4];
    r->c[5] = - a->c[5];
  }

void r6_scale (double s, r6_t *a, r6_t *r)
  { r->c[0] = s * a->c[0];
    r->c[1] = s * a->c[1];
    r->c[2] = s * a->c[2];
    r->c[3] = s * a->c[3];
    r->c[4] = s * a->c[4];
    r->c[5] = s * a->c[5];
  }

void r6_mix (double s, r6_t *a, double t, r6_t *b, r6_t *r)
  { r->c[0] = s * a->c[0] + t * b->c[0];
    r->c[1] = s * a->c[1] + t * b->c[1];
    r->c[2] = s * a->c[2] + t * b->c[2];
    r->c[3] = s * a->c[3] + t * b->c[3];
    r->c[4] = s * a->c[4] + t * b->c[4];
    r->c[5] = s * a->c[5] + t * b->c[5];
  }

void r6_mix_in (double s, r6_t *a, r6_t *r)
  { r->c[0] += s * a->c[0];
    r->c[1] += s * a->c[1];
    r->c[2] += s * a->c[2];
    r->c[3] += s * a->c[3];
    r->c[4] += s * a->c[4];
    r->c[5] += s * a->c[5];
  }

void r6_weigh (r6_t *a, r6_t *w, r6_t *r)
  { r->c[0] = a->c[0] * w->c[0];
    r->c[1] = a->c[1] * w->c[1];
    r->c[2] = a->c[2] * w->c[2];
    r->c[3] = a->c[3] * w->c[3];
    r->c[4] = a->c[4] * w->c[4];
    r->c[5] = a->c[5] * w->c[5];
  }

void r6_unweigh (r6_t *a, r6_t *w, r6_t *r)
  { r->c[0] = a->c[0] / w->c[0];
    r->c[1] = a->c[1] / w->c[1];
    r->c[2] = a->c[2] / w->c[2];
    r->c[3] = a->c[3] / w->c[3];
    r->c[4] = a->c[4] / w->c[4];
    r->c[5] = a->c[5] / w->c[5];
  }

void r6_rot_axis (r6_t *a, int32_t i, int32_t j, double ang, r6_t *r)
  {
    affirm((i >= 0) && (i < N), "r6_rot_axis: bad index {i}");
    affirm((j >= 0) && (j < N), "r6_rot_axis: bad index {j}");
    affirm(i != j, "r6_rot_axis: axes not distinct");
    (*r) = (*a);
    double c = cos(ang);
    double s = sin(ang);
    double x = + c*a->c[i] - s*a->c[j];
    double y = + s*a->c[i] + c*a->c[j];
    r->c[i] = x;
    r->c[j] = y;
  }

double r6_norm (r6_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    double a3 = a->c[3];
    double a4 = a->c[4];
    double a5 = a->c[5];
    return sqrt(a0*a0 + a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5);
  }

double r6_norm_sqr (r6_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    double a3 = a->c[3];
    double a4 = a->c[4];
    double a5 = a->c[5];
    return a0*a0 + a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5;
  }

double r6_L_inf_norm (r6_t *a)
  { double d = 0.0;
    double a0 = fabs(a->c[0]);
    double a1 = fabs(a->c[1]);
    double a2 = fabs(a->c[2]);
    double a3 = fabs(a->c[3]);
    double a4 = fabs(a->c[4]);
    double a5 = fabs(a->c[5]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    if (a2 > d) d = a2;
    if (a3 > d) d = a3;
    if (a4 > d) d = a4;
    if (a5 > d) d = a5;
    return (d);
  }

double r6_dist (r6_t *a, r6_t *b)
  { double d0 = (a->c[0] - b->c[0]);
    double d1 = (a->c[1] - b->c[1]);
    double d2 = (a->c[2] - b->c[2]);
    double d3 = (a->c[3] - b->c[3]);
    double d4 = (a->c[4] - b->c[4]);
    double d5 = (a->c[5] - b->c[5]);
    double d = sqrt(d0*d0 + d1*d1 + d2*d2 + d3*d3 + d4*d4 + d5*d5);
    return (d);
  }

double r6_dist_sqr (r6_t *a, r6_t *b)
  { double d0 = (a->c[0] - b->c[0]);
    double d1 = (a->c[1] - b->c[1]);
    double d2 = (a->c[2] - b->c[2]);
    double d3 = (a->c[3] - b->c[3]);
    double d4 = (a->c[4] - b->c[4]);
    double d5 = (a->c[5] - b->c[5]);
    return d0*d0 + d1*d1 + d2*d2 + d3*d3 + d4*d4 + d5*d5;
  }

double r6_L_inf_dist (r6_t *a, r6_t *b)
  { double d = 0.0;
    double d0 = fabs(a->c[0] - b->c[0]);
    double d1 = fabs(a->c[1] - b->c[1]);
    double d2 = fabs(a->c[2] - b->c[2]);
    double d3 = fabs(a->c[3] - b->c[3]);
    double d4 = fabs(a->c[4] - b->c[4]);
    double d5 = fabs(a->c[5] - b->c[5]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    if (d2 > d) d = d2;
    if (d3 > d) d = d3;
    if (d4 > d) d = d4;
    if (d5 > d) d = d5;
    return (d);
  }

double r6_dir (r6_t *a, r6_t *r)
  { double d = r6_norm(a);
    r->c[0] = a->c[0]/d;
    r->c[1] = a->c[1]/d;
    r->c[2] = a->c[2]/d;
    r->c[3] = a->c[3]/d;
    r->c[4] = a->c[4]/d;
    r->c[5] = a->c[5]/d;
    return (d);
  }

double r6_L_inf_dir (r6_t *a, r6_t *r)
  { double d = r6_L_inf_norm(a);
    r->c[0] = a->c[0]/d;
    r->c[1] = a->c[1]/d;
    r->c[2] = a->c[2]/d;
    r->c[3] = a->c[3]/d;
    r->c[4] = a->c[4]/d;
    r->c[5] = a->c[5]/d;
    return (d);
  }

double r6_dot (r6_t *a, r6_t *b)
  { return 
      a->c[0]*b->c[0] + 
      a->c[1]*b->c[1] + 
      a->c[2]*b->c[2] + 
      a->c[3]*b->c[3] + 
      a->c[4]*b->c[4] + 
      a->c[5]*b->c[5];
  }

double r6_cos (r6_t *a, r6_t *b)
  { return rn_cos(N, &(a->c[0]), &(b->c[0])); }

double r6_sin (r6_t *a, r6_t *b)
  { return rn_sin(N, &(a->c[0]), &(b->c[0])); }

double r6_angle (r6_t *a, r6_t *b)
  { return rn_angle(N, &(a->c[0]), &(b->c[0])); }

void r6_cross (r6_t *a, r6_t *b, r6_t *c, r6_t *d, r6_t *e, r6_t *r)
  { double *v[5] = {a->c,b->c,c->c,d->c,e->c};
    rn_cross(6, v, r->c);
  }

double r6_det (r6_t *a, r6_t *b, r6_t *c, r6_t *d, r6_t *e, r6_t *f)
  { double *v[6] = {a->c,b->c,c->c,d->c,e->c,f->c};
    return rn_det(6, v);
  }

double r6_decomp (r6_t *a, r6_t *u, r6_t *para, r6_t *perp)
  { double *paran = (para == NULL ? NULL : &(para->c[0]));
    double *perpn = (perp == NULL ? NULL : &(perp->c[0]));
    return rn_decomp(N, &(a->c[0]), &(u->c[0]), paran, perpn);
  }

bool_t r6_is_finite(r6_t *p)
  { if (fabs(p->c[0]) == INF) return FALSE;
    if (fabs(p->c[1]) == INF) return FALSE;
    if (fabs(p->c[2]) == INF) return FALSE;
    if (fabs(p->c[3]) == INF) return FALSE;
    if (fabs(p->c[4]) == INF) return FALSE;
    if (fabs(p->c[5]) == INF) return FALSE;
    return TRUE;
  }

bool_t r6_eq(r6_t *p, r6_t *q)
  { if (p->c[0] != q->c[0]) return FALSE;
    if (p->c[1] != q->c[1]) return FALSE;
    if (p->c[2] != q->c[2]) return FALSE;
    if (p->c[3] != q->c[3]) return FALSE;
    if (p->c[4] != q->c[4]) return FALSE;
    if (p->c[5] != q->c[5]) return FALSE;
    return TRUE;
  }

void r6_throw_cube (r6_t *r)
  { int32_t i;
    for (i = 0; i < N; i++)
      { r->c[i] = 2.0 * drandom() - 1.0; }
  }

void r6_throw_dir (r6_t *r)
  { /* Generate a nonzero Gaussian random vector: */
    int32_t i;
    double r2;
    do
      { r6_throw_normal(r);
        /* Discard if too close to origin: */
        r2 = 0.0;
        for (i = 0; i < N; i++) { double ci = r->c[i]; r2 += ci*ci; }
      }
    while (r2 < 1.0e-20);
    /* Normalize to unit length: */
    double m = sqrt(r2);
    for (i = 0; i < N; i++) { r->c[i] /= m; }
  }

void r6_throw_ball (r6_t *r)
  { rn_throw_ball(N, &(r->c[0])); }

void r6_throw_normal (r6_t *r)
  { int32_t i;
    for (i = 0; i < N; i++)
      { r->c[i] = dgaussrand(); }
  }

void r6_print (FILE *f, r6_t *a)
  { r6_gen_print(f, a, NULL, NULL, NULL, NULL); }

void r6_gen_print (FILE *f, r6_t *a, char *fmt, char *lp, char *sep, char *rp)
  { rn_gen_print(f, N, &(a->c[0]), fmt, lp, sep, rp); }

vec_typeimpl(r6_vec_t,r6_vec,r6_t);
