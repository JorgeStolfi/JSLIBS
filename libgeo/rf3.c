/* See rf3.h. */
/* Last edited on 2021-08-17 08:53:45 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <rf3.h>

#include <jsrandom.h>
#include <affirm.h>
#include <rfn.h>
#include <vec.h>

#define N 3

rf3_t rf3_axis (int32_t i)
  { affirm((i >= 0) && (i < N), "rf3_axis: bad index");
    rf3_t r;
    r.c[0] = 0.0;
    r.c[1] = 0.0;
    r.c[2] = 0.0;

    r.c[i] = 1.0;
    return r;
  }

rf3_t rf3_add (rf3_t *a, rf3_t *b)
  { rf3_t r;
    r.c[0] = a->c[0] + b->c[0];
    r.c[1] = a->c[1] + b->c[1];
    r.c[2] = a->c[2] + b->c[2];
    return r;
  }

rf3_t rf3_sub (rf3_t *a, rf3_t *b)
  { rf3_t r;
    r.c[0] = a->c[0] - b->c[0];
    r.c[1] = a->c[1] - b->c[1];
    r.c[2] = a->c[2] - b->c[2];
    return r;
  }

rf3_t rf3_neg (rf3_t *a)
  { rf3_t r;
    r.c[0] = - a->c[0];
    r.c[1] = - a->c[1];
    r.c[2] = - a->c[2];
    return r;
  }

rf3_t rf3_scale (double s, rf3_t *a)
  { rf3_t r;
    r.c[0] = (float)(s * a->c[0]);
    r.c[1] = (float)(s * a->c[1]);
    r.c[2] = (float)(s * a->c[2]);
    return r;
  }

rf3_t rf3_mix (double s, rf3_t *a, double t, rf3_t *b)
  { rf3_t r;
    r.c[0] = (float)(s * a->c[0] + t * b->c[0]);
    r.c[1] = (float)(s * a->c[1] + t * b->c[1]);
    r.c[2] = (float)(s * a->c[2] + t * b->c[2]);
    return r;
  }

void rf3_mix_in (double s, rf3_t *a, rf3_t *r)
  { r->c[0] = (float)(r->c[0] + s * a->c[0]);
    r->c[1] = (float)(r->c[1] + s * a->c[1]);
    r->c[2] = (float)(r->c[2] + s * a->c[2]);
  }

rf3_t rf3_weigh (rf3_t *a, rf3_t *w)
  { rf3_t r;
    r.c[0] = (float)(a->c[0] * w->c[0]);
    r.c[1] = (float)(a->c[1] * w->c[1]);
    r.c[2] = (float)(a->c[2] * w->c[2]);
    return r;
  }

double rf3_norm (rf3_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    return sqrt(a0*a0 + a1*a1 + a2*a2);
  }

double rf3_norm_sqr (rf3_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    return a0*a0 + a1*a1 + a2*a2;
  }

float rf3_L_inf_norm (rf3_t *a)
  { float d = 0.0;
    float a0 = fabsf(a->c[0]);
    float a1 = fabsf(a->c[1]);
    float a2 = fabsf(a->c[2]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    if (a2 > d) d = a2;
    return d;
  }

double rf3_dist (rf3_t *a, rf3_t *b)
  { double d0 = ((double)(a->c[0]) - b->c[0]);
    double d1 = ((double)(a->c[1]) - b->c[1]);
    double d2 = ((double)(a->c[2]) - b->c[2]);
    double d = sqrt(d0*d0 + d1*d1 + d2*d2);
    return d;
  }

double rf3_dist_sqr (rf3_t *a, rf3_t *b)
  { double d0 = ((double)(a->c[0]) - b->c[0]);
    double d1 = ((double)(a->c[1]) - b->c[1]);
    double d2 = ((double)(a->c[2]) - b->c[2]);
    return d0*d0 + d1*d1 + d2*d2;
  }

double rf3_L_inf_dist (rf3_t *a, rf3_t *b)
  { double d = 0.0;
    double d0 = fabs((double)(a->c[0]) - b->c[0]);
    double d1 = fabs((double)(a->c[1]) - b->c[1]);
    double d2 = fabs((double)(a->c[2]) - b->c[2]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    if (d2 > d) d = d2;
    return d;
  }

rf3_t rf3_dir (rf3_t *a, double *normP)
  { double d = rf3_norm(a);
    rf3_t r;
    r.c[0] = (float)(a->c[0]/d);
    r.c[1] = (float)(a->c[1]/d);
    r.c[2] = (float)(a->c[2]/d);
    if (normP != NULL) { *normP = d; }
    return r;
    
  }

rf3_t rf3_L_inf_dir (rf3_t *a, float *normP)
  { float d = rf3_L_inf_norm(a);
    rf3_t r;
    r.c[0] = (float)((double)(a->c[0])/d);
    r.c[1] = (float)((double)(a->c[1])/d);
    r.c[2] = (float)((double)(a->c[2])/d);
    if (normP != NULL) { *normP = d; }
    return r;
  }

double rf3_dot (rf3_t *a, rf3_t *b)
  { double d0 = (double)(a->c[0])*b->c[0];
    double d1 = (double)(a->c[1])*b->c[1];
    double d2 = (double)(a->c[2])*b->c[2];
    return  d0 + d1 + d2; 
  }

double rf3_cos (rf3_t *a, rf3_t *b)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];

    double b0 = b->c[0];
    double b1 = b->c[1];
    double b2 = b->c[2];

    double ab = a0*b0 + a1*b1 + a2*b2;
    double aa = a0*a0 + a1*a1 + a2*a2;
    double bb = b0*b0 + b1*b1 + b2*b2;
    return ab/(sqrt(aa)*sqrt(bb));
  }

double rf3_sin (rf3_t *a, rf3_t *b)
  { return rfn_sin(N, &(a->c[0]), &(b->c[0])); }

double rf3_angle (rf3_t *a, rf3_t *b)
  { return rfn_angle(N, &(a->c[0]), &(b->c[0])); }

rf3_t rf3_cross (rf3_t *a, rf3_t *b)
  { double d01 = (double)(a->c[0])*b->c[1] - (double)(a->c[1])*b->c[0];
    double d02 = (double)(a->c[0])*b->c[2] - (double)(a->c[2])*b->c[0];
    double d12 = (double)(a->c[1])*b->c[2] - (double)(a->c[2])*b->c[1];

    rf3_t r;
    r.c[0] = + (float)(d12);
    r.c[1] = - (float)(d02);
    r.c[2] = + (float)(d01);
    return r;
  }

double rf3_det (rf3_t *a, rf3_t *b, rf3_t *c)
  { double d01 = (double)(a->c[0])*b->c[1] - (double)(a->c[1])*b->c[0];
    double d02 = (double)(a->c[0])*b->c[2] - (double)(a->c[2])*b->c[0];
    double d12 = (double)(a->c[1])*b->c[2] - (double)(a->c[2])*b->c[1];
    return d12*c->c[0] - d02*c->c[1] + d01*c->c[2];
  }

double rf3_decomp (rf3_t *a, rf3_t *u, rf3_t *para, rf3_t *perp)
  { float *paran = (para == NULL ? NULL : &(para->c[0]));
    float *perpn = (perp == NULL ? NULL : &(perp->c[0]));
    return rfn_decomp(N, &(a->c[0]), &(u->c[0]), paran, perpn);
  }

bool_t rf3_is_finite(rf3_t *p)
  { if (fabsf(p->c[0]) == INF) return FALSE;
    if (fabsf(p->c[1]) == INF) return FALSE;
    if (fabsf(p->c[2]) == INF) return FALSE;
    return TRUE;
  }

bool_t rf3_eq(rf3_t *p, rf3_t *q)
  { if (p->c[0] != q->c[0]) return FALSE;
    if (p->c[1] != q->c[1]) return FALSE;
    if (p->c[2] != q->c[2]) return FALSE;
    return TRUE;
  }

rf3_t rf3_throw_cube (void)
  { rf3_t r;
    for (int32_t i = 0; i < N; i++)
      { r.c[i] = (float)(2.0 * frandom() - 1.0); }
    return r;
  }

rf3_t rf3_throw_dir (void)
  { /* Generate a nonzero Gaussian random vector: */
    rf3_t r;
    double r2;
    do
      { r = rf3_throw_normal();
        /* Discard if too close to origin: */
        r2 = 0.0;
        for (int32_t i = 0; i < N; i++) { double ci = r.c[i]; r2 += ci*ci; }
      }
    while (r2 < 1.0e-10);
    /* Normalize to unit length: */
    double m = sqrt(r2);
    for (int32_t i = 0; i < N; i++) { r.c[i] = (float)(r.c[i]/m); }
    return r;
  }

rf3_t rf3_throw_ball (void)
  { rf3_t r;
    rfn_throw_ball(N, &(r.c[0]));
    return r;
  }

rf3_t rf3_throw_normal (void)
  { rf3_t r;
    for (int32_t i = 0; i < N; i++)
      { r.c[i] = (float)(fgaussrand()); }
    return r;
  }

void rf3_print (FILE *f, rf3_t *a)
  { rf3_gen_print(f, a, NULL, NULL, NULL, NULL); }

void rf3_gen_print (FILE *f, rf3_t *a, char *fmt, char *lp, char *sep, char *rp)
  { rfn_gen_print(f, N, &(a->c[0]), fmt, lp, sep, rp); }

vec_typeimpl(rf3_vec_t,rf3_vec,rf3_t);
