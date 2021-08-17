/* See rf4.h. */
/* Last edited on 2021-08-17 05:12:12 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <rf4.h>

#include <jsrandom.h>
#include <affirm.h>
#include <rfn.h>
#include <vec.h>

#define N 4

rf4_t rf4_axis (int32_t i)
  { affirm((i >= 0) && (i < N), "rf4_axis: bad index");
    rf4_t r;
    r.c[0] = 0.0;
    r.c[1] = 0.0;
    r.c[2] = 0.0;
    r.c[3] = 0.0;

    r.c[i] = 1.0;
    return r;
  }

rf4_t rf4_add (rf4_t *a, rf4_t *b)
  { rf4_t r;
    r.c[0] = a->c[0] + b->c[0];
    r.c[1] = a->c[1] + b->c[1];
    r.c[2] = a->c[2] + b->c[2];
    r.c[3] = a->c[3] + b->c[3];
    return r;
  }

rf4_t rf4_sub (rf4_t *a, rf4_t *b)
  { rf4_t r;
    r.c[0] = a->c[0] - b->c[0];
    r.c[1] = a->c[1] - b->c[1];
    r.c[2] = a->c[2] - b->c[2];
    r.c[3] = a->c[3] - b->c[3];
    return r;
  }

rf4_t rf4_neg (rf4_t *a)
  { rf4_t r;
    r.c[0] = - a->c[0];
    r.c[1] = - a->c[1];
    r.c[2] = - a->c[2];
    r.c[3] = - a->c[3];
    return r;
  }

rf4_t rf4_scale (double s, rf4_t *a)
  { rf4_t r;
    r.c[0] = (float)(s * a->c[0]);
    r.c[1] = (float)(s * a->c[1]);
    r.c[2] = (float)(s * a->c[2]);
    r.c[3] = (float)(s * a->c[3]);
    return r;
  }

rf4_t rf4_mix (double s, rf4_t *a, double t, rf4_t *b)
  { rf4_t r;
    r.c[0] = (float)(s * a->c[0] + t * b->c[0]);
    r.c[1] = (float)(s * a->c[1] + t * b->c[1]);
    r.c[2] = (float)(s * a->c[2] + t * b->c[2]);
    r.c[3] = (float)(s * a->c[3] + t * b->c[3]);
    return r;
  }

void rf4_mix_in (double s, rf4_t *a, rf4_t *r)
  { r->c[0] = (float)(r->c[0] + s * a->c[0]);
    r->c[1] = (float)(r->c[1] + s * a->c[1]);
    r->c[2] = (float)(r->c[2] + s * a->c[2]);
    r->c[3] = (float)(r->c[3] + s * a->c[3]);
  }

rf4_t rf4_weigh (rf4_t *a, rf4_t *w)
  { rf4_t r;
    r.c[0] = (float)(a->c[0] * w->c[0]);
    r.c[1] = (float)(a->c[1] * w->c[1]);
    r.c[2] = (float)(a->c[2] * w->c[2]);
    r.c[3] = (float)(a->c[3] * w->c[3]);
    return r;
  }

double rf4_norm (rf4_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    double a3 = a->c[3];
    return sqrt(a0*a0 + a1*a1 + a2*a2 + a3*a3);
  }

double rf4_norm_sqr (rf4_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    double a3 = a->c[3];
    return a0*a0 + a1*a1 + a2*a2 + a3*a3;
  }

float rf4_L_inf_norm (rf4_t *a)
  { float d = 0.0;
    float a0 = fabsf(a->c[0]);
    float a1 = fabsf(a->c[1]);
    float a2 = fabsf(a->c[2]);
    float a3 = fabsf(a->c[3]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    if (a2 > d) d = a2;
    if (a3 > d) d = a3;
    return d;
  }

double rf4_dist (rf4_t *a, rf4_t *b)
  { double d0 = ((double)(a->c[0]) - b->c[0]);
    double d1 = ((double)(a->c[1]) - b->c[1]);
    double d2 = ((double)(a->c[2]) - b->c[2]);
    double d3 = ((double)(a->c[3]) - b->c[3]);
    double d = sqrt(d0*d0 + d1*d1 + d2*d2 + d3*d3);
    return (d);
  }

double rf4_dist_sqr (rf4_t *a, rf4_t *b)
  { double d0 = ((double)(a->c[0]) - b->c[0]);
    double d1 = ((double)(a->c[1]) - b->c[1]);
    double d2 = ((double)(a->c[2]) - b->c[2]);
    double d3 = ((double)(a->c[3]) - b->c[3]);
    return d0*d0 + d1*d1 + d2*d2 + d3*d3;
  }

double rf4_L_inf_dist (rf4_t *a, rf4_t *b)
  { double d = 0.0;
    double d0 = fabs((double)(a->c[0]) - b->c[0]);
    double d1 = fabs((double)(a->c[1]) - b->c[1]);
    double d2 = fabs((double)(a->c[2]) - b->c[2]);
    double d3 = fabs((double)(a->c[3]) - b->c[3]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    if (d2 > d) d = d2;
    if (d3 > d) d = d3;
    return d;
  }

rf4_t rf4_dir (rf4_t *a, double *normP)
  { double d = rf4_norm(a);
    rf4_t r;
    r.c[0] = (float)(a->c[0]/d);
    r.c[1] = (float)(a->c[1]/d);
    r.c[2] = (float)(a->c[2]/d);
    r.c[3] = (float)(a->c[3]/d);
    if (normP != NULL) { *normP = d; }
    return r;
    
  }

rf4_t rf4_L_inf_dir (rf4_t *a, float *normP)
  { float d = rf4_L_inf_norm(a);
    rf4_t r;
    r.c[0] = (float)((double)(a->c[0])/d);
    r.c[1] = (float)((double)(a->c[1])/d);
    r.c[2] = (float)((double)(a->c[2])/d);
    r.c[3] = (float)((double)(a->c[3])/d);
    if (normP != NULL) { *normP = d; }
    return r;
  }

double rf4_dot (rf4_t *a, rf4_t *b)
  { double d0 = (double)(a->c[0])*b->c[0];
    double d1 = (double)(a->c[1])*b->c[1];
    double d2 = (double)(a->c[2])*b->c[2];
    double d3 = (double)(a->c[3])*b->c[3];
    return  d0 + d1 + d2 + d3; 
  }

double rf4_cos (rf4_t *a, rf4_t *b)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    double a3 = a->c[3];

    double b0 = b->c[0];
    double b1 = b->c[1];
    double b2 = b->c[2];
    double b3 = b->c[3];

    double ab = a0*b0 + a1*b1 + a2*b2 + a3*b3;
    double aa = a0*a0 + a1*a1 + a2*a2 + a3*a3;
    double bb = b0*b0 + b1*b1 + b2*b2 + b3*b3;
    return ab/(sqrt(aa)*sqrt(bb));
  }

double rf4_sin (rf4_t *a, rf4_t *b)
  { return rfn_sin(N, &(a->c[0]), &(b->c[0])); }

double rf4_angle (rf4_t *a, rf4_t *b)
  { return rfn_angle(N, &(a->c[0]), &(b->c[0])); }

rf4_t rf4_cross (rf4_t *a, rf4_t *b, rf4_t *c)
  { double d01 = (double)(a->c[0])*b->c[1] - (double)(a->c[1])*b->c[0];
    double d02 = (double)(a->c[0])*b->c[2] - (double)(a->c[2])*b->c[0];
    double d12 = (double)(a->c[1])*b->c[2] - (double)(a->c[2])*b->c[1];
    double d03 = (double)(a->c[0])*b->c[3] - (double)(a->c[3])*b->c[0];
    double d13 = (double)(a->c[1])*b->c[3] - (double)(a->c[3])*b->c[1];
    double d23 = (double)(a->c[2])*b->c[3] - (double)(a->c[3])*b->c[2];

    rf4_t r;
    r.c[0] = (float)(- d12*c->c[3] + d13*c->c[2] - d23*c->c[1]);
    r.c[1] = (float)(+ d02*c->c[3] - d03*c->c[2] + d23*c->c[0]);
    r.c[2] = (float)(- d01*c->c[3] + d03*c->c[1] - d13*c->c[0]);
    r.c[3] = (float)(+ d01*c->c[2] - d02*c->c[1] + d12*c->c[0]);
    return r;
  }

double rf4_det (rf4_t *a, rf4_t *b, rf4_t *c, rf4_t *d)
  { double d01 = (double)(a->c[0])*b->c[1] - (double)(a->c[1])*b->c[0];
    double d02 = (double)(a->c[0])*b->c[2] - (double)(a->c[2])*b->c[0];
    double d12 = (double)(a->c[1])*b->c[2] - (double)(a->c[2])*b->c[1];
    double d03 = (double)(a->c[0])*b->c[3] - (double)(a->c[3])*b->c[0];
    double d13 = (double)(a->c[1])*b->c[3] - (double)(a->c[3])*b->c[1];
    double d23 = (double)(a->c[2])*b->c[3] - (double)(a->c[3])*b->c[2];

    double r0 = - d12*c->c[3] + d13*c->c[2] - d23*c->c[1];
    double r1 = + d02*c->c[3] - d03*c->c[2] + d23*c->c[0];
    double r2 = - d01*c->c[3] + d03*c->c[1] - d13*c->c[0];
    double r3 = + d01*c->c[2] - d02*c->c[1] + d12*c->c[0];
  
    return r0*d->c[0] + r1*d->c[1] + r2*d->c[2] + r3*d->c[3];
  }

double rf4_decomp (rf4_t *a, rf4_t *u, rf4_t *para, rf4_t *perp)
  { float *paran = (para == NULL ? NULL : &(para->c[0]));
    float *perpn = (perp == NULL ? NULL : &(perp->c[0]));
    return rfn_decomp(N, &(a->c[0]), &(u->c[0]), paran, perpn);
  }

bool_t rf4_is_finite(rf4_t *p)
  { if (fabsf(p->c[0]) == INF) return FALSE;
    if (fabsf(p->c[1]) == INF) return FALSE;
    if (fabsf(p->c[2]) == INF) return FALSE;
    if (fabsf(p->c[3]) == INF) return FALSE;
    return TRUE;
  }

bool_t rf4_eq(rf4_t *p, rf4_t *q)
  { if (p->c[0] != q->c[0]) return FALSE;
    if (p->c[1] != q->c[1]) return FALSE;
    if (p->c[2] != q->c[2]) return FALSE;
    if (p->c[3] != q->c[3]) return FALSE;
    return TRUE;
  }

rf4_t rf4_throw_cube (void)
  { rf4_t r;
    for (int32_t i = 0; i < N; i++)
      { r.c[i] = (float)(2.0 * frandom() - 1.0); }
    return r;
  }

rf4_t rf4_throw_dir (void)
  { /* Generate a nonzero Gaussian random vector: */
    rf4_t r;
    double r2;
    do
      { r = rf4_throw_normal();
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

rf4_t rf4_throw_ball (void)
  { rf4_t r;
    rfn_throw_ball(N, &(r.c[0]));
    return r;
  }

rf4_t rf4_throw_normal (void)
  { rf4_t r;
    for (int32_t i = 0; i < N; i++)
      { r.c[i] = (float)(fgaussrand()); }
    return r;
  }

void rf4_print (FILE *f, rf4_t *a)
  { rf4_gen_print(f, a, NULL, NULL, NULL, NULL); }

void rf4_gen_print (FILE *f, rf4_t *a, char *fmt, char *lp, char *sep, char *rp)
  { rfn_gen_print(f, N, &(a->c[0]), fmt, lp, sep, rp); }

vec_typeimpl(rf4_vec_t,rf4_vec,rf4_t);
