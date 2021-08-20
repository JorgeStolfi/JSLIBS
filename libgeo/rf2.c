/* See rf2.h. */
/* Last edited on 2021-08-20 16:15:07 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <rf2.h>

#include <jsrandom.h>
#include <affirm.h>
#include <rfn.h>
#include <vec.h>

#define N 2

rf2_t rf2_zero (void)
  { rf2_t r = (rf2_t){{ 0, 0 }};
    return r;
  }

rf2_t rf2_all (float v)
  { rf2_t r = (rf2_t){{ v, v }};
    return r;
  }

rf2_t rf2_axis (int32_t i)
  { affirm((i >= 0) && (i < N), "rf2_axis: bad index");
    rf2_t r = (rf2_t){{ 0, 0 }};
    r.c[i] = 1.0;
    return r;
  }

rf2_t rf2_add (rf2_t* const a, rf2_t* const b)
  { rf2_t r;
    r.c[0] = (float)(((double)a->c[0]) + b->c[0]);
    r.c[1] = (float)(((double)a->c[1]) + b->c[1]);
    return r;
  }

rf2_t rf2_sub (rf2_t* const a, rf2_t* const b)
  { rf2_t r;
    r.c[0] = (float)(((double)a->c[0]) - b->c[0]);
    r.c[1] = (float)(((double)a->c[1]) - b->c[1]);
    return r;
  }

rf2_t rf2_neg (rf2_t* const a)
  { rf2_t r;
    r.c[0] = - a->c[0];
    r.c[1] = - a->c[1];
    return r;
  }

rf2_t rf2_scale (double s, rf2_t* const a)
  { rf2_t r;
    r.c[0] = (float)(s * a->c[0]);
    r.c[1] = (float)(s * a->c[1]);
    return r;
  }

rf2_t rf2_mix (double s, rf2_t* const a, double t, rf2_t* const b)
  { rf2_t r;
    r.c[0] = (float)(s * a->c[0] + t * b->c[0]);
    r.c[1] = (float)(s * a->c[1] + t * b->c[1]);
    return r;
  }

void rf2_mix_in (double s, rf2_t* const a, rf2_t* const r)
  { r->c[0] = (float)(r->c[0] + s * a->c[0]);
    r->c[1] = (float)(r->c[1] + s * a->c[1]);
  }

rf2_t rf2_weigh (rf2_t* const a, rf2_t* const w)
  { rf2_t r;
    r.c[0] = (float)(((double)a->c[0]) * w->c[0]);
    r.c[1] = (float)(((double)a->c[1]) * w->c[1]);
    return r;
  }

rf2_t rf2_unweigh (rf2_t* const a, rf2_t* const w)
  { rf2_t r;
    r.c[0] = (float)(((double)a->c[0]) / w->c[0]);
    r.c[1] = (float)(((double)a->c[1]) / w->c[1]);
    return r;
  }

rf2_t rf2_rot (rf2_t* const a, double ang)
  { double c = cos(ang);
    double s = sin(ang);
    double x = + c*a->c[0] - s*a->c[1];
    double y = + s*a->c[0] + c*a->c[1];
    rf2_t r = (rf2_t){{ (float)x, (float)y }};
    return r;
  }

double rf2_norm (rf2_t* const a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    return sqrt(a0*a0 + a1*a1);
  }

double rf2_norm_sqr (rf2_t* const a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    return a0*a0 + a1*a1;
  }

float rf2_L_inf_norm (rf2_t* const a)
  { float d = 0.0;
    float a0 = fabsf(a->c[0]);
    float a1 = fabsf(a->c[1]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    return d;
  }

double rf2_dist (rf2_t* const a, rf2_t* const b)
  { double d0 = ((double)a->c[0]) - b->c[0];
    double d1 = ((double)a->c[1]) - b->c[1];
    double d = sqrt(d0*d0 + d1*d1);
    return d;
  }

double rf2_dist_sqr (rf2_t* const a, rf2_t* const b)
  { double d0 = ((double)a->c[0]) - b->c[0];
    double d1 = ((double)a->c[1]) - b->c[1];
    return d0*d0 + d1*d1;
  }

double rf2_L_inf_dist (rf2_t* const a, rf2_t* const b)
  { double d = 0.0;
    double d0 = fabs(((double)a->c[0]) - b->c[0]);
    double d1 = fabs(((double)a->c[1]) - b->c[1]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    return d;
  }

rf2_t rf2_dir (rf2_t* const a, double *normP)
  { double d = rf2_norm(a);
    rf2_t r;
    r.c[0] = (float)(a->c[0]/d);
    r.c[1] = (float)(a->c[1]/d);
    if (normP != NULL) { *normP = d; }
    return r;
    
  }

rf2_t rf2_L_inf_dir (rf2_t* const a, float *normP)
  { float d = rf2_L_inf_norm(a);
    rf2_t r;
    r.c[0] = (float)(((double)a->c[0])/d);
    r.c[1] = (float)(((double)a->c[1])/d);
    if (normP != NULL) { *normP = d; }
    return r;
  }

double rf2_dot (rf2_t* const a, rf2_t* const b)
  { double d0 = ((double)a->c[0])*b->c[0];
    double d1 = ((double)a->c[1])*b->c[1];
    return  d0 + d1; 
  }

double rf2_cos (rf2_t* const a, rf2_t* const b)
  { double a0 = a->c[0];
    double a1 = a->c[1];

    double b0 = b->c[0];
    double b1 = b->c[1];

    double ab = a0*b0 + a1*b1;
    double aa = a0*a0 + a1*a1;
    double bb = b0*b0 + b1*b1;
    return ab/(sqrt(aa)*sqrt(bb));
  }

double rf2_sin (rf2_t* const a, rf2_t* const b)
  { return rfn_sin(N, &(a->c[0]), &(b->c[0])); }

double rf2_angle (rf2_t* const a, rf2_t* const b)
  { return rfn_angle(N, &(a->c[0]), &(b->c[0])); }

rf2_t rf2_cross (rf2_t* const a)
  { rf2_t r;
    r.c[0] = - a->c[1];
    r.c[1] = + a->c[0];
    return r;
  }

double rf2_det (rf2_t* const a, rf2_t* const b)
  { 
    return ((double)a->c[0])*b->c[1] - ((double)a->c[1])*b->c[0];
  }

double rf2_decomp (rf2_t* const a, rf2_t* const u, rf2_t* para, rf2_t* perp)
  { float *paran = (para == NULL ? NULL : &(para->c[0]));
    float *perpn = (perp == NULL ? NULL : &(perp->c[0]));
    return rfn_decomp(N, &(a->c[0]), &(u->c[0]), paran, perpn);
  }

bool_t rf2_is_finite(rf2_t* const p)
  { if (fabsf(p->c[0]) == INF) return FALSE;
    if (fabsf(p->c[1]) == INF) return FALSE;
    return TRUE;
  }

bool_t rf2_eq(rf2_t* const p, rf2_t* const q)
  { if (p->c[0] != q->c[0]) return FALSE;
    if (p->c[1] != q->c[1]) return FALSE;
    return TRUE;
  }

rf2_t rf2_throw_cube (void)
  { rf2_t r;
    for (int32_t i = 0; i < N; i++)
      { r.c[i] = (float)(2.0 * frandom() - 1.0); }
    return r;
  }

rf2_t rf2_throw_dir (void)
  { /* Generate a nonzero Gaussian random vector: */
    rf2_t r;
    double r2;
    do
      { r = rf2_throw_normal();
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

rf2_t rf2_throw_ball (void)
  { rf2_t r;
    rfn_throw_ball(N, &(r.c[0]));
    return r;
  }

rf2_t rf2_throw_normal (void)
  { rf2_t r;
    for (int32_t i = 0; i < N; i++)
      { r.c[i] = fgaussrand(); }
    return r;
  }

void rf2_print (FILE *f, rf2_t* const a)
  { rf2_gen_print(f, a, NULL, NULL, NULL, NULL); }

void rf2_gen_print (FILE *f, rf2_t* const a, char *fmt, char *lp, char *sep, char *rp)
  { rfn_gen_print(f, N, &(a->c[0]), fmt, lp, sep, rp); }

vec_typeimpl(rf2_vec_t,rf2_vec,rf2_t);
