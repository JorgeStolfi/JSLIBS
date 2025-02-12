/* See r4.h. */
/* Last edited on 2025-02-05 15:45:44 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <r4.h>

#include <jsrandom.h>
#include <affirm.h>
#include <rn.h>
#include <vec.h>

#define N 4

void r4_zero (r4_t *r)
  { r->c[0] = 0.0;
    r->c[1] = 0.0;
    r->c[2] = 0.0;
    r->c[3] = 0.0;
  }

void r4_all (double x, r4_t *r)
  { r->c[0] = x;
    r->c[1] = x;
    r->c[2] = x;
    r->c[3] = x;
  }

void r4_axis (uint32_t i, r4_t *r)
  { affirm((i >= 0) && (i < N), "r4_axis: bad index");
    r->c[0] = 0.0;
    r->c[1] = 0.0;
    r->c[2] = 0.0;
    r->c[3] = 0.0;

    r->c[i] = 1.0;
  }

void r4_add (r4_t *a, r4_t *b, r4_t *r)
  { r->c[0] = a->c[0] + b->c[0];
    r->c[1] = a->c[1] + b->c[1];
    r->c[2] = a->c[2] + b->c[2];
    r->c[3] = a->c[3] + b->c[3];
  }

void r4_sub (r4_t *a, r4_t *b, r4_t *r)
  { r->c[0] = a->c[0] - b->c[0];
    r->c[1] = a->c[1] - b->c[1];
    r->c[2] = a->c[2] - b->c[2];
    r->c[3] = a->c[3] - b->c[3];
  }

void r4_neg (r4_t *a, r4_t *r)
  { r->c[0] = - a->c[0];
    r->c[1] = - a->c[1];
    r->c[2] = - a->c[2];
    r->c[3] = - a->c[3];
  }

void r4_scale (double s, r4_t *a, r4_t *r)
  { r->c[0] = s * a->c[0];
    r->c[1] = s * a->c[1];
    r->c[2] = s * a->c[2];
    r->c[3] = s * a->c[3];
  }

void r4_mix (double s, r4_t *a, double t, r4_t *b, r4_t *r)
  { r->c[0] = s * a->c[0] + t * b->c[0];
    r->c[1] = s * a->c[1] + t * b->c[1];
    r->c[2] = s * a->c[2] + t * b->c[2];
    r->c[3] = s * a->c[3] + t * b->c[3];
  }

void r4_mix_in (double s, r4_t *a, r4_t *r)
  { r->c[0] += s * a->c[0];
    r->c[1] += s * a->c[1];
    r->c[2] += s * a->c[2];
    r->c[3] += s * a->c[3];
  }

void r4_weigh (r4_t *a, r4_t *w, r4_t *r)
  { r->c[0] = a->c[0] * w->c[0];
    r->c[1] = a->c[1] * w->c[1];
    r->c[2] = a->c[2] * w->c[2];
    r->c[3] = a->c[3] * w->c[3];
  }

void r4_unweigh (r4_t *a, r4_t *w, r4_t *r)
  { r->c[0] = a->c[0] / w->c[0];
    r->c[1] = a->c[1] / w->c[1];
    r->c[2] = a->c[2] / w->c[2];
    r->c[3] = a->c[3] / w->c[3];
  }

void r4_rot_axis (r4_t *a, uint32_t i, uint32_t j, double ang, r4_t *r)
  {
    affirm((i >= 0) && (i < N), "r4_rot_axis: bad index {i}");
    affirm((j >= 0) && (j < N), "r4_rot_axis: bad index {j}");
    affirm(i != j, "r4_rot_axis: axes not distinct");
    (*r) = (*a);
    double c = cos(ang);
    double s = sin(ang);
    double x = + c*a->c[i] - s*a->c[j];
    double y = + s*a->c[i] + c*a->c[j];
    r->c[i] = x;
    r->c[j] = y;
  }

double r4_norm (r4_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    double a3 = a->c[3];
    return sqrt(a0*a0 + a1*a1 + a2*a2 + a3*a3);
  }

double r4_norm_sqr (r4_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    double a3 = a->c[3];
    return a0*a0 + a1*a1 + a2*a2 + a3*a3;
  }

double r4_L_inf_norm (r4_t *a)
  { double d = 0.0;
    double a0 = fabs(a->c[0]);
    double a1 = fabs(a->c[1]);
    double a2 = fabs(a->c[2]);
    double a3 = fabs(a->c[3]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    if (a2 > d) d = a2;
    if (a3 > d) d = a3;
    return d;
  }

double r4_dist (r4_t *a, r4_t *b)
  { double d0 = (a->c[0] - b->c[0]);
    double d1 = (a->c[1] - b->c[1]);
    double d2 = (a->c[2] - b->c[2]);
    double d3 = (a->c[3] - b->c[3]);
    double d = sqrt(d0*d0 + d1*d1 + d2*d2 + d3*d3);
    return d;
  }

double r4_dist_sqr (r4_t *a, r4_t *b)
  { double d0 = (a->c[0] - b->c[0]);
    double d1 = (a->c[1] - b->c[1]);
    double d2 = (a->c[2] - b->c[2]);
    double d3 = (a->c[3] - b->c[3]);
    return d0*d0 + d1*d1 + d2*d2 + d3*d3;
  }

double r4_L_inf_dist (r4_t *a, r4_t *b)
  { double d = 0.0;
    double d0 = fabs(a->c[0] - b->c[0]);
    double d1 = fabs(a->c[1] - b->c[1]);
    double d2 = fabs(a->c[2] - b->c[2]);
    double d3 = fabs(a->c[3] - b->c[3]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    if (d2 > d) d = d2;
    if (d3 > d) d = d3;
    return d;
  }

double r4_dir (r4_t *a, r4_t *r)
  { double d = sqrt
      ( a->c[0]*a->c[0] + 
        a->c[1]*a->c[1] + 
        a->c[2]*a->c[2] + 
        a->c[3]*a->c[3]
      );
    if (d == 0)
      { (*r) = (r4_t){{ NAN, NAN, NAN, NAN }}; }
    else
      { r->c[0] = a->c[0]/d;
        r->c[1] = a->c[1]/d;
        r->c[2] = a->c[2]/d;
        r->c[3] = a->c[3]/d;
      }
    return d;
  }

double r4_L_inf_dir (r4_t *a, r4_t *r)
  { double d = 0.0;
    double a0 = fabs(a->c[0]); if (a0 > d) d = a0;
    double a1 = fabs(a->c[1]); if (a1 > d) d = a1;
    double a2 = fabs(a->c[2]); if (a2 > d) d = a2;
    double a3 = fabs(a->c[3]); if (a3 > d) d = a3;
    if (d == 0)
      { (*r) = (r4_t){{ NAN, NAN, NAN, NAN }}; }
    else
      { r->c[0] = a->c[0]/d;
        r->c[1] = a->c[1]/d;
        r->c[2] = a->c[2]/d;
        r->c[3] = a->c[3]/d;
      }
    return d;
  }

double r4_dot (r4_t *a, r4_t *b)
  { return a->c[0]*b->c[0] + a->c[1]*b->c[1] + a->c[2]*b->c[2] + a->c[3]*b->c[3]; }

double r4_cos (r4_t *a, r4_t *b)
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

double r4_sin (r4_t *a, r4_t *b)
  { return rn_sin(N, &(a->c[0]), &(b->c[0])); }

double r4_angle (r4_t *a, r4_t *b)
  { return rn_angle(N, &(a->c[0]), &(b->c[0])); }

void r4_cross (r4_t *a, r4_t *b, r4_t *c, r4_t *r)
  { double d01 = a->c[0]*b->c[1] - a->c[1]*b->c[0];
    double d02 = a->c[0]*b->c[2] - a->c[2]*b->c[0];
    double d12 = a->c[1]*b->c[2] - a->c[2]*b->c[1];
    double d03 = a->c[0]*b->c[3] - a->c[3]*b->c[0];
    double d13 = a->c[1]*b->c[3] - a->c[3]*b->c[1];
    double d23 = a->c[2]*b->c[3] - a->c[3]*b->c[2];

    r->c[0] = - d12*c->c[3] + d13*c->c[2] - d23*c->c[1];
    r->c[1] =   d02*c->c[3] - d03*c->c[2] + d23*c->c[0];
    r->c[2] = - d01*c->c[3] + d03*c->c[1] - d13*c->c[0];
    r->c[3] =   d01*c->c[2] - d02*c->c[1] + d12*c->c[0];
  }

double r4_det (r4_t *a, r4_t *b, r4_t *c, r4_t *d)
  { r4_t p;
    r4_cross(a, b, c, &p);
    return p.c[0]*d->c[0] + p.c[1]*d->c[1] + p.c[2]*d->c[2] + p.c[3]*d->c[3];
  }

double r4_decomp (r4_t *a, r4_t *u, r4_t *para, r4_t *perp)
  { double *paran = (para == NULL ? NULL : &(para->c[0]));
    double *perpn = (perp == NULL ? NULL : &(perp->c[0]));
    return rn_decomp(N, &(a->c[0]), &(u->c[0]), paran, perpn);
  }

bool_t r4_is_finite(r4_t *p)
  { if (fabs(p->c[0]) == INF) return FALSE;
    if (fabs(p->c[1]) == INF) return FALSE;
    if (fabs(p->c[2]) == INF) return FALSE;
    if (fabs(p->c[3]) == INF) return FALSE;
    return TRUE;
  }

bool_t r4_eq(r4_t *p, r4_t *q)
  { if (p->c[0] != q->c[0]) return FALSE;
    if (p->c[1] != q->c[1]) return FALSE;
    if (p->c[2] != q->c[2]) return FALSE;
    if (p->c[3] != q->c[3]) return FALSE;
    return TRUE;
  }
void r4_barycenter(uint32_t np, r4_t p[], double w[], r4_t *bar)
  { r4_t sum_wp = (r4_t){{ 0, 0, 0, 0 }};
    double sum_w = 0.0;
    for (uint32_t k = 0;  k < np; k++) 
      { r4_t *pk = &(p[k]);
        double wk = (w != NULL ? w[k] : 1.0);
        r4_mix(1.0, &sum_wp, wk, pk, &sum_wp);
        sum_w += wk;
      }
    r4_scale(1.0/sum_w, &sum_wp, bar);
  }

void r4_bbox(uint32_t np, r4_t p[], interval_t B[], bool_t finite)
  { double cmin[N], cmax[N];
    for (uint32_t j = 0;  j < N; j++) { cmin[j] = +INF; cmax[j] = -INF; }
    for (uint32_t ip = 0;  ip < np; ip++)
      { r4_t *pi = &(p[ip]);
        if ((! finite) || r4_is_finite(pi))
          { for (uint32_t j = 0;  j < N; j++) 
              { double cij = pi->c[j];
                if (cij < cmin[j]) { cmin[j] = cij; }
                if (cij > cmax[j]) { cmax[j] = cij; }
              }
          }
      }
    for (uint32_t j = 0;  j < N; j++)
      { B[j] = (interval_t){{ cmin[j], cmax[j] }}; }
  }

void r4_throw_cube (r4_t *r)
  { for (uint32_t i = 0;  i < N; i++)
      { r->c[i] = 2.0 * drandom() - 1.0; }
  }

void r4_throw_dir (r4_t *r)
  { /* Generate a nonzero Gaussian random vector: */
    double r2;
    do
      { r4_throw_normal(r);
        /* Discard if too close to origin: */
        r2 = 0.0;
        for (uint32_t i = 0;  i < N; i++) { double ci = r->c[i]; r2 += ci*ci; }
      }
    while (r2 < 1.0e-20);
    /* Normalize to unit length: */
    double m = sqrt(r2);
    for (uint32_t i = 0;  i < N; i++) { r->c[i] /= m; }
  }

void r4_throw_ball (r4_t *r)
  { rn_throw_ball(N, &(r->c[0])); }

void r4_throw_normal (r4_t *r)
  { for (uint32_t i = 0;  i < N; i++)
      { r->c[i] = dgaussrand(); }
  }
  
double r4_throw_ortho(r4_t *u, r4_t *r)
  { r4_t d = (*u);
    double umag = r4_dir(&d, &d);
    if (umag < 1.0e-320)
      { r4_zero(r); }
    else
      { /* Get an {r} that is not zero nor too close to {±d}: */
        double dvdot; 
        do 
          { r4_throw_ball(r);
            double vnorm = r4_dir(r, r);
            if (vnorm < 0.05) { continue; }
            dvdot = r4_dot(&d, r);
          }
        while (fabs(dvdot) > 0.8);
        /* Make {r} orthogonal to {d} and unit-norm: */
        r4_mix_in(-dvdot, &d, r);
        double vnorm = r4_dir(r, r);
        assert(vnorm > 0.1); /* Paranoia. */
        /* Once more, just to be sure: */
        dvdot = r4_dot(&d, r);
        r4_mix_in(-dvdot, &d, r);
        vnorm = r4_dir(r, r);
        assert(vnorm > 0.999); /* Paranoia. */
        /* Rescale {r} with original length of {u}: */
        if (fabs(umag - 1) > 1.0e-14) { r4_scale(umag, r, r); }
      }
    return umag;
  }

void r4_print (FILE *f, r4_t *a)
  { r4_gen_print(f, a, NULL, NULL, NULL, NULL); }

void r4_gen_print (FILE *f, r4_t *a, char *fmt, char *lp, char *sep, char *rp)
  { rn_gen_print(f, N, &(a->c[0]), fmt, lp, sep, rp); }

vec_typeimpl(r4_vec_t, r4_vec, r4_t);
