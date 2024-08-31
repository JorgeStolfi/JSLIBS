/* See r3.h */
/* Last edited on 2024-08-30 17:57:59 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <r3.h>

#include <jsrandom.h>
#include <affirm.h>
#include <sign.h>
#include <rn.h>
#include <vec.h>

#define N 3

void r3_zero(r3_t *r)
  { r->c[0] = 0.0;
    r->c[1] = 0.0;
    r->c[2] = 0.0;
  }

void r3_all(double x, r3_t *r)
  { r->c[0] = x;
    r->c[1] = x;
    r->c[2] = x;
  }

void r3_axis(int32_t i, r3_t *r)
  { affirm((i >= 0) && (i < N), "r3_axis: bad index");
    r->c[0] = 0.0;
    r->c[1] = 0.0;
    r->c[2] = 0.0;

    r->c[i] = 1.0;
  }

void r3_add(r3_t *a, r3_t *b, r3_t *r)
  { r->c[0] = a->c[0] + b->c[0];
    r->c[1] = a->c[1] + b->c[1];
    r->c[2] = a->c[2] + b->c[2];
  }

void r3_sub(r3_t *a, r3_t *b, r3_t *r)
  { r->c[0] = a->c[0] - b->c[0];
    r->c[1] = a->c[1] - b->c[1];
    r->c[2] = a->c[2] - b->c[2];
  }

void r3_neg(r3_t *a, r3_t *r)
  { r->c[0] = - a->c[0];
    r->c[1] = - a->c[1];
    r->c[2] = - a->c[2];
  }

void r3_scale(double s, r3_t *a, r3_t *r)
  { r->c[0] = s * a->c[0];
    r->c[1] = s * a->c[1];
    r->c[2] = s * a->c[2];
  }

void r3_mix(double s, r3_t *a, double t, r3_t *b, r3_t *r)
  { r->c[0] = s * a->c[0] + t * b->c[0];
    r->c[1] = s * a->c[1] + t * b->c[1];
    r->c[2] = s * a->c[2] + t * b->c[2];
  }

void r3_mix_in(double s, r3_t *a, r3_t *r)
  { r->c[0] += s * a->c[0];
    r->c[1] += s * a->c[1];
    r->c[2] += s * a->c[2];
  }

void r3_weigh(r3_t *a, r3_t *w, r3_t *r)
  { r->c[0] = a->c[0] * w->c[0];
    r->c[1] = a->c[1] * w->c[1];
    r->c[2] = a->c[2] * w->c[2];
  }

void r3_unweigh(r3_t *a, r3_t *w, r3_t *r)
  { r->c[0] = a->c[0] / w->c[0];
    r->c[1] = a->c[1] / w->c[1];
    r->c[2] = a->c[2] / w->c[2];
  }

void r3_rot_axis(r3_t *a, int32_t i, int32_t j, double ang, r3_t *r)
  {
    affirm((i >= 0) && (i < N), "r3_rot_axis: bad index {i}");
    affirm((j >= 0) && (j < N), "r3_rot_axis: bad index {j}");
    affirm(i != j, "r3_rot_axis: axes not distinct");
    (*r) = (*a);
    double c = cos(ang);
    double s = sin(ang);
    double x = + c*a->c[i] - s*a->c[j];
    double y = + s*a->c[i] + c*a->c[j];
    r->c[i] = x;
    r->c[j] = y;
  }

double r3_norm(r3_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    return sqrt(a0*a0 + a1*a1 + a2*a2);
  }

double r3_norm_sqr(r3_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    return a0*a0 + a1*a1 + a2*a2;
  }

double r3_L_inf_norm(r3_t *a)
  { double d = 0.0;
    double a0 = fabs(a->c[0]);
    double a1 = fabs(a->c[1]);
    double a2 = fabs(a->c[2]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    if (a2 > d) d = a2;
    return d;
  }

double r3_dist(r3_t *a, r3_t *b)
  { double d0 = (a->c[0] - b->c[0]);
    double d1 = (a->c[1] - b->c[1]);
    double d2 = (a->c[2] - b->c[2]);
    double d = sqrt(d0*d0 + d1*d1 + d2*d2);
    return d;
  }

double r3_dist_sqr(r3_t *a, r3_t *b)
  { double d0 = (a->c[0] - b->c[0]);
    double d1 = (a->c[1] - b->c[1]);
    double d2 = (a->c[2] - b->c[2]);
    return d0*d0 + d1*d1 + d2*d2;
  }

double r3_L_inf_dist(r3_t *a, r3_t *b)
  { double d = 0.0;
    double d0 = fabs(a->c[0] - b->c[0]);
    double d1 = fabs(a->c[1] - b->c[1]);
    double d2 = fabs(a->c[2] - b->c[2]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    if (d2 > d) d = d2;
    return d;
  }

double r3_dir(r3_t *a, r3_t *r)
  { double d = sqrt(a->c[0]*a->c[0] + a->c[1]*a->c[1] + a->c[2]*a->c[2]);
    r->c[0] = a->c[0]/d;
    r->c[1] = a->c[1]/d;
    r->c[2] = a->c[2]/d;
    return d;
  }

double r3_L_inf_dir(r3_t *a, r3_t *r)
  { double d = 0.0;
    double a0 = fabs(a->c[0]);
    double a1 = fabs(a->c[1]);
    double a2 = fabs(a->c[2]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    if (a2 > d) d = a2;
    r->c[0] = a->c[0]/d;
    r->c[1] = a->c[1]/d;
    r->c[2] = a->c[2]/d;
    return d;
  }

double r3_dot(r3_t *a, r3_t *b)
  { return a->c[0]*b->c[0] + a->c[1]*b->c[1] + a->c[2]*b->c[2]; }

double r3_cos(r3_t *a, r3_t *b)
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

double r3_sin(r3_t *a, r3_t *b)
  { return rn_sin(N, &(a->c[0]), &(b->c[0])); }

double r3_angle(r3_t *a, r3_t *b)
  { return rn_angle(N, &(a->c[0]), &(b->c[0])); }

void r3_cross(r3_t *a, r3_t *b, r3_t *r)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];

    double b0 = b->c[0];
    double b1 = b->c[1];
    double b2 = b->c[2];

    r->c[0] = (a1 * b2) - (a2 * b1);
    r->c[1] = (a2 * b0) - (a0 * b2);
    r->c[2] = (a0 * b1) - (a1 * b0);
  }

double r3_det(r3_t *a, r3_t *b, r3_t *c)
  { r3_t p;
    r3_cross(a, b, &p);
    return p.c[0]*c->c[0] + p.c[1]*c->c[1] + p.c[2]*c->c[2];
  }

double r3_decomp(r3_t *a, r3_t *u, r3_t *para, r3_t *perp)
  { double u0 = u->c[0];
    double u1 = u->c[1];
    double u2 = u->c[2];

    double sau = a->c[0]*u0 + a->c[1]*u1 + a->c[2]*u2;
    if (sau == 0.0) 
      { if (para != NULL)
          { para->c[0] = 0.0; 
            para->c[1] = 0.0;
            para->c[2] = 0.0; 
          }

        if (perp != NULL)
          { perp->c[0] = a->c[0]; 
            perp->c[1] = a->c[1];
            perp->c[2] = a->c[2]; 
          } 
        return 0.0;
      }
    else
      { double suu = u0*u0 + u1*u1 + u2*u2;
        double c = sau / suu;
        double p0 = c * u0; 
        double p1 = c * u1; 
        double p2 = c * u2; 

        if (para != NULL)
          { para->c[0] = p0;
            para->c[1] = p1;
            para->c[2] = p2;
          }

        if (perp != NULL)
          { perp->c[0] = a->c[0] - p0; 
            perp->c[1] = a->c[1] - p1;
            perp->c[2] = a->c[2] - p2; 
          }
        return c;
      }
  }

bool_t r3_is_finite(r3_t *p)
  { if (fabs(p->c[0]) == INF) return FALSE;
    if (fabs(p->c[1]) == INF) return FALSE;
    if (fabs(p->c[2]) == INF) return FALSE;
    return TRUE;
  }

bool_t r3_eq(r3_t *p, r3_t *q)
  { if (p->c[0] != q->c[0]) return FALSE;
    if (p->c[1] != q->c[1]) return FALSE;
    if (p->c[2] != q->c[2]) return FALSE;
    return TRUE;
  }
  
void r3_barycenter(int32_t np, r3_t p[], double w[], r3_t *bar)
  { r3_t sum_wp = (r3_t){{ 0, 0, 0 }};
    double sum_w = 0.0;
    for (int32_t k = 0; k < np; k++) 
      { r3_t *pk = &(p[k]);
        double wk = (w != NULL ? w[k] : 1.0);
        r3_mix(1.0, &sum_wp, wk, pk, &sum_wp);
        sum_w += wk;
      }
    r3_scale(1.0/sum_w, &sum_wp, bar);
  }

void r3_bbox(int32_t np, r3_t p[], interval_t B[], bool_t finite)
  { double cmin[N], cmax[N];
    for (int32_t j = 0; j < N; j++) { cmin[j] = +INF; cmax[j] = -INF; }
    int32_t ip;
    for (ip = 0; ip < np; ip++)
      { r3_t *pi = &(p[ip]);
        if ((! finite) || r3_is_finite(pi))
          { for (int32_t j = 0; j < N; j++) 
              { double cij = pi->c[j];
                if (cij < cmin[j]) { cmin[j] = cij; }
                if (cij > cmax[j]) { cmax[j] = cij; }
              }
          }
      }
    for (int32_t j = 0; j < N; j++)
      { B[j] = (interval_t){{ cmin[j], cmax[j] }}; }
  }

sign_t r3_orient(r3_t *a, r3_t *b, r3_t *c, r3_t *d)
  { r3_t ab, ac, ad;
    r3_sub(b, a, &ab);
    r3_sub(c, a, &ac);
    r3_sub(d, a, &ad);
    double det = r3_det(&ab, &ac, &ad);
    if (det > 0)
      { return +1; }
    else if (det < 0)
      { return -1; }
    else
      { return 0; }
  }

r3_t r3_circumcenter(r3_t *a, r3_t *b, r3_t *c, r3_t *d)
  { 
    affirm(FALSE, "** !!! {r3_circumcenter} not implemented yet !!!");
  }

bool_t r3_insphere(r3_t *a, r3_t *b, r3_t *c, r3_t *d, r3_t *e)
  {
    affirm(FALSE, "** !!! {r3_insphere} not implemented yet !!!");
  }

void r3_throw_cube(r3_t *r)
  { r->c[0] = 2.0 * drandom() - 1.0;
    r->c[1] = 2.0 * drandom() - 1.0;
    r->c[2] = 2.0 * drandom() - 1.0;
  }

void r3_throw_dir(r3_t *r)
  { double z = 2.0 * drandom() - 1.0;
    double m = sqrt(1 - z*z);
    double theta = 2*M_PI*drandom();
    r->c[0] = m*cos(theta);
    r->c[1] = m*sin(theta);
    r->c[2] = z;
  }

void r3_throw_ball(r3_t *r)
  { double x, y, z;
    do
      { x = 2.0 * drandom() - 1.0; r->c[0] = x;
        y = 2.0 * drandom() - 1.0; r->c[1] = y;
        z = 2.0 * drandom() - 1.0; r->c[2] = z;
      }
    while (x*x + y*y + z*z >= 1.0);
  }

void r3_throw_normal(r3_t *r)
  { r->c[0] = dgaussrand();
    r->c[1] = dgaussrand();
    r->c[2] = dgaussrand();
  }

void r3_print(FILE *f, r3_t *a)
  { r3_gen_print(f, a, NULL, NULL, NULL, NULL); }

void r3_gen_print(FILE *f, r3_t *a, char *fmt, char *lp, char *sep, char *rp)
  { rn_gen_print(f, N, &(a->c[0]), fmt, lp, sep, rp); }

vec_typeimpl(r3_vec_t, r3_vec, r3_t);
