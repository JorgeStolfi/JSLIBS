/* See r2.h */
/* Last edited on 2020-10-16 01:55:12 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#include <r2.h>

#include <jsrandom.h>
#include <interval.h>
#include <affirm.h>
#include <sign.h>
#include <rn.h>
#include <vec.h>

#define N 2

void r2_zero (r2_t *r)
  { r->c[0] = 0.0;
    r->c[1] = 0.0;
  }

void r2_all (double x, r2_t *r)
  { r->c[0] = x;
    r->c[1] = x;
  }

void r2_axis (int i, r2_t *r)
  { affirm((i >= 0) && (i < N), "r2_axis: bad index");
    r->c[0] = 0.0;
    r->c[1] = 0.0;

    r->c[i] = 1.0;
  }

void r2_add (r2_t *a, r2_t *b, r2_t *r)
  { r->c[0] = a->c[0] + b->c[0];
    r->c[1] = a->c[1] + b->c[1];
  }

void r2_sub (r2_t *a, r2_t *b, r2_t *r)
  { r->c[0] = a->c[0] - b->c[0];
    r->c[1] = a->c[1] - b->c[1];
  }

void r2_neg (r2_t *a, r2_t *r)
  { r->c[0] = - a->c[0];
    r->c[1] = - a->c[1];
  }

void r2_scale (double s, r2_t *a, r2_t *r)
  { r->c[0] = s * a->c[0];
    r->c[1] = s * a->c[1];
  }

void r2_mix (double s, r2_t *a, double t, r2_t *b, r2_t *r)
  { r->c[0] = s * a->c[0] + t * b->c[0];
    r->c[1] = s * a->c[1] + t * b->c[1];
  }

void r2_mix_in (double s, r2_t *a, r2_t *r)
  { r->c[0] += s * a->c[0];
    r->c[1] += s * a->c[1];
  }

void r2_weigh (r2_t *a, r2_t *w, r2_t *r)
  { r->c[0] = a->c[0] * w->c[0];
    r->c[1] = a->c[1] * w->c[1];
  }

double r2_norm (r2_t *a)
  { return hypot(a->c[0], a->c[1]); }

double r2_norm_sqr (r2_t *a)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    return a0*a0 + a1*a1;
  }

double r2_L_inf_norm (r2_t *a)
  { double d = 0.0;
    double a0 = fabs(a->c[0]);
    double a1 = fabs(a->c[1]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    return (d);
  }

double r2_dist (r2_t *a, r2_t *b)
  { double d0 = (a->c[0] - b->c[0]);
    double d1 = (a->c[1] - b->c[1]);
    return hypot(d0, d1);
  }

double r2_dist_sqr (r2_t *a, r2_t *b)
  { double d0 = (a->c[0] - b->c[0]);
    double d1 = (a->c[1] - b->c[1]);
    return d0*d0 + d1*d1;
  }

double r2_L_inf_dist (r2_t *a, r2_t *b)
  { double d = 0.0;
    double d0 = fabs(a->c[0] - b->c[0]);
    double d1 = fabs(a->c[1] - b->c[1]);
    if (d0 > d) d = d0;
    if (d1 > d) d = d1;
    return (d);
  }

double r2_dir (r2_t *a, r2_t *r)
  { double d = sqrt(a->c[0]*a->c[0] + a->c[1]*a->c[1]);
    r->c[0] = a->c[0]/d;
    r->c[1] = a->c[1]/d;
    return (d);
  }

double r2_L_inf_dir (r2_t *a, r2_t *r)
  { double d = 0.0;
    double a0 = fabs(a->c[0]);
    double a1 = fabs(a->c[1]);
    if (a0 > d) d = a0;
    if (a1 > d) d = a1;
    r->c[0] = a->c[0]/d;
    r->c[1] = a->c[1]/d;
    return (d);
  }

double r2_dot (r2_t *a, r2_t *b)
  { return a->c[0]*b->c[0] + a->c[1]*b->c[1]; }

double r2_cos (r2_t *a, r2_t *b)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    double b0 = b->c[0];
    double b1 = b->c[1];
    double ab = a0*b0 + a1*b1;
    double aa = a0*a0 + a1*a1;
    double bb = b0*b0 + b1*b1;
    return ab/(sqrt(aa)*sqrt(bb));
  }

double r2_sin (r2_t *a, r2_t *b)
  { return rn_sin(N, &(a->c[0]), &(b->c[0])); }

double r2_angle (r2_t *a, r2_t *b)
  { return rn_angle(N, &(a->c[0]), &(b->c[0])); }

void r2_cross (r2_t *a, r2_t *r)
  { double a0 = a->c[0];
    double a1 = a->c[1];
    r->c[0] = - a1;
    r->c[1] =   a0;
  }

double r2_det (r2_t *a, r2_t *b)
  { return a->c[0]*b->c[1] - a->c[1]*b->c[0]; }

double r2_decomp (r2_t *a, r2_t *u, r2_t *para, r2_t *perp)
  { double u0 = u->c[0];
    double u1 = u->c[1];
    double sau = a->c[0]*u0 + a->c[1]*u1;
    if (sau == 0.0) 
      { if (para != NULL)
          { para->c[0] = 0.0; 
            para->c[1] = 0.0;
          }
        if (perp != NULL)
          { perp->c[0] = a->c[0]; 
            perp->c[1] = a->c[1];
          }
        return 0.0;
      }
    else
      { double suu = u0*u0 + u1*u1;
        double c = sau / suu;
        double p0 = c*u0; 
        double p1 = c*u1; 

        if (para != NULL)
          { para->c[0] = p0;
            para->c[1] = p1;
          }

        if (perp != NULL)
          { perp->c[0] = a->c[0] - p0; 
            perp->c[1] = a->c[1] - p1;
          }
        return c;
      }
  }

int r2_is_finite(r2_t *p)
  { return (isfinite(p->c[0]) && isfinite(p->c[1]));
  }

int r2_eq(r2_t *p, r2_t *q)
  { return (p->c[0] == q->c[0]) && (p->c[1] == q->c[1]);
  }
  
void r2_barycenter(int32_t np, r2_t p[], double w[], r2_t *bar)
  { r2_t sum_wp = (r2_t){{ 0, 0 }};
    double sum_w = 0.0;
    for (int32_t k = 0; k < np; k++) 
      { r2_t *pk = &(p[k]);
        double wk = (w != NULL ? w[k] : 1.0);
        r2_mix(1.0, &sum_wp, wk, pk, &sum_wp);
        sum_w += wk;
      }
    r2_scale(1.0/sum_w, &sum_wp, bar);
  }

void r2_bbox(int32_t np, r2_t p[], interval_t B[], bool_t finite)
  { double xmin = INF, xmax = -INF;
    double ymin = INF, ymax = -INF;
    int ip;
    for (ip = 0; ip < np; ip++)
      { r2_t *pi = &(p[ip]);
        if ((! finite) || r2_is_finite(pi))
          { double xi = pi->c[0];
            if (xi < xmin) { xmin = xi; }
            if (xi > xmax) { xmax = xi; }
            double yi = pi->c[1];
            if (yi < ymin) { ymin = yi; }
            if (yi > ymax) { ymax = yi; }
          }
      }
    B[0] = (interval_t){{xmin, xmax}};
    B[1] = (interval_t){{ymin, ymax}};
  }

r2_t r2_circumcenter(r2_t *a, r2_t *b, r2_t *c)
  { double xa = a->c[0], ya = a->c[1];
    double xb = b->c[0], yb = b->c[1];
    double xc = c->c[0], yc = c->c[1];
    double A00 = 2*(xa-xb), A01 = 2*(ya-yb), B0 = (xa-xb)*(xa+xb)+(ya-yb)*(ya+yb);
    double A10 = 2*(xa-xc), A11 = 2*(ya-yc), B1 = (xa-xc)*(xa+xc)+(ya-yc)*(ya+yc);
    double det = A00*A11 - A01*A10;
    r2_t v;
    v.c[0] = (B0*A11 - B1*A01)/det;
    v.c[1] = (A00*B1 - A10*B0)/det;
    return v;
  }

sign_t r2_orient(r2_t *a, r2_t *b, r2_t *c)
  { double x1 = a->c[0], y1 = a->c[1];
    double x2 = b->c[0], y2 = b->c[1];
    double x3 = c->c[0], y3 = c->c[1];
    double det = (x2*y3-y2*x3) - (x1*y3-y1*x3) + (x1*y2-y1*x2);
    if (det > 0)
      { return +1; }
    else if (det < 0)
      { return -1; }
    else
      { return 0; }
  }

bool_t r2_incircle(r2_t *a, r2_t *b, r2_t *c, r2_t *d)
  {
    double x1 = a->c[0], y1 = a->c[1];
    double x2 = b->c[0], y2 = b->c[1];
    double x3 = c->c[0], y3 = c->c[1];
    double x4 = d->c[0], y4 = d->c[1];
    
    double da = ((y4-y1)*(x2-x3)+(x4-x1)*(y2-y3))*((x4-x3)*(x2-x1)-(y4-y3)*(y2-y1));
    double db = ((y4-y3)*(x2-x1)+(x4-x3)*(y2-y1))*((x4-x1)*(x2-x3)-(y4-y1)*(y2-y3));

    return da > db;
  }

void r2_throw_cube (r2_t *r)
  { r->c[0] = 2.0 * drandom() - 1.0;
    r->c[1] = 2.0 * drandom() - 1.0;
  }

void r2_throw_dir (r2_t *r)
  { double theta = 2*M_PI*drandom();
    r->c[0] = cos(theta);
    r->c[1] = sin(theta);
  }

void r2_throw_ball (r2_t *r)
  { double x, y;
    do
      { x = 2.0 * drandom() - 1.0; r->c[0] = x;
        y = 2.0 * drandom() - 1.0; r->c[1] = y;
      }
    while (x*x + y*y >= 1.0);
  }

void r2_throw_normal (r2_t *r)
  { r->c[0] = dgaussrand();
    r->c[1] = dgaussrand();
  }

void r2_print (FILE *f, r2_t *a)
  { r2_gen_print(f, a, NULL, NULL, NULL, NULL); }

void r2_gen_print (FILE *f, r2_t *a, char *fmt, char *lp, char *sep, char *rp)
  { rn_gen_print(f, N, &(a->c[0]), fmt, lp, sep, rp); }

vec_typeimpl(r2_vec_t,r2_vec,r2_t);
