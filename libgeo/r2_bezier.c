/* See r2_bezier.h */
/* Last edited on 2023-10-02 08:43:45 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <affirm.h>

#include <r2_bezier.h>

double r2_bezier_length_estimate(r2_t *p0, r2_t *p1, r2_t *p2, r2_t *p3, int32_t order)
  { 
    if (order <= 0)
      { return r2_dist(p0,p3); }
    else
      { r2_t p01, p12, p23, p012, p123, p0123; /* Intermediate points. */
        r2_bezier_split(0.0, 1.0, p0, p1, p2, p3, 0.5, &p01, &p12, &p23, &p012, &p123, &p0123, NULL);      
        double d0 = r2_bezier_length_estimate(p0, &p01, &p012, &p0123, order-1);
        double d1 = r2_bezier_length_estimate(&p0123, &p123, &p23, p3, order-1);
        return d0 + d1;
      }
  }

void r2_bezier_split
  ( double t0, 
    double t1,
    r2_t *p0, 
    r2_t *p1, 
    r2_t *p2, 
    r2_t *p3, 
    double t,
    r2_t *p01, 
    r2_t *p12,
    r2_t *p23, 
    r2_t *p012, 
    r2_t *p123,
    r2_t *p0123,
    r2_t *v
  )
  {
    /* The DeCasteljau algorithm: */
    double dt = t1 - t0;
    demand(dt != 0, "path duration is zero");
    double b = (t - t0)/dt;
    double a = 1 - b;
    
    r2_mix(a, p0, b, p1, p01);
    r2_mix(a, p1, b, p2, p12);
    r2_mix(a, p2, b, p3, p23);

    r2_mix(a, p01, b, p12, p012);
    r2_mix(a, p12, b, p23, p123);

    r2_mix(a, p012, b, p123, p0123);

    if (v != NULL) 
      { r2_sub(p123, p012, v);
        r2_scale(3/dt, v, v);
      }
  }

void r2_bezier_eval
  ( double t0, 
    double t1,
    r2_t *p0, 
    r2_t *p1, 
    r2_t *p2, 
    r2_t *p3, 
    double t,
    r2_t *p,
    r2_t *v
  )
  {
    r2_t p01, p12, p23, p012, p123;
    r2_bezier_split(t0, t1, p0, p1, p2, p3, t, &p01, &p12, &p23, &p012, &p123, p, v);
  }

void r2_bezier_from_bend(r2_t *p0, r2_t *p3, double bend, r2_t *p1, r2_t *p2)
  { 
    /* The straight edge's coord system: */
    r2_t o; r2_mix(0.5, p0, 0.5, p3, &o); /* Midpoint of {org,dst}. */
    r2_t u; 
    r2_sub(p3, p0, &u);
    double h = r2_dir(&u, &u)/2;
    demand (fabs(h) > 0.001, "edge too short");
    r2_t v = (r2_t){{ -u.c[1], +u.c[0] }};
    (void)r2_dir(&v, &v);
    
    double cx = (2*bend*bend + h*h)/(3*h);
    double cy = 4*bend/3;
    
    (*p1) = (r2_t){{ - cx*u.c[0] + cy*v.c[0], - cx*u.c[1] + cy*v.c[1] }};
    r2_add(&o, p1, p1);
    (*p2) = (r2_t){{ + cx*u.c[0] + cy*v.c[0], + cx*u.c[1] + cy*v.c[1] }};
    r2_add(&o, p2, p2);
  }
