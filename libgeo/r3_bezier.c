/* See r3_bezier.h */
/* Last edited on 2023-10-01 19:22:30 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r3.h>
#include <affirm.h>

#include <r3_bezier.h>

double r3_bezier_length_estimate(r3_t *p0, r3_t *p1, r3_t *p2, r3_t *p3, int32_t order)
  { 
    if (order <= 0)
      { return r3_dist(p0,p3); }
    else
      { 
        r3_t p01, p12, p23, p012, p123, p0123; /* Intermediate points. */
        r3_bezier_split(0.0, 1.0, p0, p1, p2, p3, 0.5, &p01, &p12, &p23, &p012, &p123, &p0123, NULL);      
        double d0 = r3_bezier_length_estimate(p0, &p01, &p012, &p0123, order-1);
        double d1 = r3_bezier_length_estimate(&p0123, &p123, &p23, p3, order-1);
        return d0 + d1;
      }
  }

void r3_bezier_split
  ( double t0, 
    double t1,
    r3_t *p0, 
    r3_t *p1, 
    r3_t *p2, 
    r3_t *p3, 
    double t,
    r3_t *p01, 
    r3_t *p12,
    r3_t *p23, 
    r3_t *p012, 
    r3_t *p123,
    r3_t *p0123,
    r3_t *v
  )
  {
    /* The DeCasteljau algorithm: */
    double dt = t1 - t0;
    demand(dt != 0, "path duration is zero");
    double b = (t - t0)/dt;
    double a = 1 - b;
    
    r3_mix(a, p0, b, p1, p01);
    r3_mix(a, p1, b, p2, p12);
    r3_mix(a, p2, b, p3, p23);

    r3_mix(a, p01, b, p12, p012);
    r3_mix(a, p12, b, p23, p123);

    r3_mix(a, p012, b, p123, p0123);

    if (v != NULL) 
      { r3_sub(p123, p012, v);
        r3_scale(3/dt, v, v);
      }
  }

void r3_bezier_eval
  ( double t0, 
    double t1,
    r3_t *p0, 
    r3_t *p1, 
    r3_t *p2, 
    r3_t *p3, 
    double t,
    r3_t *p,
    r3_t *v
  )
  {
    r3_t p01, p12, p23, p012, p123;
    r3_bezier_split(t0, t1, p0, p1, p2, p3, t, &p01, &p12, &p23, &p012, &p123, p, v);
  }
