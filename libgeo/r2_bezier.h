/* r2_bezier.h --- Bezier arcs in {\RR^3 */
/* Last edited on 2024-12-05 10:27:47 by stolfi */

#ifndef r2_bezier_H
#define r2_bezier_H

#include <stdint.h>

#include <bool.h>
#include <r2.h>

double r2_bezier_length_estimate(r2_t *p0, r2_t *p1, r2_t *p2, r2_t *p3, int32_t order);
  /* Estimates the length of the cubic curve arc defined by the Bezier 
     control points {p0,p1,p2,p3}, by bisecting it recursively {order} times. */

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
  );
  /* Evaluates a Bezier arc defined by {*p0,*p1,*p2,*p3} 
    at time {t}, assumed to vary between {t0} and {t1}.
    Returns the result in {*p}.  If {v} is not {NULL}, returns in {*v} the velocity of the 
    curve (the time derivative of {*p}) at time {t}. */
    
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
  );
  /* Subdivides a Bezier arc defined by {*p0,*p1,*p2,*p3} 
    into two arcs at time {t}, assumed to vary between {t0} and {t1}.
    The control points of the first arc are {*p0,*p01,*p012,*p0123}; those
    of the second arc are {*p0123,*p123,*p23,*p3}.  
    If {v} is not {NULL}, returns in {*v} the velocity of the 
    curve (the time derivative of {*p0123}) at time {t}. */

void r2_bezier_from_bend(r2_t *p0, r2_t *p3, double bend, r2_t *p1, r2_t *p2);
  /* Computes the Bézier control points {p1} and {p2} for an approximately circular
    arc that goes from {p0} to (p3} and deviates aboyt {bend} from the straight line. */

#endif

