/* r3_bezier.h --- voxel-based modeling of antialiased 3D objects */
/* Last edited on 2021-06-09 19:57:20 by jstolfi */

#ifndef r3_bezier_H
#define r3_bezier_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <r3_path.h>

void r3_bezier_from_path_states(r3_path_state_t *S, r3_path_state_t *T, r3_t *p1, r3_t *p2);
  /* Conmputes the two middle Bezier control points of a cubic arc that starts 
    at time {S.t} and position {S.p} with velocity {S.v} and ends at time {T.t} 
    and position {T.p} with velocity {T.v}. The control points
    of the arc will be {p0=S.p}, {p1}, {p2}, and {p3=T.p}. */

double r3_bezier_length_estimate(r3_t *p0, r3_t *p1, r3_t *p2, r3_t *p3, int32_t order);
  /* Estimates the length of the cubic curve arc defined by the Bezier 
     control points {p0,p1,p2,p3}, by bisecting it recursively {order} times. */

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
  );
  /* Subdivides a Bezier arc defined by {*p0,*p1,*p2,*p3} 
    into two arcs at time {t}, assumed to vary between {t0} and {t1}.
    The control points of the first arc are {*p0,*p01,*p012,*p0123}; those
    of the second arc are {*p0123,*p123,*p23,*p3}.  
    If {v} is not {NULL}, returns in {*v} the velocity of the 
    curve (the time derivative of {*p0123}) at time {t}. */

#endif

