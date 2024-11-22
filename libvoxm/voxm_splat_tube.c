/* See voxm_splat_tube.h */
/* Last edited on 2024-11-11 07:37:51 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <r3x3.h>
#include <affirm.h>
#include <ppv_array.h>

#include <voxm_obj.h>
#include <r3_path.h>
#include <r3_bezier.h>
#include <voxm_splat.h>

#include <voxm_splat_tube.h>

void voxm_splat_tube_round_helix
  ( ppv_array_t *A, 
    double t0,
    double t1,
    r3_motion_state_t *S, 
    double len, 
    double ang,
    double hht, 
    double inR, 
    double otR, 
    double fuzzR, 
    bool_t sub,
    r3_motion_state_t *S0,
    r3_motion_state_t *S1
  )
  {
    bool_t debug = TRUE;
    
    demand((otR >= inR) && (inR >= 0), "invalid tube radii");
    demand(t1 >= t0, "invalid time interval");

    auto r3_motion_state_t helix_path(double t);
      /* Function that computes the position and orientation at time {t} along
        the helix defined by {len,ang,hht}, modified by state {S}. */
    
    auto double fuzzy_donut(r3_t *p);
      /* Indicator function for a donut of major radius {(iR+oR)/2},
        minor radius {(oR-iR)/2}, with a fuzzy layer of thickness
        {2*fuzzR}. The donut is rotationally symmetric around
        the {X} axis and symmetric across the {YZ} plane. */
       
    auto double fuzzy_ball(r3_t *p);
      /* Indicator function for a ball of radius {(iR+oR)/2},
        minor radius {(oR-iR)/2}, with a fuzzy layer of thickness
        {2*fuzzR}. */
       
    /* Estimate the total arc length {totL} of outer edge from {t0} to {t1}: */
    double totL = (t1 - t0)*(hypot(len,hht) + ang*otR);
    
    /* Select the number of samples {ns} : */
    double step = 0.25;     /* Linear step length for path enumeration. */
    int32_t ns = 1 + (int32_t)ceil(totL/step); /* Number of sample points to use. */
    if (debug) { fprintf(stderr, "  sampling path at %d points\n", ns); }
    
    /* Sample the path states {sp[0..ns-1]}: */
    r3_motion_state_t *sp = notnull(malloc(ns*sizeof(r3_motion_state_t)), "no mem");
    bool_t mids = FALSE;
    r3_motion_sample_uniform(helix_path, t0, t1, ns, mids, NULL, sp);
    if (S0 != NULL) { (*S0) = sp[0]; }
    if (S1 != NULL) { (*S1) = sp[ns-1]; }
    
    if (!sub)
      { /* Splat the tube: */
        voxm_splat_object_multi(A, fuzzy_donut, ns, sp, otR + fuzzR, FALSE);
      }
    else
      { /* Subtract the bore: */
        voxm_splat_object_multi(A, fuzzy_ball, ns, sp, inR + fuzzR, TRUE);
      }
      
    if (debug) { fprintf(stderr, "\n"); }
    return;
    
    /* INTERNAL IMPLEMENTATION */
    
    r3_motion_state_t helix_path(double t)
      { 
        r3_motion_state_t T;
        r3_motion_helix(t, len, ang, hht, &T);
        r3_motion_state_compose(&T, S, &T);
        return T;
      }
      
    double fuzzy_donut(r3_t *p)
      { /* Major and minor donut radii: */
        double minR = (otR - inR)/2;
        double majR = (otR + inR)/2;
        return voxm_obj_donut(p, minR, majR, 0, fuzzR);
      }
      
    double fuzzy_ball(r3_t *p)
      { return voxm_obj_ball(p, inR, fuzzR); }
  }

void voxm_splat_tube_round_bezier
  ( ppv_array_t *A, 
    r3_t *p0, 
    r3_t *p1, 
    r3_t *p2, 
    r3_t *p3, 
    double inR, 
    double otR, 
    double fuzzR, 
    bool_t sub
  )
  {
    bool_t debug = TRUE;
    
    demand((otR >= inR) && (inR >= 0), "invalid tube radii");

    auto r3_motion_state_t bezier_path(double t);
      /* Function that computes the state {S} (position and orientation)
        at time {t} along the Bezier arc defined by {p0,p1,p2,p3}.
        The arc is assumed to be at {p0} when {t=0} and at {p3} when {t=1}. 
        The matrix {S.M} is orthonormal; 
        the vector {S.u} will be tangent to the curve, while {S.v} and {S.w} will
        be arbitrary vectors orthogonal to each other and to {S.u}. */
    
    auto double fuzzy_donut(r3_t *p);
      /* Indicator function for a donut of major radius {(inR+otR)/2},
        minor radius {(otR-inR)/2}, with a fuzzy layer of thickness
        {2*fuzzR}. The donut is rotationally symmetric around
        the {X} axis and mirror-symmetric across the {YZ} plane. */
       
    auto double fuzzy_ball(r3_t *p);
      /* Indicator function for a donut of major radius {(inR+otR)/2},
        minor radius {(otR-inR)/2}, with a fuzzy layer of thickness
        {2*fuzzR}. The donut is rotationally symmetric around
        the {Z} axis and symmetric across the {XY} plane. */
       
    /* Estimate the arc length {len} and select the number of samples {ns} : */
    double step = 0.25;     /* Linear step length for path enumeration. */
    double len = r3_dist(p0,p1) + r3_dist(p1,p2) + r3_dist(p2,p3);
    int32_t ns = 1 + (int32_t)ceil(len/step); /* Number of sample points to use. */
    if (debug) { fprintf(stderr, "sampling path at %d points\n", ns); }
    
    /* Sample the path times and states {st[0..ns-1],sp[0..ns-1]}: */
    double *t = notnull(malloc(ns*sizeof(double)), "no mem");
    r3_motion_state_t *S = notnull(malloc(ns*sizeof(r3_motion_state_t)), "no mem");
    bool_t mids = FALSE;
    r3_motion_sample_uniform(bezier_path, 0.0, 1.0, ns, mids, t, S);
    
    if (!sub)
      { /* Splat the tube: */
        voxm_splat_object_multi(A, fuzzy_donut, ns, S, otR + fuzzR, FALSE);
      }
    else
      { /* Subtract the bore: */
        voxm_splat_object_multi(A, fuzzy_ball, ns, S, otR + fuzzR, TRUE);
      }
      
    if (debug) { fprintf(stderr, "\n"); }
    return;
    
    /* INTERNAL IMPLEMENTATION */
    
    r3_motion_state_t bezier_path(double t)
      { 
        /* Compute the position and velocity of the Bezier arc at {t}: */
        r3_t p01, p12, p23, p012, p123; /* Intermediate points. */
        r3_path_state_t S; /* Position and velocity of the Bezier arc at time {t}. */
        S.t = t;
        r3_bezier_split(0.0, 1.0, p0, p1, p2, p3, t, &p01, &p12, &p23, &p012, &p123, &(S.p), &(S.v));      
        
        /* Convert to a path state: */
        r3_motion_state_t T = r3_path_state_to_r3_motion_state(&S);
        return T;
      }
      
    double fuzzy_donut(r3_t *p)
      { /* Major and minor donut radii: */
        double minR = (otR - inR)/2;
        double majR = (otR + inR)/2;
        return voxm_obj_donut(p, minR, majR, 0, fuzzR);
      }
      
    double fuzzy_ball(r3_t *p)
      { return voxm_obj_ball(p, inR, fuzzR); }
  }

void voxm_splat_tube_round_segment
  ( ppv_array_t *A, 
    r3_path_state_t *S, 
    r3_path_state_t *T, 
    double inR, 
    double otR, 
    double fuzzR, 
    bool_t sub
  )
  {
    /* Temporary: convert to Bezier arc. */
    /* Compute Bézier points: */
    r3_t *p0 = &(S->p);
    r3_t p1, p2;
    r3_t *p3 = &(T->p);
    r3_path_bezier_from_states(S, T, &p1, &p2);
    voxm_splat_tube_round_bezier(A, p0, &p1, &p2, p3, inR, otR, fuzzR, sub);
  }
  
