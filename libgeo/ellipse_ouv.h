#ifndef ellipse_ouv_H
#define ellipse_ouv_H

/* ellipse_ouv.h -- tools for ellipses in the dir-radius form. */
/* Last edited on 2021-06-09 20:14:44 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <r2.h>
#include <sign.h>
#include <interval.h>

/* ELLIPSE GEOMETRY */

typedef struct ellipse_ouv_t 
  { r2_t u;     /* Unit direction of long axis (major semidiameter). */
    double a;   /* Major semi-diameter. */
    r2_t v;     /* Unit direction of short axis (minor semidiameter). */
    double b;   /* Minor semi-diameter. */
  } ellipse_ouv_t;
  /* Defines an ellipse centered at the origin. The ellipse has major
    semidiameter {a*u}, minor semidiameter {b*v}, with {u,v}
    orthogonal and unitary, and {a >= b >= 0}.

    The vectors {u,v} and the origin {(0,0)} define the ellipse's 
    /UV system/.  The /UV coordinates/ of {p} are {p\cdot u,p\cdot v}. 
    
    In the UV system, the ellipse is defined by the parametric equation
    {F(th) = (a*cos(th), b*sin(th))}, where {th} is the /angular
    argument/ of the ellipse, varying in {[0:2*PI]}.
    
    Alternatively, a point with UV coordinates {up,vp} is inside, on
    the boundary of, or outside the ellipse depending on whether
    {F(up,vp) = (up/a)^2 + (vp/b)^2 - 1} is negative, zero, or
    positive.
    
    Consider doing all your computations in the ellipse's UV
    coordinate system and using the operations in {ellipse_aligned.h},
    which are faster than the ones below. */

void ellipse_ouv_bbox(ellipse_ouv_t *F, interval_t bbox[]);
  /* Stores in {bbox[0..1]} an axis-aligned
    bounding box for the ellipse {F}. */

void ellipse_ouv_int_bbox
  ( r2_t *ctr,  /* Center of ellipse. */
    ellipse_ouv_t *F,
    double mrg, /* Extra margin. */
    int32_t *xLoP,  /* (OUT) Min X of clip area. */
    int32_t *xHiP,  /* (OUT) Max X of clip area. */
    int32_t *yLoP,  /* (OUT) Min Y of clip area. */
    int32_t *yHiP   /* (OUT) Max Y of clip area. */
  );
  /* Stores in {*xLoP,*xHiP,*yLoP,*yHiP} a bounding box with integer coordinates for the
    ellipse {F} displaced by {ctr}. The box is at least {mrg}
    units away from the ellipse, on each side. */

bool_t ellipse_ouv_inside(ellipse_ouv_t *F, r2_t *p);
  /* Returns TRUE if {p} is inside {F}, FALSE if outside.  
    May return either value if {p} is very close to the boundary. */

double ellipse_ouv_position(ellipse_ouv_t *F, r2_t *p, r2_t *csp);
  /* Returns the radial and angular position of {p} relative to the
    ellipse {F}. 
    
    The returned result is the position of {p} relative to the
    ellipse's boundary, along the ray from the origin to {p}.
    This value is less than 1 if {p} is inside, and greater than 1 if {p}
    is outside. It is 0 if {p} is at the origin, and 1 if it is on the
    boundary. 
    
    If {csp} is not NULL, the procedure also stores in {*csp} the
    cosine and sine of the angular argument of {p}, measured from the
    vector {F.u}. */

double ellipse_ouv_nearest_point(ellipse_ouv_t *F, r2_t *p, r2_t *q);
  /* Finds the signed distance from {p} to the boundary of the ellipse {F}.
    If {q} is not NULL, the procedure also stores in {*q} the point on
    the boundary of the ellipse {F} that is closest to {p}. 
    
    The return value {+dist(p,q)} if {p} is outside, and {-dist(p,q)}
    if {p} is inside. */

double ellipse_ouv_border_position(ellipse_ouv_t *F, double hwd, r2_t *p);
  /* If the distance from {p} to {F} is less than {hwd}, returns that
    distance divided by {hwd}. Otherwise returns {+1} if {p} is
    outside the ellipse, {-1} if it is inside.
    
    This is the relative position of {p} in the boundary of {F}, when
    it is painted with a round brush o radius {hwd}. This procedure is
    much faster than {ellipse_ouv_nearest_point} when {p} is not close
    to the stroked region. */

double ellipse_ouv_compute_t(double A, double B);
  /* Solves the polynomial equation {P(t) == 0} that occurs in the
    nearest-point computation, where
    
      { P(t) = A*(1+t**2)*(1-t**2) - 2*t*((B+1)*t**2 - (B-1)) }
      
    Requires {A >= 0, B >= 0}.  The result is in {[0 _ 1]}. */

void ellipse_ouv_print(FILE *wr, ellipse_ouv_t *F, char *fmt);
  /* Writes to {wr} the geometric parameters of {F}.  Each parameter
    is printed with format {fmt} (which should be appropriate 
    for a {double} value). */

#endif
 
