#ifndef ellipse_aligned_H
#define ellipse_aligned_H

/* ellipse_aligned.h -- tools for axis-aligned, origin-centered ellipses. */
/* Last edited on 2021-06-09 19:55:19 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <sign.h>
#include <bool.h>
#include <interval.h>

/* ELLIPSE GEOMETRY

  This moduel has tools for axis-aligned ellipses centered at the origin.
  
  Such an ellipse is defined by the X and Y semidiameters (X and Y radii)
  {rx,ry}.

  The ellipse is defined by the parametric equation {G(th) =
  (rx*cos(th), ry*sin(th))}, where {th} is the /angular argument/ of the
  ellipse, varying in {[0:2*PI]}.

  Alternatively, a point {p=(x,y)} is inside, on the boundary of, or
  outside the ellipse depending on whether {G(x,y) = (x/rx)^2 + (y/ry)^2
  - 1} is negative, zero, or positive.

  The /canonical coordinates/ of a point {p=(x,y)} are the 
  quantities {(x/rx,y/ry)}.  In the canonical coordinates,
  the ellipse is the unit origin-centered circle.
  
  
  !!! Currently we require {rx >= ry >= 0}; remove this restriction! !!!
  */

void ellipse_aligned_bbox(double rx, double ry, interval_t bbox[]);
  /* Stores in {bbox[0..1]} a bounding box for the ellipse {G}. */

void ellipse_aligned_int_bbox
  ( double cx, double cy,  /* Center of ellipse. */
    double rx, double ry,
    double mrg, /* Extra margin. */
    int32_t *xLoP,  /* (OUT) Min X of clip area. */
    int32_t *xHiP,  /* (OUT) Max X of clip area. */
    int32_t *yLoP,  /* (OUT) Min Y of clip area. */
    int32_t *yHiP   /* (OUT) Max Y of clip area. */
  );
  /* Stores in {*xLoP,*xHiP,*yLoP,*yHiP} a bounding box with integer coordinates
    for the ellipse with center {cx,xy} and X,Y radii {rx,ry}. The box is at least
    {mrg} units away from the ellipse, on each side. */

bool_t ellipse_aligned_inside(double rx, double ry, double xp, double yp);
  /* Returns TRUE if {p} is inside the ellipse centered at the origin
    with radii {rx,ry}, FALSE if outside. Beware of roundoff errors. */

double ellipse_aligned_position(double rx, double ry, double xp, double yp, double *cosP, double *sinP);
  /* Returns the radial and angular position of {p=(xp,yp)} relative to the
    ellipse centered at the origin with radii {rx,ry}. 
    
    The returned result is the position of {p} relative to the
    ellipse's boundary, along the ray from the origin to {p}.
    This value is less than 1 if {p} is inside, and greater than 1 if {p}
    is outside. It is 0 if {p} is at the origin, and 1 if it is on the
    boundary. 
    
    If {cosP,sinP} are not NULL, the procedure also stores through them
    the cosine and sine of the angular argument of {p}, measured from the
    X axis. */

double ellipse_aligned_nearest_point(double rx, double ry, double xp, double yp, double *xqP, double *yqP);
  /* Fids the signed distance from {p=(xp,yp)} to the boundary of the ellipse
    centered at the origin with radii {rx,ry}.
    
    If {xqP,yqP} is not NULL, the procedure also stores throgh them the coordinates
    of the points on the boundary of the ellipse with radii {rx,ry} that is closest to {p}. 
    
    The return value is {+dist(p,q)} if {p} is outside, and {-dist(p,q)}
    if {p} is inside. */

double ellipse_aligned_border_position(double rx, double ry, double hwd, double xp, double yp);
  /* If the distance from {p} to the boundary of the ellipse centered at the origin 
    with radii {rx,ry} is less than {hwd}, returns the signed distance
    (positive outside, negative inside) divided by {hwd}. Otherwise
    returns {+1} if {p} is outside the ellipse, {-1} if it is inside.
    
    This is the relative transversal position of {p} in the strip
    created when the ellipse's boundary is stroked with a round brush
    of radius {hwd}. This procedure is much faster than
    {ellipse_aligned_nearest_point} when {p} is far from the
    stroked region. */

void ellipse_aligned_print(FILE *wr, double rx, double ry, char *fmt);
  /* Writes to {wr} the ellipse radii {rx,ry}.  Each parameter
    is printed with format {fmt} (which should be appropriate 
    for a {double} value). */

/* SPECIALIZED TOOLS */

double ellipse_aligned_compute_t(double A, double B);
  /* Solves the polynomial equation {P(t) == 0} that occurs in the
    nearest-point computation, where
    
      { P(t) = A*(1+t**2)*(1-t**2) - 2*t*((B+1)*t**2 - (B-1)) }
      
    Requires {A >= 0, B >= 0}.  The result is in {[0 _ 1]}. */

#endif
 
