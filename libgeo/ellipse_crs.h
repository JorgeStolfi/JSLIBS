#ifndef ellipse_crs_H
#define ellipse_crs_H

/* ellipse_crs.h -- tools for ellipses in the center-radius-stretch form. */
/* Last edited on 2010-03-18 22:11:39 by stolfi */

#include <r2.h>
#include <interval.h>
#include <ellipse_ouv.h>

/* ELLIPSE GEOMETRY */

typedef struct ellipse_crs_t /* Parameters of a sphere's image: */
  { r2_t ctr;    /* Center (pixels). */
    double rad;  /* Radius (length of shortest semidiameter). */
    r2_t str;    /* Stretch vector (longest semidiameter minus {rad}). */
  } ellipse_crs_t;
  /* An {ellipse_crs_t} describes an ellipse on the plane. The length
    of the smallest radius (minor semidiameter) is {rad}. The longest
    radius (major semidiameter) is parallel to the vector {str} and
    has length {rad+len}, where {len} is the Euclidean length of
    {str}.
    
    If you expect to do many operations from this module with the same
    {ellipse_crs_t}, consider converting it to an [ellipse_ouv_t} and
    using the analogous operations in {ellipse_ouv.h}, which are
    generally faster.
    
    Consider also doing all your computations in the ellipse's UV
    coordinate system (whose axes are the major and minor diameters of
    the ellipse) and using the operations in {ellipse_aligned.h},
    which are faster still. */

void ellipse_crs_bbox(ellipse_crs_t *E, interval_t bbox[]);
  /* Stores in {bbox[0..1]} a bounding box for the ellipse {E}. */

void ellipse_crs_int_bbox
  ( ellipse_crs_t *E, 
    double mrg, /* Extra margin. */
    int *xLoP,  /* (OUT) Min X of clip area. */
    int *xHiP,  /* (OUT) Max X of clip area. */
    int *yLoP,  /* (OUT) Min Y of clip area. */
    int *yHiP   /* (OUT) Max Y of clip area. */
  );
  /* Stores in {*xLoP,*xHiP,*yLoP,*yHiP} a bounding box for the
    ellipse {E} with integer coordinates. The box is at least {mrg}
    units away from the ellipse, on each side. */

void ellipse_crs_to_ouv(ellipse_crs_t *E, ellipse_ouv_t *F);
  /* Stores in {*F} the dir-radius form of the ellipse {*E},
    except that {*F} will be centered at the origin. 
    
    If {*E} is a circle, {F.u,F.v} are set to two arbitrary
    orthonormal vectors, and {F.F.a,F.F.b} are set to {E.rad}. */
  
void ellipse_ouv_to_crs(ellipse_ouv_t *F, r2_t *ctr, ellipse_crs_t *E);
  /* Stores in {*E} the center-radius-stretch form of the ellipse {*F},
    with {*ctr} as the center. */
  
bool_t ellipse_crs_inside(ellipse_crs_t *E, r2_t *p);
  /* Returns TRUE if {p} is inside {E}, FALSE if outside.  
    May return either value if {p} is very close to the boundary. */

r2_t ellipse_crs_relative_coords(ellipse_crs_t *E, r2_t *p);
  /* Returns the coordinates of {p} relative to the ellipse {E}.
    Namely, the image of {p} by the transformation that takes the
    center to the origin and makes the major and minor semidiameters
    of {E} equal to 1, preserving their directions.  Thus,
    {p} is inside {E} if an only if the resut is inside the 
    unit disk. */

double ellipse_crs_position(ellipse_crs_t *E, r2_t *p, r2_t *csp);
  /* Returns the position of {p} relative to the boundary of the
    ellipse {E}, along the ray from the center of {E} to {p}. The
    result is less than 1 if {p} is inside, and greater than 1 if {p}
    is outside. It is 0 if {p} is at the center, and 1 if it is on the
    boundary.
    
    If {csp} is not NULL, the procedure also stores in {*csp} the
    cosine and sine of the angular argument of {p}, measured from the
    vector {E.str}. */

double ellipse_crs_nearest_point(ellipse_crs_t *E, r2_t *p, r2_t *q);
  /* Finds the signed distance from {p} to the boundary of the ellipse {E}.
    If {q} is not NULL, the procedure also stores in {*q} the point on
    the boundary of the ellipse {E} that is closest to {p}. 
    
    The return value {+dist(p,q)} if {p} is outside, and {-dist(p,q)}
    if {p} is inside.
    
    If you have many points to test against the same elipse, consider
    using {ellipse_crs_to_ouv} and then {ellipse_ouv_nearest_point},
    which is slightly faster. */

double ellipse_crs_border_position(ellipse_crs_t *E, double hwd, r2_t *p);
  /* If the distance from {p} to {E} is less than {hwd}, returns that
    distance divided by {hwd}. Otherwise returns {+1} if {p} is
    outside the ellipse, {-1} if it is inside.
    
    This is the relative position of {p} in the boundary of {E}, when
    it is painted with a round brush o radius {hwd}. This procedure is
    much faster than {ellipse_ouv_nearest_point} when {p} is not close
    to the stroked region. */

void ellipse_crs_print(FILE *wr, ellipse_crs_t *E, char *fmt);
  /* Writes to {wr} the geometric parameters {E}.  Each parameter
    is printed with format {fmt} (which should be appropriate 
    for a {double} value). */

#endif
 
