/* voxm_obj.h --- basic clipped distance functions for voxel-based modeling */
/* Last edited on 2021-06-06 03:20:27 by jstolfi */

#ifndef voxm_obj_H
#define voxm_obj_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>

/* FUZZY OBJECT CLASSIFICATION FUNCTIONS

  These functions classify a point {p} relative to solid object,
  returning an indicator value that is greater than 0.5 inside the
  object, less than 0.5 outside. The indicator is fuzzyfied for proper
  antialiasing when extracting the object surface (e.g. by marching
  tetrahedra).
  
  Specifically, if {p} is at least {fuzzR} away from the boundary the
  object, the functions return 1 if {p} is inside the object and 0 if it
  is outside it. If {p} is less than {fuzzR} away from the boundary, the
  result varies linearly with that distance, from 0.5 to 1.0 inside,
  from 0.5 to 0.0 outside; so that the indicator is a continuous
  function of {p}. The 0.5 isosurface is the surface of the
  corresponding non-fuzzy object.
  
  If {fuzzR=0}, in particular, the function simply returns 
  0 outside, 1 inside the respective object. */

double voxm_obj_cube(r3_t *p, double R, double fillR, double fuzzR);
  /* A fuzzy cube with center at the origin and side {2*R}, 
    with sides parallel to the axes and a round fillet of radius {fillR} 
    around edges and corners.  */

double voxm_obj_box(r3_t *p, double RX, double RY, double RZ, double fillR, double fuzzR);
  /* A fuzzy box with center at the origin and sides {2*RX} by {2*RY} by {2*RZ}, 
    with sides parallel to the axes and a round fillets of radius {fillR} along
    edges and corners.  */

double voxm_obj_rounded_box(r3_t *p, double RX, double RY, double RZ, double roundR, double fillR, double fuzzR);
  /* A fuzzy box with center at the origin, 
    and sides {2*RX} by {2*RY} by {2*RZ} parallel to the axes, 
    with the vertical edges rounded off with radius {roundR},
    and a fillets of radius {fillR} along the top and bottom 
    edges.  */

double voxm_obj_ball(r3_t *p, double R, double fuzzR);
  /* A fuzzy ball with center at the origin and radius {R}.  */

double voxm_obj_donut(r3_t *p, double minR, double majR, int axis, double fuzzR);
  /* A fuzzy donut (filled torus) with midline radius {majR}, and dough
    radius {minR}. The donut is centered at the origin, and is
    rotationally symmetric around the coordinate axis {ax} 
    (0, 1, or 2). */

double voxm_obj_rod(r3_t *p, double H, double R, double fillR, double fuzzR);
  /* A cylindrical rod with total length {2*H}, radius {R}, and fillets along
    both edges with radius {fillR}.  The rod is centered at the origin
    and has the {Z} axis as the cylinder axis. */

double voxm_obj_tube(r3_t *p, double H, double Ri, double Ro, double fillR, double fuzzR);
  /* A cylindrical tube with total length {2*H}, inner radius {Ri}, outer radius {Ro}, 
    and fillets along all four edges with radius {fillR}.  The tube is centered at the origin
    and has the {Z} axis as the cylinder axis. */

double voxm_obj_axes(r3_t *p, double L, double fuzzR);
  /* A representation of the coordinate axes {X,Y,Z} as three sticks
    of length {L} with {0,1,2} balls at the tip, respectively.  */

double voxm_obj_round_cup
  ( r3_t *p, 
    double H,
    double R, 
    double thk,
    double fillR,
    double fuzzR
  );
  /* A cylindrical cup with total height {2*H}, outer radius {R}, and
    walls of thickness {thk}, open at the top, with a flat closed bottom,
    centered at the origin. The cup has a fillet around the base with
    radius {fillR} and a rounded top edge.  Requires {fillR >= thk}. */

/* CLIPPED RELATIVE DISTANCE */

double voxm_obj_classify(double d, double fuzzR);
  /* Given the signed distance {d} from some point to the boundary of an object
    (positive outside, negative inside), returns the corresponding fuzzified indicator
    value, with a fuzzy layer {fuzzR} thick on each side of the boundary.
    
    Specifically, the resut is 1 if {d >= fuzzR}, 0 if {d < -fuzzR}, and interpolates
    linearly between twose two points when {d} is in the range {(-fuzzR _ + fuzzR)}.
    
    Note that the distance {d} need not be accurate for points that are
    more than {fuzzR} away from the surface, as long as those outside
    have {d >= +fuzzR} and those inside have {d <= -fuzzR}. */

#endif

