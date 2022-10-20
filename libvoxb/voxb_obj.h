/* voxb_obj.h --- basic shape predicates for voxel-based modeling */
/* Last edited on 2022-10-20 05:48:22 by stolfi */

#ifndef voxb_obj_H
#define voxb_obj_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>

/* OBJECT PREDICATES

  These functions classify a point {p} relative to solid object,
  returning {TRUE} if the point is inside, {FALSE} if it is ouside. */

bool_t voxb_obj_ball(r3_t *p, double R);
  /* A ball with center at the origin and radius {R}.  */

bool_t voxb_obj_cube(r3_t *p, double R, double fillR);
  /* A cube with center at the origin and side {2*R}, 
    with sides parallel to the axes and a round fillet of radius {fillR} 
    around edges and corners.  */

bool_t voxb_obj_box(r3_t *p, double RX, double RY, double RZ, double fillR);
  /* A box with center at the origin and sides {2*RX} by {2*RY} by {2*RZ}, 
    with sides parallel to the axes and round fillets of radius {fillR} along
    edges and corners.  */

bool_t voxb_obj_rounded_box(r3_t *p, double RX, double RY, double RZ, double roundR, double fillR);
  /* A box with center at the origin, and sides {2*RX} by {2*RY} by
    {2*RZ} parallel to the axes, with the vertical edges rounded off
    with radius {roundR}, and fillets of radius {fillR} along the top
    and bottom edges. */

bool_t voxb_obj_donut(r3_t *p, double minR, double majR, int32_t axis);
  /* A donut (filled torus) with midline radius {majR}, and dough radius
    {minR}. The donut is centered at the origin, and is rotationally
    symmetric around the coordinate axis {ax} (0, 1, or 2). */

bool_t voxb_obj_rod(r3_t *p, double H, double R, double fillR);
  /* A cylindrical rod with total length {2*H}, radius {R}, and fillets
    along both edges with radius {fillR}. The rod is centered at the
    origin and has the {Z} axis as the cylinder axis. */

bool_t voxb_obj_tube(r3_t *p, double H, double Ri, double Ro, double fillR);
  /* A cylindrical tube with total length {2*H}, inner radius {Ri},
    outer radius {Ro}, and fillets along all four edges with radius
    {fillR}. The tube is centered at the origin and has the {Z} axis as
    the cylinder axis. */

bool_t voxb_obj_axes(r3_t *p, double L);
  /* A representation of the coordinate axes {X,Y,Z} as three sticks
    of length {L} with {0,1,2} balls at the tip, respectively.  */

bool_t voxb_obj_round_cup(r3_t *p, double H, double R, double thk, double fillR);
  /* A cylindrical cup with total height {2*H}, outer radius {R}, and
    walls of thickness {thk}, open at the top, with a flat closed bottom,
    centered at the origin. The cup has a fillet around the base with
    radius {fillR} and a rounded top edge.  Requires {fillR >= thk}. */

#endif

