/* Axis-aligned rectangular boxes in R^d. */
/* Last edited on 2024-12-31 00:53:08 by stolfi */

#ifndef ia_box_H
#define ia_box_H

#include <stdint.h>

#include <ia.h>

#define ia_box_MAX_DIM 16

typedef Interval ia_box_Interval;
typedef uint8_t ia_box_Dim;
typedef uint8_t ia_box_Axis;
typedef uint16_t ia_box_Corner;
typedef enum{ LO = -1, MD = 0, HI = +1 } ia_box_Dir;

void ia_box_corner(ia_box_Dim d, ia_box_Interval *b, ia_box_Corner c, Float *x);
  /* Stores in {x} the corner of box {b} with corner number {c}. */

void ia_box_center(ia_box_Dim d, ia_box_Interval *b, Float *x);
  /* Stores in {x} the center of box {b}, rounding to nearest.
    The result always lies in the (closed) box. */
  
ia_box_Axis ia_box_max_axis(ia_box_Dim d, ia_box_Interval *b);
  /* The axis along which box {b} has maximum extent. */

Float ia_box_max_rad(ia_box_Dim d, ia_box_Interval *b);
  /* Half (rounded up) of the maximum extent of box {b} along any axis. */

double ia_box_radius(ia_box_Dim d, ia_box_Interval *b);
  /* The Euclidean radius of {b}, rounded up. May return {+\oo}
    even for finite boxes. */

void ia_box_from_corners(ia_box_Dim d, Float *x, Float *y, ia_box_Interval *b);
  /* Stores in {b} the box that encloses the two given points {x} and {y}
    of {R^d}, namely {b_i = [u_i _ v_i]} where {u_i = \min{x_i,y_i}},
    {v_i = \max{x_i,y_i}}. */

void ia_box_split(ia_box_Dim d, ia_box_Axis k, ia_box_Dir dir, ia_box_Interval *b, ia_box_Interval *h);
  /* Splits the box {b} along axis {k}, and stores in {h} the 
    (closed) half specified by the {dir} argument: {LO} for the
    low half, {HI} for the high half.  

    Either half may be degenerate (zero-width) if {b} is too thin along
    axis {k}. If {dir = MD}, stores in {h} the degenerate box which is
    the intersection of the two closed halves. */

#endif
