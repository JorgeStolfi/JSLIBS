/* Plots dimensional specs as two lines, two arrows, and a number. */
/* Last edited on 2020-10-27 16:39:35 by jstolfi */

#ifndef epswr_dim_H
#define epswr_dim_H

#include <epswr.h>
  
void epswr_dim_linear
  ( epswr_figure_t *eps, 
    double xa, double ya,
    double xb, double yb,
    double *dabP,
    double agap,  double bgap, 
    double elen,
    bool_t inner,
    double dpos, double dlen, 
    double hoff, double voff,
    double *xrP, double *yrP,
    double *rotP
  );
  /* Draws a grahical indication of the distance between points
    {a=(xa,ya)} and {b=(xb,yb)}.  
    
    If {dabP} is not {NULL}, the procedure sets {*dabP} to the Euclidean
    distance between {a} and {b}.
    
    Specifically, let {A,B} be the lines that are perpendicular to the
    line {a--b} through those two points. The procedure draws an
    /extension segment/ on each of those two lines, then a /dimension
    segment/ parallel to {a--b} with tips on those lines.
    
    The extension segment on {A} will start {agap} away from point {a}. The
    extension segment on {B} will start {bgap} away from point {b}. Both
    {agap} and {bgap} are in CLIENT units.  In each case, the extension start 
    will be to the left of the {a--b} line if the parametes is positive,
    and to the right if negative.
    
    If {elen} is positive, the extension segments will extend for
    further {elen} millimeters to the left beyond the leftmost of the
    two extension start points. If {elen} is negative, they will extend
    {|elen|} millimeters to the right, beyond the rightmost extension
    start point. Note that {elen} is in MILLIMETERS, irrespective of the
    plotting scale.
    
    The dimension segment will be paralell to the {a--b} line and
    will be placed {dpos} millimeters before the end of the extension
    segments.  If negative, the dimension segment will be placed beyond the
    extension segments.
    
    If {inner} is true, the dimension segment will be drawn spanning the
    strip between {A} and {B}. In that case, if {dlen} is greater than
    half the millimeter distance between the lines, the dimension line will be a
    single two-headed arrow, otherwise it will be two arrows, each with
    total length {dlen} millimeters.  If {inner} is false, the dimension segment will be
    two converging arrows outside that strip, each with total length
    {dlen} millimeters. 
    
    Note that the arrowheads will have the same fixed
    size in any case. If {dlen} is zero or negative, the dimension segment
    and the arrowheads will be omitted altogether.
    
    If {xrP,yrP} are not {NULL}, the procedure returns in {*xrP,*yrP}
    the coordinates of a reference point that can be used to place a
    label for this dimension. Let {m} be the point on the dimension segment
    that is halfway between lines {A} and {B}. The reference point {r}
    will be displaced from {m} by {hoff} in the direction from {a} to
    {b}, and {voff} in the perpendicular direction. Also, if {rotP} is
    not {NULL}, the procedure returns in {*rotP} the CCW angle from the
    {X}-axis to the vector from {a} to {b}, in the range {[-180 _
    +180]}.
    
    The coordinates {xa,ya,xb,yb} and the result {*dabP} are in Client
    coordinates.  All other dimensions and displacements are in millimeters,
    irrespective of the current Client scale. */

#endif
