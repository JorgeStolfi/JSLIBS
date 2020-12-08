/* Plots dimensional specs as two lines, two arrows, and a number. */
/* Last edited on 2020-10-27 17:02:19 by jstolfi */

#ifndef epswr_dev_dim_H
#define epswr_dev_dim_H

#include <epswr.h>
  
void epswr_dev_dim_linear
  ( epswr_figure_t *eps, 
    double psxa, double psya,
    double psxb, double psyb,
    double psagap,
    double psbgap, 
    double pselen,
    bool_t inner,
    double psdpos, double psdlen, 
    double pshoff, double psvoff,
    double *psxrP, double *psyrP,
    double *rotP
  );
  /* Draws a grahical indication of the distance between points
    {a=(psxa,psya)} and {b=(psxb,psyb)} (in absolute Device coordinates).
    
    Specifically, let {A,B} be the lines that are perpendicular to the
    line {a--b} through those two points. The procedure draws an
    /extersion segment/ on each of those two lines, then a /dimension segment/
    parallel to {a--b} with tips on those lines.
    
    The extension segment on {A} will start {psagap} away from {a}. The
    extension segment on {B} will start {psbgap} away from point {b}.  In each case, the extension start 
    will be to the left of the {a--b} line if the parametes is positive,
    and to the right if negative.
    
    If {pselen} is positive, the extension segments will extend for
    further {pselen} distance to the left beyond the leftmost of the
    two extension start points. If {pselen} is negative, they will extend
    {|elen|} distance to the right, beyond the rightmost extension
    start point. 
    
    The dimension segment will be paralell to the {a--b} line and
    will be placed at distance {psdpos} before the end of the extension
    segments.  If negative, the dimension segment will be placed beyond the
    extension segments.
    
    If {inner} is true, the dimension segment will be drawn spanning the
    strip between {A} and {B}. In that case, if {psdlen} is greater than
    half the distance between the lines, the dimension segment will be a
    single two-headed arrow, otherwise it will be two arrows, each with
    total length {psdlen}. If {inner} is false, the dimension segment will be
    two converging arrows outside that strip, each with total length
    {psdlen}.
    
    Note that the arrowheads will have the same fixed size in any case.
    If {psdlen} is zero or negative, the dimension segment and the
    arrowheads will be omitted altogether.
    
    If {psxrP,psyrP} are not {NULL}, the procedure returns in {*psxrP,*psyrP}
    the coordinates of a reference point that can be used to place a
    label for this dimension. Let {m} be the point on the dimension segment
    that is halfway between lines {A} and {B}. The reference point {r}
    will be displaced from {m} by {pshoff} in the direction from {a} to
    {b}, and {psvoff} in the perpendicular direction. Also, if {rotP} is
    not {NULL}, the procedure returns in {*rotP} the CCW angle from the
    {X}-axis to the vector from {a} to {b}, in the range {[-180 _
    +180]}.
    
    The signs of the displacements {psagap,psalen,psbgap,psblen,psdpos,psvoff} 
    are interpreted relative to the direction from {a} to {b}: to the left
    if positive, to the right if negative.
    
    All coordinates, lengths, and displacements are in points, in the absolute
    Device coordinate system. */

#endif
