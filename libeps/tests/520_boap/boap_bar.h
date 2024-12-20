/* Bar elements for the Boaretto 113 side gate drawings. */
/* Last edited on 2024-12-05 10:16:02 by stolfi */

#ifndef boap_bar_H
#define boap_bar_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <sign.h>
#include <vec.h>
#include <frgb.h>
#include <r3.h>
#include <r2.h>
#include <i3.h>

#include <boap_plate.h>

/* All dimensions are in mm. */

typedef struct boap_bar_t
  { char type;               /* Bar type. */
    /* World coordinates: */
    r3_t ctr;           /* World coordinates of nominal center of bar. */
    i3_t wax;           /* World axes that correspond to local axes. */
    /* Cross section parameters (in the local coordinate system): */
    r3_t size;          /* Bar extent along each local axis. */
    r2_t thk;           /* Thickness of cross-section elements. */
    r2_t emin_trim;     /* Length trimming of each part at {emin} end. */
    r2_t emax_trim;     /* Length trimming of each part at {emax} end. */
    r2_t emin_acut;     /* Angle of cut of bar at the {emin}  end. */
    r2_t emax_acut;     /* Angle of cut of bar at the {emax}  end. */
    double irad;        /* Roundoff radius of inner corners. */
    double orad;        /* Roundoff radius of outer corners. */
    frgb_t *color;      /* Fill color for plotting. */
  } boap_bar_t;
  /* A record that defines a steel bar.

    CROSS-SECTON TYPES
    
    The shape of the cross-section is defined by the {type}. Currently,
    valid types are 'F' for a flat iron, 'T' for a 'T'-iron, 'L' for a
    corner iron, 'U' for an 'U'-iron, 'M' for a metalon (rectangular-section tube). 
    
    The bar has its own coordinate system with axes {A,B,C}, where {C}
    (axis 2) is along the bar's length, and {A,B} (axes 0 and 1) are
    perpendcular to it. The origin of the local coordinate
    system is the point with world coordinates {(ctr[0],ctr[1],ctr[2])}.
    Axis {k} of the local coordinate system is parallel to
    world coordinate axis {wax[k]}, for {k} in {0..2}.
    
    In principle, the bar extends along the local {C} axis
    from coordinate {-size[2]/2} to {+size[2]}, except as 
    modified by trimming and cutting (see below).
    
    The cross-section is a two-dimensional figure in a generic {A,B}
    plane. In the descriptions below, the {A} axis (breadth, width) is
    assumed to be horizontal, and the {B} axis (height, tallness)
    vertical, pointing up.
    
    The fields {size[0]} and {size[1]} define the extent of the
    cross-section along the {A} and {B} axes, respectively. The fields
    {thk[0]} and {thk[1]} define the thickness of the cross-section
    elements along the SAME axes.

      If {type} is 'F', the cross-section is a rectangle that extends
      from {-size[0]/2} to {+size[0]/2} along the {A} axis and from
      {-thk[1]/2} to {+thk[1]/2} along the {B} axis. The fields
      {size[1],thk[0]} are ignored.

      If {type} is 'T', the cross-section is a symmetrical UPSIDE-DOWN 'T'
      shape. The horizontal part (the "arms" of the 'T') extends from
      {-size[0]/2} to {+size[0]/2} along the {A} axis and from 0 to
      {thk[1]} along the {B} axis. The vertical part (the "leg" of the
      'T') extends from 0 to {size[1]} allong the {B} axis and from {-thk[0]/2} to
      {+thk[0]/2} along the {A} axis.

      If {type} is 'L', the cross-section is an 'L' shape. The horizontal
      part (the "foot" of the 'L') extends from {0} to {size[0]} along the
      {A} axis and from 0 to {thk[1]} along the {B} axis. The
      vertical part (the "leg" of the 'L') extends from 0 to {size[1]}
      along the {B} axis and from 0 to {thk[0]} along the {A} axis.

      If {type} is 'U', the cross-section is an 'U' shape with square
      bottom. The horizontal part (the "base" of the 'U') extends from
      {-r} to {+r} along the {A} axis, where {r} is {size[0]/2} and from
      0 to {thk[1]} along the {B} axis. The vertical parts (the "horns"
      of the 'U') extend from 0 to {size[1]} along the {B} axis and from
      {-e} to {-e+thk[0]}, and from {+e-thk[0]} to {+e}, along the {A}
      axis.

      If {type} is 'M', the cross-section is a hollow rectangle that
      extends from {-size[0]/2} to {+size[0]/2} along the {A} axis and
      from {-size[1]/2} to {+size[1]/2} along the {B} axis. The vertical
      walls (parallel to the {B} axis) are {thk[0]} thick, and the
      horizontal walls (parallel to the {A} axis) are {thk[1]} thick.

    Any of the parameters {size[0..1]} and {thk[0..1]} may be negative
    to specify that the corresponding part of the cross-section is directed
    towards negative coordinates. However, {size[0]} is expected to have
    the same sign are {thk[0]}, and ditto for {size[1]} and {thk[1]}.
    
    PARTS OF THE BAR
    
    The "{A} parts" of the bar are the flat iron bars that corresponds
    to the horizontal parts of the cross-section: the "arms" of a 'T',
    the "foot" of the 'L', the "base" of the 'U', or the two the
    horizontal walls of an 'M'. For an 'F' bar, it is the bar itself.
    
    Similarly, the "{B} parts" of the bar are the flat iron bars that
    corresponds to the vertical parts of the cross-section: the "leg" of
    the 'T' or 'L', or the two the vertical walls of an 'U' or 'M'. An 'F' bar
    has are no {B} parts.
    
    NOMINAL BAR AXIS AND CENTER
    
    The "nominal axis" of the bar is the {C} axis of the local coordinate system.
    In world coordinates, it is the line parallel to axis {wax[2]}
    that goes through the nominal center {ctr[0..2]}.  Note that:
    
      For 'F' and 'M' bars, the axis line goes through the center of the 
      cross section.
      
      For 'L' bars, the axis line runs along the outer corner edge, that is,
      at the "heel" of the cross-section.
      
      For 'T' bars, the axis line runs along the midline of the flat
      surface that corresponds to the edge of "arms" of the
      cross-section, on the side opposite to the "leg".
      
      For 'U' bars, the axis line runs along the midline of the flat
      surface that corresponds to the edge of "base" of the
      cross-section, on the side opposite to the "horns".
    
    TRIMMING 
    
    If {emin_trim[0]} is non-zero, the {A} parts of the bar will be
    trimmed off by that amount at the beginning end (the end with local
    {C} coordinate {-size[2]/2}. Similarly if {emax_trim[0]} is
    non-zero, the {A} parts will be trimmed at the opposite end 
    (with local {C} coordinate {+size[2]/2}) by that
    amount.  Thus the {A} parts of the bar will actually extend from
    {-size[2]/2+emin_trim[0]} to {+size[2]/2-emax_trim[0]}.
    
    The parameters {emin_trim[1]} and {emax_trim[1]} have the 
    similar meaning for the {B} parts of the bar.
    
    ANGULAR CUTTING
    
    If {emin_acut[0]} is not 0, the bar will be cut by a plane
    that makes that angle with the {A} axis, is parallel to 
    the {B} axis, and goes through the nominal start of the bar,
    that is, the point with local coordinate {-size[2]/2} on the {C} axis.
    Similarly, {emax_acut[0]} specifies an oblique cut at the other end,
    the point {+size[2]/2} on the {C} axis. 
    
    The parameters {emin_acut[1]} and {emax_acut[1]} specify similar
    oblique cuts, but with planes parallel to the {A} axis that 
    make those specified angles with the {B} axis. 
    
    All cut angles are in degrees. */

boap_bar_t *boap_bar_new
  ( char type,
    double Xctr,
    double Yctr,
    double Zctr,
    int8_t axA, double Asize, double Athk,
    int8_t axB, double Bsize, double Bthk,
    int8_t axC, double Csize,
    frgb_t *color
  );
  /* Allocates and returns a bar element of given {type},
    with the given nominal center, axes, extents, and thicknesses along each axis.
    The ends will be untrimmed and cut perpendicularly. 
    
    Note that the {*color} record is shared, not copied. */
    
void boap_bar_trim
  ( boap_bar_t *bar,
    int8_t bax,
    int8_t end,
    double amount
  );
  /* Trims the specified parts of the bar by the specified {amount}.
    
    The {bax} argument specifies the parts to be trimmed, the {A} parts ({bax=0})
    or the {B} parts ({bax=1}).
    
    The {end} argument specifies which end of the bar should be trimmed.
    If {end} is {-1}, trims the end of that part with local coordinate {C =
    -size[2]/2} and if {end} is {+1}, trims the end with local {C = +size[2]/2}. */

typedef boap_bar_t *boap_bar_ref_t;

vec_typedef(boap_bar_ref_vec_t, boap_bar_ref_vec, boap_bar_ref_t);

void boap_bar_smash 
  ( boap_bar_t  *bar,
    boap_plate_ref_vec_t *plv,
    int32_t *npP
  );
  /* Decomposes the plate {bar} into one or more plates and
    computes their world coordinate ranges.
    
    Assumes that the vector {plv} has {np} plates in {plv.e[0..np-1]},
    where {np} is the input value of {*npP}. Appends the new
    plates, reallocating {plv.e} as needed, and updating {*npP}. */

void boap_bar_print(FILE *wr, boap_bar_t *bar);
  /* Prints to {wr} a legible description of the given {bar}. */

#endif 
