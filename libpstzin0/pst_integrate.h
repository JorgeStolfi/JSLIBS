#ifndef pst_integrate_H
#define pst_integrate_H

/* pst_integrate.h -- simple procedures for integratings slope maps. */
/* Last edited on 2025-01-08 00:29:25 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_imgsys_solve.h>
#include <pst_slope_map.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>
   
typedef void pst_integrate_report_data_proc_t(uint32_t level, float_image_t *IG, float_image_t *IW); 
  /* Type of a client-given procedure that may be called
    by recursive integrators to report the reduced slope map and weight map used at each scale. */   
          
typedef void pst_integrate_report_heights_proc_t(uint32_t level, uint32_t iter, double change, bool_t final, float_image_t *OZ); 
  /* Type of a client-given procedure that may be called
    by slope integration procedures to report the current height map in each scale.
    The argument {iter} should be the number of iterations already done
    (0 = initial guess) and {change} should be the max height change from the 
    previous iteration (irrelevant for the initial guess).  The {final} arg should
    be true if the iteration has stopped. */   

/* THE INTEGRATION SYSTEM */

pst_imgsys_t* pst_integrate_build_system
  ( float_image_t *G,
    float_image_t *W,
    bool_t verbose
  );
  /*
    Allocates and returns an {pst_imgsys_t} system {S} whose equations define
    the height map {Z} given a slope map {G}.

    The value of {G[0,x,y]} must be the X derivative {dZ/dX} of the heigh map,
    averaged over the pixel {[x,y]} The value of {G[1,x,y]} must
    be the Y derivative {dZ/dY} averaged over the same pixel.
    
    The samples of the {Z} map are heights at the CORNERS of
    the {G} pixels. Therefore, if the slope map {G} has {NX} columns
    and {NY} rows, the {Z} map has {NX+1} columns and {NY+1} rows.
    Thus {Z[x,y]} is the height at the integer point
    {(x,y)}, not at the pixel center {(x+0.5,y+0.5)}, for all {x} in
    {0..NX} and all {y} in {0..NY}.

    The image {W}, if not {NULL}, must have the same number of rows
    and columns as the slope map {G}. The value of {W[0,x,y]} is
    interpreted as the relative weight or reliability of the slope
    data {G[c,x,y]}. The absolute value of the weights is not
    relevant, only their ratios. The smaller {W[0,x,y]} is, compared
    to the weights of neighboring pixels, the smaller will be the
    influence of the slopes {G[c,x,y]} on the computed heights.
    The weights must be non-negative.  If {x,y} are outside the
    image domain then {W[0,x,y]} is assumed to be 0.

    The linear system has the form {A h = b} where {A} is a known matrix, {h}
    is the vector of unknown height values, and {b} is a known vector. The
    system has {S.N} equations and {S.N} unknowns.
    The size {S.N} will be {NZ = (NX+1)*(NY+1)}.
    
    Depending on the weight map, the equation that defines the height
    {Z[x,y]} may be {0*Z[x,y]=0}, leaving {Z[x,y]} indeterminate. The
    client may exclude those heights and equations from the system (see
    {pst_imgsys_remove_holes}), or fudge them so that they
    provide some definite value for the heights (see
    {pst_imgsys_fill_holes}).
    
    In either case the procedure stores in {S.col,S.row,S.uid} three
    tables that tell the correspondence between the height values and
    equations of the system and the elements of the {Z} array. Namely,
    for each {k} in {0..S.N-1}, height value {h[k]} is the height map element
    {Z[S.col[k],S.row[k]]}. Conversely, for each {x} in {0..NX} and each
    {y} in {0..NY}, the height map element {Z[x,y]} corresponds to the height value
    {h[k]} where {k = S.uid[x + y*(NX+1)]}, if {k >= 0}; os is omitted
    from the system, if {k<0}.
    
    If all weights are equal, the equations essentially state that the
    Laplacian of the height, as computed from the (unknown) height map
    {Z}, is equal to the Laplacian computed from the (known) slope
    maps {G}. Both Laplacians are computed by weights finite-difference
    formulas. Other equations, based on equality of slopes, are used
    along the edges and where the weights are zero.
    
    If all weights are non-zero, the system has a one-dimensional family
    of solutions, since the homogeneous system {A h = 0} is satisfied by
    {h[k] == 1} for all {k}. These solutions differ by adding the same
    constant to all height values. Thus the system only determines the
    differences between height values, not the actual heights.
    
    If enough weight elements {W[x,y]} are zero, the equations of the
    system may be partitioned into two or more independent connected
    components, each referring to a separate sets of height values, with
    no equation relating heights from different sets. In that case, each
    of these sub-systems has a one-dimensional family of solutions,
    independent from the others. That is, the system determines the
    differences in heights within each component, but not between
    different components. */

#endif
