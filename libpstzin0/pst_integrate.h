#ifndef pst_integrate_H
#define pst_integrate_H

/* pst_integrate.h -- simple procedures for integratings slope maps. */
/* Last edited on 2025-01-25 09:25:02 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_imgsys_solve.h>
#include <pst_slope_map.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>
   
typedef void pst_integrate_report_data_proc_t
  ( int32_t level,
    float_image_t *G,
    float_image_t *H,
    float_image_t *Z
  ); 
  /* Type of a client-given procedure that may be called
    by recursive integrators to report the reduced slope map {G}, 
    the reduced external height map estimate {H},
    and the initial guess {Z} used at each scale. */   
          
typedef void pst_integrate_report_heights_proc_t
  ( int32_t level,
    int32_t iter,
    double change,
    bool_t final,
    float_image_t *Z
  ); 
  /* Type of a client-given procedure that may be called
    by slope integration procedures to report the current height map in each scale.
    The argument {iter} should be the number of iterations already done
    (0 = initial guess) and {change} should be the max height change from the 
    previous iteration (irrelevant for the initial guess).  The {final} arg should
    be true if the iteration has stopped.
    
    If the {Z} image has two ro more channels, channel 1 i set to thee
    total reliability {eq.wtot} of the equation {eq} that defines each
    height value. It will be the same in all calls to this report
    procedure for a given integration problem. */   

/* THE INTEGRATION SYSTEM */

pst_imgsys_t* pst_integrate_build_system(float_image_t *G, float_image_t *H, bool_t verbose);
  /*
    Allocates and returns a {pst_imgsys_t} system {S} whose equations define
    the height map {Z} given a slope map {G}.

    The slope map {G} must have three channels. The value of {G[0,x,y]}
    must be the X derivative {dZ/dX} of the heigh map, averaged over the
    pixel {[x,y]} The value of {G[1,x,y]} must be the Y derivative
    {dZ/dY} averaged over the same pixel.
    
    The value of {G[2,x,y]} must be non-negative, and is interpreted as
    the relative weight or reliability of the slope data {G[c,x,y]}. The
    smaller the weight, the smaller will be the influence of the slopes
    {G[c,x,y]} on the computed heights.
    
    The samples of a height map {Z} are heights at the CORNERS of the
    grid of {G} pixels. Therefore, if the slope map {G} has {NX_G}
    columns and {NY_G} rows, the unknown height map {Z} will have {NX_Z
    = NX_G+1} columns and {NY_Z = NY_G+1} rows. Thus {Z[0,x,y]} is the
    height at the integer point {(x,y)}, not at the pixel center
    {(x+0.5,y+0.5)}, for all {x} in {0..NX_Z-1} = {0..NX_G} and all {y}
    in {0..NY_Z-1} = {0..NY_G}.

    If {H} is not {NULL}, it must be a height map with two channels. 
    The value of {H[0,x,y]} is taken to be independent estimate of the height 
    {Z[0,x,y]}; and {H[1,x,y]} is taken to be the reliability weight of that
    estimate.  The equations will include a term that pulls the solution {Z[0,x,y]}
    towards {H[0,x,y]} with that weight.

    The linear system has {S.N} equations and {S.N} variables. It is {A
    z = b} where {A} is a known square matrix of size {S.N}, {z} is the
    column vector of the {S.N} unknown height values, and {b} is a known
    column vector of size {S.N}. The size {S.N} will be {NXY_Z =
    NX_Z*NY_Z}. The height map sample {Z[0,x,y]} will correspond to the
    variable {z[x + y*NX_Z]}.
    
    In either case the procedure stores in {S.col,S.row} two tables that
    tell the correspondence between the height values and equations of
    the system and the elements of the height map {Z}. Namely, for each
    {k} in {0..S.N-1}, height value {z[k]} is the height map element
    {Z[S.col[k],S.row[k]]}.
    
    If there is no {H} estimate, and all {G} weights are equal, the
    equations essentially state that the Laplacian of the height, as
    computed from the (unknown) height map {Z}, is equal to the
    Laplacian computed from the (known) slope maps {G}. Both Laplacians
    are computed by weighted finite-difference formulas, properly
    adjusted along the edges and at the places where the weights are zero.
    
    The slope map {G} specifies only differences between the height values {Z[0,x,y]}.
    Therefore, if {H} is null, or some of the weights {H[1,x,y]} are zero, the
    system would be indeterminate, since adding a constant to all
    {Z[0,x,y]} values would not affect the equations.  In fact, 
    depending on the {G} weights, the unknown variables {Z[0,x,y]} could be
    partitioned into {m >= 2} /connected components/ so that no equation 
    relates variables from two distinct components.  In that case,
    the system would have an {m}-dimensional family of solutions,
    that correspond to adding a different constant to all variables
    in the same component. 
    
    To make the system determinate, if {H} is null or {H[1,x,y]} is
    zero, the procedure implicitly assumes {H[0,x,y] = 0} with a 
    very small but positive weight.  This has the effect of pulling
    the average of the solution {Z[0,x,y]} in each connected 
    component towards zero. */

#endif
