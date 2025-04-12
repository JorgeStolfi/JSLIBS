#ifndef pst_integrate_H
#define pst_integrate_H

/* pst_integrate.h -- simple procedures for integratings slope maps. */
/* Last edited on 2025-04-03 17:35:31 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_imgsys_solve.h>
#include <pst_slope_map.h>
#include <pst_height_map.h>
   
typedef void pst_integrate_report_data_proc_t
  ( int32_t level,
    float_image_t *G,
    float_image_t *H,
    float_image_t *R
  ); 
  /* Type of a client-given procedure that may be called by recursive
    integrators to report the reduced slope map {G}, the reduced
    hints (external height estimate) map {H}, and the reference height
    map {R} used at each scale. */   
          
typedef void pst_integrate_report_heights_proc_t
  ( int32_t level,
    int32_t iter,
    double change,
    bool_t final,
    float_image_t *Z,
    float_image_t *R
  ); 
  /* Type of a client-given procedure that may be called
    by slope integration procedures to report the current height map {Z} in each scale.
    The argument {iter} should be the number of iterations already done
    (0 = initial guess) and {change} should be the max height change from the 
    previous iteration (irrelevant for the initial guess).  The {final} arg should
    be true if the iteration has stopped.
    
    If the {Z} image has two or more channels, {Z[1,X,Y]} will be some
    finite non-negative reliability weight of the height {Z[0,X,Y]},
    such as the weight {eq.wtot} of the corresponding equation in the
    system. If {Z[1,X,Y]} is zero, the height {Z[0,X,Y]} should be
    considered undefined.
    
    The parameter {R}, if not {NULL}, will be a reference height map
    with same col and row counts as {Z}, that the procedure may use
    to evaluate the error of the current solution.  If {R} has two channels,
    channel 1 will be interpreted as a reliability weight, as that of {Z} above. */   

/* THE INTEGRATION SYSTEM */

pst_imgsys_t* pst_integrate_build_system
  ( float_image_t *G,
    float_image_t *H,
    double hintsWeight,
    bool_t extrapolate,
    int32_t indent,
    bool_t verbose
  );
  /*
    Allocates and returns a {pst_imgsys_t} system {S} whose equations define
    the height map {Z} given a slope map {G}.  See {pst_integrate_build_system_ARGS_INFO("hintsWeight")}
    and {pst_integrate_build_system_RESULT_INFO} for details. 
    
    The {extrapolate} flag is passed to {pst_slope_map_get_edge_data} (q. v.). It 
    authorizes the procedure to use linear extrapolation when 
    
    Messages to {stderr} are indented by {indent} spaces. */
    
#define pst_integrate_build_system_ARGS_INFO(hintsWeight) \
  "The input slope map {G} must have two or three channels. The value of {G[0,x,y]}  must" \
  " be the X derivative {dZ/dX} of the heigh map, averaged over the domain pixel with" \
  " indices {[x,y]}; that is, with opposite corners {(x,y)} and {(x+1,y+1)}.  The value of {G[1,x,y]} must" \
  " be the Y derivative {dZ/dY} averaged over the same pixel.\n" \
  "\n" \
  "  If the slope map has three channels, the value of {G[2,x,y]} must be a finite" \
  " non-negative number which is interpreted as the relative weight or reliability of the" \
  " slope data {G[c,x,y]}.  The smaller the weight, the smaller will be the influence" \
  " of the slopes {G[c,x,y]} on the computed heights.  If the slope map has only two" \
  " channels, all reliability weights are assumed to be 1.  If either {G[0,x,y]} or {G[1,x,y]} are" \
  " not finite, the weight is assumed to be zero.  If the weight is zero, the" \
  " values of {G[0,x,y]} or {G[1,x,y]} are ignored.\n" \
  "\n" \
  "  The desired height map {Z} will be a two-channel" \
  " image. Channel 0 will contain the {Z} values of the terrain computed from the" \
  " given slopes.   Channel 1 will have reliability weights for those height values.  If" \
  " the slope map {G} has {NX_G} columns and {NY_G} rows, the" \
  " desired height map {Z} will have {NX_Z=NX_G+1} columns" \
  " and {NY_Z=NY_G+1} rows. Specifically, for all {x} in {0..NX_Z-1} = {0..NX_G} and" \
  " all {y} in {0..NY_Z-1} = {0..NY_G},  sample {Z[0,x,y]} will be the height at" \
  " the integer point {(x,y)}; which is the lower CORNER of pixel {[x,y]} of {G}, and" \
  " not the pixel center {(x+0.5,y+0.5)}.\n" \
  "\n" \
  "  If the hints map {H} is provided, it must be a height map with one or two channels and" \
  " same size as {Z}. The value of {H[0,x,y]} is taken to be an independent" \
  " estimate of the height {Z[0,x,y]}.  If {H} has two channels, {H[1,x,y]} should" \
  " be a finite non-negative number which is taken to be the" \
  " reliability weight of that estimate; otherwise all the reliability weights are" \
  " assumed to be 1.  If {H[0,x,y]} is not finite, its weight is assumed to be zero. " \
  " If the weight of {H[0,x,y]} is nonzero, the equation for the height at grid" \
  " vertex {(x,y)} will include a term that pulls" \
  " the solution {Z[0,x,y]} towards {H[0,x,y]} with that weight times the overall weight" \
  " factor {" hintsWeight "}.\n" \
  "\n" \
  "  If there is no {H} hint for the height of a grid vertex {(x,y)}, and" \
  " all {G} weights on the edges around" \
  " that vertex are equal, the equations essentially state that the Laplacian of the height, as" \
  " computed from the (unknown) height map {Z}, shall be equal to the Laplacian computed" \
  " from the (known) slope map {G}.  Both Laplacians are computed by weighted finite-difference" \
  " formulas, properly adjusted along the edges and at the places where the weights are zero.\n" \
  "\n" \
  "   If there is no {H} hint term for some variable, and the {G} weights of the edge" \
  " terms around that vertex turn out to be all zero, the equation would be {0 = 0} which" \
  " could cause problems in the solution of the system.  In such cases, the procedure" \
  " sets the equation to {Z[0,x,y] = 0}.\n" \
  "\n" \
  "  Note that the slope map {G} specifies only differences between the height" \
  " values {Z[0,x,y]}. Therefore, if there are no {H} hints, the system of equations" \
  " described above is indeterminate, since adding a constant to all {Z[0,x,y]} values" \
  " would not affect the equations.\n" \
  "\n" \
  "  In fact, depending on the {G} weights, there may exist a proper subset {U} of unknown" \
  " variables {Z[0,x,y]} so that no equation relates variables of {U} to variables not" \
  " in {U}.  In that case, if there is no {H} hint for those variables one could add" \
  " an arbitrary constant to all of them without breaking the equations.  Thus, if the" \
  " variables can be partitioned into {m} such subsets, the whole system will have" \
  " family of solutions with dimension up to {m}."
  
#define pst_integrate_build_system_RESULT_INFO \
   "The linear system {S} will have {S.N} equations and {S.N} variables, where {S.N} will" \
  " be {NXY_Z=NX_Z*NY_Z}.   It is {A z = b} where {A} is a known square matrix of" \
  " size {S.N}, {z} is the column vector of the {S.N} unknown height values, and {b} is" \
  " a known column vector of size {S.N}.  Specifically, the height map sample {Z[0,x,y]} will" \
  " correspond to the variable {z[x + y*NX_Z]}.\n" \
  "\n" \
  "  In any case, the procedure stores in {S.col,S.row} two tables that tell the" \
  " correspondence between the height values and equations of the system and the elements" \
  " of the height map {Z}. Namely, for each {k} in {0..S.N-1}, height value {z[k]} is" \
  " the height map element {Z[0,S.col[k],S.row[k]]}."

#endif
