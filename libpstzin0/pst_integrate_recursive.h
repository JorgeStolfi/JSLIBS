#ifndef pst_integrate_recursive_H
#define pst_integrate_recursive_H

/* procedures for integratings slope maps by smultiscale ieration. */
/* Last edited on 2025-01-15 14:34:26 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_imgsys_solve.h>
#include <pst_slope_map.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>

#include <pst_integrate.h>

void pst_integrate_recursive
  ( float_image_t *G, 
    float_image_t *W, 
    bool_t keepNull,
    float_image_t *Z, 
    float_image_t *U,
    int32_t level,
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    bool_t verbose,
    pst_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  );
  /* Computes a depth map {Z} given the slope map {G}, where channel
    0 is {dZ/dX} and channel 1 is {dZ/dY}. The {Z} image must be
    bigger than {G} by one row and one column.
    
    If {W} is not NULL, it is taken to be a single-channel weight map, 
    with the same dimensions as the slope map {G}.
    
    The image {Z} is computed by solving a hierarchy of systems of linear 
    equations.  Each linear system is solved by the Gauss-Jordan iterative method.
    The initial guess for this method is obtained by scaling down the
    slope maps {G} by 1/2 in all three axes, computing the
    corresponding height field recursively, and un-scaling the result.
    The {level} parameter indicates the current depth of the recursion.
    The recursion stops when the maps are reduced to a single pixel.
    
    At each level, the iteration will stop when the maximum change in
    any height value is less than {convTol}, or after {maxIter}
    iterations, whichever happens first. If {topoSort} is TRUE,
    solves the equations in order of increasing equation weight {wtot}.  
    
    If {keepNull} is true, each system is fudged so that any height map
    pixel that is completely surrounded by zero-weight data (so that its
    equation has only one term with weight zero) will be ultimately be
    set to the average of their neighbors. Thus the regions of weight
    zero in the input data will be filled with a Laplace equilibrium
    surface. This fudging does not affect the final result for pixels
    with non-null equations. If {keepNull} is false, those pixels will
    be excluded from the system, and will be set to zero in the final
    height map.
    
    If {U} is not null, it must be a single-channel image 
    with the same weight map with
    the same size as the height map {Z}. Each sample of {U} will be
    set to the weight {wtot} of the equation that defines the
    corresponding sample of {Z}.  Pixels of {Z} that end up with no equation 
    (if {keepHoles} is false) are set to zero.

    If {verbose} is TRUE, prints a report of the interations
    at each scale. 
    
    The procedures {reportData,reportSys,reportHeights}, if not NULL,
    are called at each level of the iteration. If {reportData} is not
    null, the procedure calls {reportData(level,G,W)} at the
    beginning of each level of the recursion, on the way up
    (increasing {level}). Note that {W} may be null.
    
    If {reportSys} is not null, the procedure calls
    {reportSys(level,S)}, once at each level on the way down
    (decreasing {level}), to report the system {S}. Note that {S} may
    be NULL.  
    
    If {reportHeights} is not null, the procedure calls
    {reportHeights(level,iter,change,final,OZ)} one or more times for
    each level, on the way out of the recursion (decreasing {level}).
    Here {iter} is the number of complete Gauss-Seidel or Gauss-Jacobi
    iterations performed before the call; and {changeZ} is the max
    absolute change in any {OZ} sample since the previous iteration
    (meaningless when {iter=0}). The procedure {reportHeights} is always
    called once with {final=TRUE} after the last iteration, and once
    with {iter=0} and {final=FALSE} before the first iteration, if any. If
    {reportStep>0}, it is also called with {final=FALSE} before each
    iteration whose index {iter} is less than {reportStep} or a
    positive multiple thereof. In particular, if the 
    {Z} map is computed directly without iteration,
    {reportHeights} is called just once with {iter=0} and {final=TRUE}.
    Therefore it is never called twice with the same {level} and
    {iter}. */
    

#endif
