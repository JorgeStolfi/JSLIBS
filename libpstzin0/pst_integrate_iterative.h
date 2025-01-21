#ifndef pst_integrate_iterative_H
#define pst_integrate_iterative_H

/* procedures for integratings slope maps by smultiscale ieration. */
/* Last edited on 2025-01-15 14:30:54 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_imgsys_solve.h>
#include <pst_slope_map.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>

#include <pst_integrate.h>

void pst_integrate_iterative
  ( float_image_t *G, 
    float_image_t *W, 
    bool_t keepNull,
    float_image_t *Z, 
    float_image_t *U,
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    bool_t verbose,
    int32_t level,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  );
  /* Fills {Z} with the depth map that best matches the slope map {G}, where channel
    0 is {dZ/dX} and channel 1 is {dZ/dY}.  The image {Z} must have one channel and
    must be bigger than {G} by one col and one row
    
    If {W} is not NULL, it is taken to be a single-channel weight map, 
    with the same dimensions as the slope map {G}.
    
    The height map {Z} is computed by solving a set of linear equations,
    constructed from {G} and {W} using {pst_integrate_build_system}
    (q.v.).  
    
    If {keepNull} is true, the system is fudged so that any height map
    pixel that is completely surrounded by zero-weight data (so that its
    equation has only one term with weight zero) will be ultimately be
    set to the average of their neighbors. Thus the regions of weight
    zero in the input data will be filled with a Laplace equilibrium
    surface. This fudging does not affect the final result for pixels
    with non-null equations. If {keepNull} is false, those pixels will
    be excluded from the system, and will be set to zero in the final
    height map.
    
    The linear system is solved by the Gauss-Jordan iterative method.
    On entry, the contents of, which will be used as the initial guess
    for the iteration. Otherwise the initial guess will be all zeros.
    The iteration will stop when the maximum change in
    any height value is less than {convTol}, or after {maxIter}
    iterations, whichever happens first. If {topoSort} is TRUE,
    solves the equations in order of increasing equation weight {wtot}.
    
    IF {U} is not {NULL}, it must be a single-channel weight map with
    the same size as the height map {Z}. Each sample of {U} will be
    set to the weight {wtot} of the equation that defines the
    corresponding sample of {Z}.  Pixels of {Z} that end up with no equation 
    (if {keepHoles} is false) are set to zero.

    If {verbose} is TRUE, prints a report of the interations. 
    
    if {reportSys} is not {NULL}, the procedure calls {reportSys}
    once after building the system.
    
    If {reportHeights} is not null, the procedure calls
    {reportHeights(level,iter,change,final,OZ)} one or more times. Here
    {iter} is the number of complete Gauss-Seidel or Gauss-Jacobi
    iterations performed before the call; and {changeZ} is the max
    absolute change in any {OZ} sample since the previous iteration
    (meaningless when {iter=0}).  The {level} argument is used only here.
    
    The procedure {reportHeights} is always
    called once with {final=TRUE} after the last iteration, and once
    with {iter=0} and {final=FALSE} before the first iteration. If
    {reportStep} is not zero, {reportHeight} is also called with
    {final=FALSE} before each iteration whose index {iter} is less than 
    {reportStep} or a positive multiple thereof.  */

#endif
