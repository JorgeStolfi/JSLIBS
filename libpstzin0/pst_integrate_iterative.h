#ifndef pst_integrate_iterative_H
#define pst_integrate_iterative_H

/* procedures for integratings slope maps by smultiscale ieration. */
/* Last edited on 2025-01-08 00:17:58 by stolfi */

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
  ( float_image_t *IG, 
    float_image_t *IW, 
    bool_t keepNull,
    float_image_t *IZ, 
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    float_image_t **OZP,
    float_image_t **OWP,
    bool_t verbose,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  );
  /* Computes a depth map {OZ} given the slope map {IG}, where channel
    0 is {dZ/dX} and channel 1 is {dZ/dY}. The {OZ} image is returned
    in {*OZP}; it will be allocated by the procedure and will be
    bigger than {IG} by one row and one column.
    
    If {IW} is not NULL, it is taken to be a single-channel weight map, 
    with the same dimensions as the slope map.
    
    The image {OZ} is computed by solving a set of linear equations,
    constructed from {IG} and {IW} using {pst_integrate_build_system}
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
    If {IZ} is not {NULL}, it must be a height map, bigger than {IG}
    by one col and one row, which will be used as the initial guess
    for the iteration. Otherwise the initial guess will be all zeros.
    The iteration will stop when the maximum change in
    any height value is less than {convTol}, or after {maxIter}
    iterations, whichever happens first. If {topoSort} is TRUE,
    solves the equations in order of increasing equation weight {wtot}.
    
    IF {IW} and {OWP} are not {NULL}, the procedure also creates a weight map {OW}
    for the result {OZ}, returned in {*OWP}. It will be a
    single-channel map with the same size as {OZ}. Each sample of {OW}
    will be set to the weight {wtot} of the equation that defines the
    corresponding sample of {OZ}.

    If {verbose} is TRUE, prints a report of the interations. 
    
    If {reportHeights} is not null, the procedure calls
    {reportHeights(0,iter,change,final,OZ)} one or more times. Here
    {iter} is the number of complete Gauss-Seidel or Gauss-Jacobi
    iterations performed before the call; and {changeZ} is the max
    absolute change in any {OZ} sample since the previous iteration
    (meaningless when {iter=0}). The procedure {reportHeights} is always
    called once with {final=TRUE} after the last iteration, and once
    with {iter=0} and {final=FALSE} before the first iteration. If
    {reportStep} is not zero, {reportHeight} is also called with
    {final=FALSE} before each iteration whose index {iter} is a positive
    multiple of {reportStep}. */

#endif
