#ifndef pst_integrate_recursive_H
#define pst_integrate_recursive_H

/* procedures for integratings slope maps by smultiscale ieration. */
/* Last edited on 2025-01-08 00:34:47 by stolfi */

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
  ( float_image_t *IG, 
    float_image_t *IW, 
    uint32_t level,
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    float_image_t **OZP,
    float_image_t **OWP,
    bool_t verbose,
    uint32_t reportStep,
    pst_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    pst_integrate_report_heights_proc_t *reportHeights
  );
  /* Computes a depth map {OZ} given the slope map {IG}, where channel
    0 is {dZ/dX} and channel 1 is {dZ/dY}. The {OZ} image is returned
    in {*OZP}; it will be allocated by the procedure and will be
    bigger than {IG} by one row and one column.
    
    If {IW} is not NULL, it is taken to be a single-channel weight map, 
    with the same dimensions as the slope map.
    
    The image {OZ} is computed by solving a set of linear equations.
    
    The linear system is solved by the Gauss-Jordan iterative method.
    The initial guess for this method is obtained by scaling down the
    slope maps {IG} by 1/2 in all three axes, computing the
    corresponding height field recursively, and un-scaling the result.
    The {level} parameter indicates the current depth of the recursion.
    The recursion stops when the maps are reduced to a single pixel.
    
    At each level, the iteration will stop when the maximum change in
    any height value is less than {convTol}, or after {maxIter}
    iterations, whichever happens first. If {topoSort} is TRUE,
    solves the equations in order of increasing equation weight {wtot}.
    
    If {IW} and {OWP} are not null, the procedure also creates a weight map {OW}
    for the result {OZ}, returned in {*OWP}. It will be a
    single-channel map with the same size as {OZ}. Each sample of {OW}
    will be set to the weight {wtot} of the equation that defines the
    corresponding sample of {OZ}.

    If {verbose} is TRUE, prints a report of the interations
    at each scale. 
    
    The procedures {reportData,reportSys,reportHeights}, if not NULL,
    are called at each level of the iteration. If {reportData} is not
    null, the procedure calls {reportData(level,IG,IW)} at the
    beginning of each level of the recursion, on the way up
    (increasing {level}). Note that {IW} may be null.
    
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
    called once with {final=TRUE} after the last iteration. If
    {reportStep>0}, it is also called with {final=FALSE} before each
    iteration whose index {iter} is multiple of {reportStep}; including
    before the first iteration, if any. In particular, if the {OZ} map
    is computed directly without iteration, {reportHeights} is called
    just once with {iter=0} and {final=TRUE}. Therefore it is never
    called twice with the same {level} and {iter}. */
    

#endif
