#ifndef pst_integrate_iterative_H
#define pst_integrate_iterative_H

/* procedures for integratings slope maps by smultiscale ieration. */
/* Last edited on 2025-03-05 15:06:32 by stolfi */

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
    float_image_t *H, 
    double hintsWeight,
    float_image_t *Z, 
    float_image_t *R, 
    uint32_t maxIter,
    double convTol,
    bool_t sortSys,
    bool_t verbose,
    int32_t level,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  );
  /* Fills the height map {Z} with the height values that best match the
    slope map {G} and the optional height hints map {H}.
    
    The image {G} must have 3 channels, where channel 0 is {dZ/dX},
    channel 1 is {dZ/dY}, and channel 2 is the reliability weight for
    that data.
    
    The image {Z} must have two channels, and must be bigger than {G} by
    one col and one row. Channel 0 will be the integrated height values,
    and channel 1 will be a reliability weight for them.

    If {H} is not {NULL}, it must be a height map with one or two
    channels and same size as {Z}, which is taken to be an externally
    determined hint for {Z}, with weights multiplied by {hintsWeight}. The
    height map {Z} is computed by solving a set of linear equations,
    constructed by {pst_integrate_build_system(G,H,hintsWeight)} (q.v.).
    
    See {pst_integrate_iterative_INFO("maxIter","convTol","sortSys")} for
    details of the solution method. The initial guess for the iteration
    will be the value of {Z[0,x,y]} on entry. On exit, each sample
    {Z[1,x,y]} will be set to the weight {wtot} of the main equation
    that defines the value {Z[0,x,y]}.
    
    The image {R} is a reference height map, mostly intended for
    debugging and performance mesurements, and may be {NULL}. If not
    {NULL}, it must be a height map with the same col and row count as
    {Z}, which is passed unchanged to {reportHeights} below.

    If {verbose} is TRUE, prints a report of the interations. 
    
    if {reportSys} is not {NULL}, the procedure calls {reportSys}
    once after building the system.
    
    If {reportHeights} is not null, the procedure calls
    {reportHeights(level,iter,change,final,Z',R)} one or more times.
    Here {Z'} is the current computed height map, {iter} is the number
    of complete Gauss-Seidel or Gauss-Jacobi iterations performed before
    the call; and {change} is the max absolute change in any {Z'} sample
    since the previous iteration (meaningless when {iter=0}). The
    {level} argument is used only here.
    
    The procedure {reportHeights} is always called once with
    {final=TRUE} after the last iteration, and once with {iter=0} and
    {final=FALSE} before the first iteration. If {reportStep} is not
    zero, {reportHeight} is also called with {final=FALSE} before each
    iteration whose index {iter} is less than {reportStep} or a positive
    multiple thereof. */

#define pst_integrate_iterative_INFO(maxIter,convTol,sortSys) \
  "The linear system is solved by the Gauss-Seidel iterative" \
  " method. The iteration will stop when the maximum change" \
  " in any height value is less than {" convTol "}, or" \
  " after {" maxIter "} iterations, whichever happens first.\n" \
  "\n" \
  "  If {" sortSys "} is true, the equations are rearranged" \
  " so that the Gauss-Seidel method considers them in" \
  " decreasing order of their total weight {eq.wtot}."

#endif
