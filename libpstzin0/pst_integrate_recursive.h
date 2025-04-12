#ifndef pst_integrate_recursive_H
#define pst_integrate_recursive_H

/* procedures for integratings slope maps by smultiscale ieration. */
/* Last edited on 2025-04-03 17:30:02 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_imgsys_solve.h>
#include <pst_slope_map.h>
#include <pst_height_map.h>

#include <pst_integrate.h>

void pst_integrate_recursive
  ( int32_t level,
    float_image_t *G, 
    float_image_t *H,
    double hintsWeight,
    bool_t extrapolate, 
    float_image_t *Z, 
    float_image_t *R,
    uint32_t maxLevel,
    uint32_t maxIter,
    double convTol,
    bool_t sortSys,
    bool_t verbose,
    pst_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  );
  /* Fills the height map {Z} with the height values that best match the
    slope map {G} and the optional height hints map {H}.  The 
    fitting is done by recursively solving a hierarchy of problems of 
    smaller and smaller sizes; see {pst_integrate_recursive_INFO} for details.
  
    The parameters {G,H,hintsWeight,extrapolate} as well as
    {Z,R,maxIter,convTol,sortSys} have the same requirements and
    meanings as in {pst_integrate_iterative} (q.v.), but apply to each
    level of the recursion.

    The {level} parameter indicates the current depth of the recursion.
    
    The {maxLevel} parameter specifies the maximum value of {level}. The
    recursion will stop when {level=maxLevel} or when the slope map is
    reduced to a single pixel, whichever happens first. If zero, the
    procedure is essentially equivalent to {pst_integrate_iterative}}.

    If {verbose} is TRUE, prints a report of the interations at each
    scale.
    
    The procedures {reportData,reportSys,reportHeights}, if not NULL,
    are called at each level of the iteration. If {reportData} is not
    null, the procedure calls {reportData(level,G',H',Z',R')} at the
    beginning of each level of the recursion, on the way in (in order of
    increasing {level}); where {G',H',Z',R'} are the input {G,H,R}
    reduced to the current level.
    
    If {reportSys} is not null, the procedure calls
    {reportSys(level,S')}, once at each level, on the way out of the
    recursion (in order of decreasing {level}), to report the system
    {S'} to be solved at that level. Note that {S'} may be {NULL} at the
    depest {level}.
    
    If {reportHeights} is not null, the procedure calls
    {reportHeights(level,iter,change,final,Z',R')} one or more times for
    each level, on the way out of the recursion (in order of decreasing
    {level}). Here {Z'} is the current computed height map, {R'} is the
    given reference maps {R} reduced to the scale of the current level,
    {iter} is the number of complete Gauss-Seidel or Gauss-Jacobi
    iterations performed before the call; and {change} is the max
    absolute change in any {Z} sample since the previous iteration
    (meaningless when {iter=0}).
    
    The procedure {reportHeights} is always called once per level with
    {final=TRUE} after the last iteration, and once with {iter=0} and
    {final=FALSE} before the first iteration, if any. If {reportStep0}
    is not zero, it is also called with {final=FALSE} before each
    iteration whose index {iter} is less than {reportStep} or a positive
    multiple thereof.  In particular, if the {Z} map gets computed
    directly without iteration, {reportHeights} is called just once with
    {iter=0} and {final=TRUE}. Therefore it is never called twice with
    the same {level} and {iter}. */

#define pst_integrate_recursive_INFO(hintsWeight,initialOpt,maxIter,convTol,sortSys) \
  "The height map is computed by solving systems of linear equations, created" \
  " from a hierarchy of slope maps which are reduced copies of the" \
  " given map {G} with" \
  " decreasing sizes.  Level zero of the hierarchy is the given" \
  " slope map {G}, and the slope map" \
  " at each subsequent level {r+1}" \
  " is a copy of the map of level {r} by a contant scale factor.  The" \
  " hints map {H} and the reference map {R}, if given, are also reduced" \
  " for each level in a consistent way.\n" \
  "\n" \
  "  At each level, starting from the highest" \
  " one (with the smallest map), a heightmap is computed from the corresponding" \
  " slope map, by solving a linear equation system.\n" \
  "\n" \
  "  " pst_integrate_build_system_ARGS_INFO(hintsWeight) "\n" \
  "\n" \
  "  At each level of the hierarchy, the integration consists in solving the" \
  " corresponding linear equation" \
  " system.  " pst_integrate_iterative_INFO(maxIter,convTol,sortSys) "\n" \
  "\n" \
  "   In general, the initial guess for this method at level {r} is" \
  " obtained by scaling down the maps {G,H,Z,R} by 1/2 in all three" \
  " axes, computing the height field {Z'} of level {r+1} recursively" \
  " from this scaled data, and un-scaling {Z'} to obtain an initial" \
  " height map for level {r}.   The iteration at the highest (smallest" \
  " map) level starts with a height map" \
  " may be either all zeros, or the reduced hints map {H}, or the" \
  " reduced reference" \
  " height map {R}, depending on {" initialOpt "}.  In either" \
  " case, a specified amount of noise is added to the initial guess."
    
#endif
