#ifndef pst_slope_map_H
#define pst_slope_map_H

/* pst_slope_map.h -- procedures for working with slope maps. */
/* Last edited on 2024-12-22 23:11:41 by stolfi */

#include <bool.h>
#include <r2.h>

#include <float_image.h>
#include <pst_imgsys.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>

/* SLOPE MAPS
  
  A /slope map/ is a two-channel float-valued image of some
  height field {Z(X,Y)}, where the value of each pixel is the 
  gradient of that height --- that is, its derivatives
  {dZ/dX} and {dZ/dY}.
  
  !!! TO DO !!!
  
  !!! Make {pst_slope_map_build_integration_system_1} accept a weight map, too.
  
  !!! Generalize to use a 2-channel weight map, with separate
  reliability weights for the X and Y component of the gradient. 
  
*/

r2_t pst_slope_map_get_pixel(float_image_t *G, int32_t x, int32_t y);
  /* Extracts the derivatives {dZ/dX} and {dZ/dY} from channels 0 and
    1 of the pixel in column {x}, row {y} of {IG}, and returns them
    as a {r2_t} gradient vector. */

void pst_slope_map_set_pixel(float_image_t *G, int32_t x, int32_t y, r2_t *grd);
  /* Stores the derivatives {dZ/dX} and {dZ/dY}, taken from the
    gradient vector {grd} into channels 0 and 1 of the pixel in column
    {x}, row {y} of {IG}. */

void pst_slope_and_weight_map_shrink
  ( float_image_t *IG, 
    float_image_t *IW, 
    float_image_t **SG, 
    float_image_t **SW
  );
  /* Given a slope map {IG}, containing the derivative of a height
    function {Z} along X and Y axes, returns another slope map {SG},
    with half the size as {IG}, containing the derivatives of a
    version {SZ} of {IZ} with both dimensions and heights scaled by
    half.  If the given image has size {NX} by {NY},
    the result has size {NX/2} by {NY/2}, rounded up.
    
    A pixel in column {x} and row {y} of the result is conceptually
    centered at the vertex point {(2x+1,2y+1)} of {IG}'s domain.
    
    If {IW} is not null, it should be a monochromatic image with the
    same size as {IG}. Each element {IW[0,x,y]} (which must be
    non-negative) is interpreted as the relative reliability of the
    slopes {IG[c,x,y]}, and is used as a weight in the averaging
    of the slopes. In particular, pixels of {IG} that have zero reliability in {IW} are ignored in
    the averaging.  If {IW} is not null, the preocedure returns also a 
    reduced version {SW} of {IW}.
    
    The images {SG} and {SW} are allocated by the procedure. */

/* COMPUTING HEIGHTS FROM SLOPES */
           
typedef void pst_slope_map_report_proc_t(uint32_t level, float_image_t *IG, float_image_t *IW); 
  /* Type of a client-given procedure that may be called
    by recursive integrators to report the reduced slope map and weight map used at each scale. */   

void pst_slope_map_to_depth_map_recursive
  ( float_image_t *IG, 
    float_image_t *IW, 
    uint32_t level,
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    float_image_t **OZP,
    float_image_t **OWP,
    bool_t verbose,
    uint32_t reportIter,
    pst_slope_map_report_proc_t *reportData,
    pst_imgsys_report_proc_t *reportSys,
    pst_height_map_report_proc_t *reportHeights
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
    
    If {IW} is not null, the procedure also creates a weight map {OW}
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
    each level, on the way down (decreasing {level}). Here {iter} is
    the number of complete Gauss-Seidel or Gauss-Jacobi iterations
    performed before the call; and {changeZ} is the max absolute
    change in any {OZ} sample since the previous iteration
    (meaningless when {iter=0}). The procedure {reportHeights} is
    always called once with {final=TRUE} after the last iteration. If
    {reportIter>0}, it is also called with {final=FALSE} before each
    iteration whose index {iter} is multiple of {reportIter};
    including before the first iteration, if any. In particular, if
    the {OZ} map is computed directly without iteration,
    {reportHeights} is called just once with {iter=0} and
    {final=TRUE}. Therefore it is never called twice with the same
    {level} and {iter} */

/* THE INTEGRATION SYSTEM */

pst_imgsys_t* pst_slope_map_build_integration_system(float_image_t *G, float_image_t *W, bool_t full);
  /*
    Allocates and returns an {pst_imgsys_t} system {S} whose equations define
    the height map {Z} given a slope map {G}.

    The value of {G[0,x,y]} must be the X derivative {dZ/dX} of the heigh map,
    averaged over the pixel {[x,y]} The value of {G[1,x,y]} must
    be the Y derivative {dZ/dY} averaged over the same pixel.
    
    The samples of the {Z} map are heoghts at the CORNERS of
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
    is the unknown solution vector, and {b} is a known vector. The
    system has {S.N} equations and {S.N} unknowns.
    The size {S.N} is at most {NZ = (NX+1)*(NY+1)}, but may be
    much less, depending on the weight map.
    
    If {full} is FALSE, then any height {Z[x,y]} which is between four
    pixels with zero weight is considered to be indeterminate and is
    excluded from the system. The client must compute those heights
    separately.
    
    If {full} is true, those heights are included too. In this case
    the equation of that height relates it to its four neighbors with
    equal weight, but from it does not enter into the equations of those
    neighbors unless they too are indeterminate.  The effect is to fill
    in any indeterminate regions with a smooth surface that is attached 
    at the edges to the determinate heights.  In this case {S.N}
    will be equal to {NZ}.
    
    In either case the procedure stores in {S.col,S.row,S.ix}
    three tables that tell the correspondence between the unknowns and equations of
    the system and the elements of the {Z} array.  Namely,
    for each {k} in {0..S.N-1}, {unknown h[k]} is the height
    element {Z[S.col[k],S.row[k]]}.  Conversely, for each {x} in {0..NX}
    and each {y} in {0..NY}, the sample {Z[x,y]} corresponds
    to the unknown {h[k]} where {k = S.ix[x + y*(NX+1)]}, if {k >= 0};
    os is omitted from the system, if {k<0}.
    
    If all weights are equal, the equations essentially state that the
    Laplacian of the height, as computed from the (unknown) height map
    {Z}, is equal to the Laplacian computed from the (known) slope
    maps {G}. Both Laplacians are computed by finite-difference
    formulas. Other equations, based on equality of slopes, are used
    along the edges and where the weights are zero.
    
    If all weights are non-zero, the system has a one-dimensional
    family of solutions, since the homogeneous system {A h = 0} 
    is satisfied by {h[k]==1} for all {k}. */
  
void pst_slope_map_solve_system
  ( pst_imgsys_t *S, 
    float_image_t *OZ,
    uint32_t ord[], 
    uint32_t maxIter, 
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose, 
    uint32_t level, 
    uint32_t reportIter, 
    pst_height_map_report_proc_t *reportHeights
  );
  /* 
    Solves system {S} by Gauss-Seidel iteration, and stores the solution
    into the height map {OZ}.  On input {OZ} should contain the initial
    guess. 
    
    The vector {ord[0..S->N-1]} specifies the order in which
    the equations are to be solved.
      
    The image {OZ} must be single-channel. The procedure assumes that the
    unknowns of {S} are related to the samples of {OZ}, as defined by
    the tables {ix[0..NX*NY-1]}, {col[0..S->N-1]}, and {row[0..S->N-1]}.

    If {reportHeights} is not null and {reportIter} is not zero, will call
    {reportZ(level,iter,change,final,OZ)} before the first iteration, after
    the last iteration, and after any number of iterations that is
    divisible by {reportIter}. */

void pst_slope_map_copy_height_map_to_sol_vec(pst_imgsys_t *S, float_image_t *IZ, double VZ[]);      
void pst_slope_map_copy_sol_vec_to_height_map(pst_imgsys_t *S, double VZ[], float_image_t *IZ);
  /* Copies the values of the unknowns {Z[i]} from or to the
    elements of the height map {IZ[0,x,y]}. */

#endif
