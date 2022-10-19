#ifndef r2_align_H
#define r2_align_H

/* General tools and concepts for translational alignment of 2D objects. */
/* Last edited on 2021-12-19 09:33:39 by stolfi */ 

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>

/* 
  ALIGNMENT VECTOR
  
  A (two-dimensional) /alignment vector/ for {ni} objects is a list of
  points {p[0..ni-1]} of the plane, that is, a vector of {\RC = (\RR^2)^ni}.
  It says that point {p[i]} of each object {i} corresponds in some sense
  to object {p[j]} of some other object {j}.
  
  For example, the objects could be images, and the alignment {p[0..ni-1]} could be 
  claiming that some neighborhood of point {p[i]} of image {i} looks like
  the similar neighborhood of point {p[j]} of image {j}.
  
  The number of objects {ni} is assumed fixed in the rest of these comments.
  A pair {i,j} with {i} in {0..ni-1} and {j} in {0..1} identifies a coordinate
  axis of {\RC} and a coordinate {p[i].c[j]} of any alignment vector {p}.
  
  We denote by {\BX} and {\BY} the alignment vectors such that {\BX[i]} = (1,0)}
  sand {\BY[i] = (0,1)} for all {i} in {0..ni-1}. */

double r2_align_norm_sqr (int32_t ni, r2_t p[]);
  /* Given an alignment vector {p} with {ni} components, returns the square of its total
    Euclidean norm. Namely, returns the sum of squares of all {p[i].c[j]}
    for {i} in {0..ni-1]} and {j} in {0..1}. */

double r2_align_dot (int32_t ni, r2_t p[], r2_t q[]);
  /* Returns the dot product of the alignment vectors {p} and {q} with {ni}
    components each. Namely, returns the sum of products {p[i].c[j]*q[i].c[j]}
    for {i} in {0..ni-1]} and {j} in {0..1}. */

double r2_align_dist_sqr (int32_t ni, r2_t p[], r2_t q[]);
  /* Given two alignment vectors {p,q} with {ni} components, returns the square of its total
    Euclidean distance. Namely, returns the sum of squares of {p[i].c[j] - q[i].c[j]}
    for {i} in {0..ni-1]} and {j} in {0..1}. */

void r2_align_throw_ball_vector (int32_t ni, double rmin, double rmax, r2_t p[]);
  /* Fills {p[0..ni-1]} with an alignment vector uniformly distributed in the
    ball centered at the origin, whose norm is in {[rmin _ rmax]}. */

/* 
  BASIC DOMAIN ELLIPSOID
  
  Searches for an optimal alignment vector are confined to a /basic
  ellipsoid/ {\RE \sub \RC}, defined by two alignment vectors
  {ctr[0..ni-1]} (the /domain center/) and {arad[0..ni-1]} (the /radius
  vector/). The coordinates {arad[i].c[j]} of the radius vector must be
  non-negative. In the rest of these comments, the vectors {ctr} and
  {arad} are assumed fixed.
  
  The ellipsoid {\RE} has the main axes aligned with the coordinate axes
  of {\RC}, and its radius along the axis {i,j} is {arad[i].c[j]}.
  Namely, it consists of all alignment vectors {p} such that the sum of
  the squares of {(p[i].c[j] - ctr[i].c[j])/arad[i].c[j]}, over all
  {i,j} such that {arad[i].c[j]!=0}, is at most 1. */

double r2_align_rel_dist_sqr (int32_t ni, r2_t p[], r2_t q[], r2_t arad[]);
  /* Given two alignment vectors {p,q} with {ni} components, returns the
    total squared adjustment between them, relative to the radius
    vector {arad}. 
    
    If {q} is NULL, assumes it is the alignment vector with all zeros.
    
    That is, if {p-q} is not conformal to {arad}, returns {+INF}.
    Otherwise returns the sum of ((p[i].c[j] - q[i].c[j])/arad[i].c[j])^2 
    for all {i} in {0..ni-1} and {j} in {0..1} such that {arad[i].c[j]} is non-zero.
    
    Thus the basic domain ellipsoid {\CE} consists of all alignment
    vectors {p} such that {r2_align_rel_dist_sqr(ni, p, ctr,arad) <= 1}. */

/* 
  CONFORMAL ADJUSTMENT VECTORS
  
  We will denote by {\RB} the affine span of the ellipsoid {\RE}, and by
  {\RV} the linear subspace of {\RC} parallel to {\RB}. Their dimension
  is {nv = 2*ni - nz}, where {nz} is the number of coordinates of {arad}
  that are zero.
  
  Each element of {\RV} is called a/conformal adjustment/. It is an
  alignment vector {d} such that {d[i].c[j]} is zero whenever
  {arad[i].c[j]} is zero. */

bool_t r2_align_coord_is_variable (r2_t arad[], int32_t i, int32_t j);
  /* Returns {TRUE} iff {arad[i].c[j]} is nonzero. 
  
    Namely, returns {TRUE} if the coordinate axis {i,j} of {\RC} is in
    the subspace {\RV}; that is, if coordinate {p[i].c[j]} of a point
    {p} in the basic ellipsoid {\CE} can be different from
    {ctr[i].c[j]}. */

i2_t r2_align_count_variable_coords (int32_t ni, r2_t arad[]);
  /* Returns a integer pair {nv} such that {nv.c[j]} is the number of
    coordinates {arad[i].c[j]} that are not zero.
    
    That is, {nv.c[0]} and {nv.c[1]} are the counts of axes of the
    conformal subspace {\CV} that are parallel to {\BX} and {\BY},
    respectively. */

void r2_align_points_to_vars (int32_t ni, r2_t p[], r2_t arad[], r2_t ctr[], int32_t nv, double y[]);
  /* Stores the non-fixed coordinates of the adjustment from
    {p[0..ni-1]} to {ctr[0..ni-1]}, into the vector {y[0..nv-1]}. 
    
    The number {nv} must be the total number of non-fixed coordinates ,
    that is, the number of non-zero coordinates {arad[i].c[j]}. If {ctr}
    is {NULL}, assumes it is all zeros. */

void r2_align_vars_to_points (int32_t nv, double y[], int32_t ni, r2_t arad[], r2_t ctr[], r2_t p[]);
  /* Stores {ctr[0..ni-1]} plus {y[0..nv-1]} into the non-fixed 
    coordinates of the alignment vector {p[0..ni-1]}. 
    
    Namely sets {p[i].c[j]} to {ctr[i].c[j] + y[k]} whenever
    {arad[i].c[j]} is nonzero. The number {nv} must be the total number
    of non-fixed coordinates , that is, the number of non-zero
    coordinates {arad[i].c[j]}. If {ctr} is {NULL}, assumes it is all
    zeros. */
  
/*   
  BALANCED ADJUSTMENT VECTORS
  
  A /balanced adjustment/ is an alignment vector {d} such that
  the sum of {d[i].c[j]} over all {i} in {0..ni-1} is zero, for each {j}.
  
  That is, a balanced adjustment is an element of {\RC} that is
  orthogonal to the alignment vectors {\BX} and {\BY}.
  
  We will denote by {\RU} the subspace of {\RV} (hence of {\RC})
  consisting of all balanced and conformal adjustments; and by
  {\RA} the affine subspace of {\RC} that contains the center {ctr} and
  is parallel to {\RC}.
  
  THE SEARCH ELLIPSOID
  
  Procedures that look for an optimal alignment vector will consider only
  alignments {p} in {\CE} that such that the difference {p-ctr}
  is a balanced adjustment. Since {p} is in {\CE}, that difference 
  will also be a conformal adjustment.  
  
  Thus those procedures will consider only adjustment vectors
  {p} in the /search ellipsoid/ {\CF} that is the intersection of 
  the basic ellipsoid {\CE} and the affine subspace {\CB}.  */
  
int32_t r2_align_count_degrees_of_freedom (int32_t ni, r2_t arad[]);
  /* Returns the max number of linearly independent balanced adjustment
    vectors that are conformal with the radius vector {arad}.  That is,
    the dimension of the spaces {\RU}, {\RA} and 
    
    Namely, let {nv} be {r2_align_count_variable_coords(ni,arad)}.
    Each axis {j} contributes {nv.c[j]-1} degrees of freedom,
    unless {nv.c[j]} is zero, in which case it contributes none. */

void r2_align_compute_search_ellipsoid (int32_t ni, r2_t arad[], int32_t nd, r2_t U[], double urad[]);
  /* Stores into {U} a set of {nd} orthonormal adjustment vectors with {ni} points each, and in {urad} 
    a list of {nd} real radii, that together define the search ellipsoid {\RF}.
    
    Specifically, for each {k} in {0..nd-1}, the adjustment vector
    {u[k]} will be points {u[k][i] = U[k*ni + i]}, for {i} in {0..ni-1}.
    That vector will have the following properties:
      
      * It will be conformal to {arad}, in the sense that {u[k][i].c[j]} will be
      zero whenever {arad[i].c[j]} is zero.
    
      * It will be a balanced adjustment, meaning that for each {j} in
      {0..1}, the coordinates {u[k][i].c[j]} summed over all {i}, will add to
      zero.
        
      * It will be normalized, meaning that the squares of the coordinates
      {u[k][i].c[j]}, summed over all {i,j}, will add to 1.
      
      * It will be orthogonal to every other adjustment vector {u[r]} with
      {r!=k}, meaning that the products {u[k][i].c[j] * u[r][i].c[j]},
      summed over all {i,j}, add to zero.
      
      * It will be the direction of a major axis of the search ellipsoid {\CF}.
      
    Moreover, {urad[k]} will be the radius of the ellipsoid {\CF} in the direction {u[k]}.
      
    If {arad[i].c[j]} is non-zero, its value must be positive, but it is otherwise
    ignored by the procedure.
    
    The parameter {ni} must be at least 1, and the number of directions
    {nd} must be the dimension of the balanced conformal adjustment space {\RU}. */

/* MISMATCH FUNCTIONS */

typedef double r2_align_mismatch_t (int32_t ni, r2_t p[]); 
  /* Type of a function that evaluates the mismatch of {ni} objects.
    in some neighborhood of points {p[0..ni-1]}.  The size of the neighborhood
    is client-defined.
    
    The client must ensure that the sampling is reliable, e.g. by
    interpolating, smoothing, or shrinking the object according to
    Nyquist's criterion.
    
    The function had better be C2 near the optimum alignment for
    quadratic optimization (see below). */

/* DEBUGGING/TESTING */

void r2_align_print_vector (FILE *wr, int32_t ni, char *name, int32_t ix, r2_t p[], bool_t ctvars);
  /* Prints the alignment vector {p[0..ni-1]} on {wr}, one point per line. 
     Prints "{name}[{ix}][{i}]" = " before each {p[i]}, or just "{name}[i] = " if {ix}
     is snegative. 
     
     If {ctvars} is true, assumes that {p} is the {arad} ellipsoid size vector, 
     and also prints the number {nz} of zero coordinates and the number {nv=2*ni-2-nz} 
     of degrees of freedom of balanced adjustments. */

void r2_align_plot_mismatch
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    double step,
    r2_align_mismatch_t *F2
  );
  /* Writes to file {wr} a data for a 3D plot of the function {*F2}, assumed
    defined over the ellipsoid of {\RC = (\RR^2)^ni} with center {ctr[0..ni-1]} 
    and semi-axes {arad[0..ni-1]}.
    
    Namely, finds the directions {u0,u1} of the two longest axes of the search ellipsoid {\RF}
    and enumerates a grid of points with step {tol} on the plane defined 
    by them.  For each sample point {(s0,s1)} of that grid that lies inside the
    ellipsoid, computes the corresponding alignment vector {p}
    and writes to {wr} a line with {s0}, {s1}, and {F2(ni,p)}.*/
    
void r2_align_throw_arad (int32_t ni, r2_t zfrac, double rmin, double rmax, r2_t arad[]);
  /* Stores into {arad[0..ni-1]} a random radius vector. For each {j}, approximately {zfrac.c[j]*ni}
    coordinates {arad[0..ni-1].c[j]} will be zero; the other cordinates 
    will be randomly chosen in the range {[rmin _ rmax]}. */

#endif
