#ifndef r2_align_H
#define r2_align_H

/* Tools for lists of points of the plane describing translational alignment of 2D objects. */
/* Last edited on 2021-12-17 15:25:14 by stolfi */ 

#include <bool.h>
#include <r2.h>
#include <i2.h>

/* 
  ALIGNMENT VECTOR
  
  A (two-dimensional) /alignment vector/ for {ni} objects is a list of
  points {p[0..ni-1]} of the plane, that is, a vector of {(\RR^2)^ni}.
  It says that point {p[i]} of each object {i} corresponds in some sense
  to object {p[j]} of some other object {j}.
  
  For example, the objects could be images, and the alignment {p[0..ni-1]} could be 
  claiming that some neighborhood of point {p[i]} of image {i} looks like
  the similar neighborhood of point {p[j]} of image {j}.
  
  BALANCED SEARCH
  
  An alignment vector {p} is /balanced/ with respect to a reference
  alignment {p0} if the differences {p[i].c[j] - p0[i].c[j]} add to
  zero along each coordinate {j}.  */

double r2_align_norm_sqr(int32_t ni, r2_t p[]);
  /* Given an alignment vector {p} with {ni} components, returns the square of its total
    Euclidean norm. Namely, returns the sum of squares of all {p[i].c[j]}
    for {i} in {0..ni-1]} and {j} in {0..1}. */

double r2_align_dot(int32_t ni, r2_t p[], r2_t q[]);
  /* Returns the dot product of the alignment vectors {p} and {q} with {ni}
    components each. Namely, returns the sum of products {p[i].c[j]*q[i].c[j]}
    for {i} in {0..ni-1]} and {j} in {0..1}. */

double r2_align_dist_sqr(int32_t ni, r2_t p[], r2_t q[]);
  /* Given two alignment vectors {p,q} with {ni} components, returns the square of its total
    Euclidean distance. Namely, returns the sum of squares of {p[i].c[j] - q[i].c[j]}
    for {i} in {0..ni-1]} and {j} in {0..1}. */

double r2_align_rel_dist_sqr(int32_t ni, r2_t p[], r2_t q[], r2_t arad[]);
  /* Given two alignment vectors {p,q} with {ni} components, returns the total squared 
    displacement between them, relative to the radius vector {arad}, that is,
    
      { SUM { (p[i].c[j] - q[i].c[j])/arad[i].c[j]|^2 : i=0..ni-1,j=0..1 } }
    
    for all {i} in {0..ni-1} and {j} in {0..1}; except that terms where
    {arad[i].c[j]} is zero are ignored. 
    
    When minimizing some mismatch function within the ellipsoid with
    center {q} and radius {arad}, clients may want to add an appropriate
    multiple of this function to the goal function in order to bias the
    search towards the neighborhood of {q}, in case the function is
    mostly constant in that region. */

void r2_align_throw_ball_vector(int32_t ni, r2_t p[], double rmin, double rmax);
  /* Fills {p[0..ni-1]} with an alignment vector uniformly distributed in the
    ball centered at the origin, whose norm is in {[rmin _ rmax]}. */

void r2_align_throw_ortho_disp_vectors(int32_t ni, int32_t nd, r2_t U[], r2_t arad[]);
  /* Stores into {U} a set of {nd} balanced displacements for alignment vectors 
    with {ni} points each.
    
    Specifically, for each {k} in {0..nd-1}, the displacement vector
    {u[k]} will be points {u[k][i] = U[k*ni + i]}, for {i} in {0..ni-1}.
    That vector will have the following properties:
    
      * It is balanced in both axes, meaning that for each {j} in
      {0..1}, the coordinates {u[k][i].c[j]} summed over all {i}, add to
      zero.
        
      * It is normalized, meaning that the squares of the coordinates
      {u[k][i].c[j]}, summed over all {i,j}, add to 1.
      
      * It is orthogonal to every other displacement vector {u[r]} with
      {r!=k}, maeaing that the products {u[k1][i].c[j] * u[k2][i].c[j]},
      summed over all {i,j}, add to zero.
      
      * It is conformal to {arad}, in the sense that {u[h][i].c[j]} is
      zero whenever {arad[i].c[j]} is zero.
      
    If {arad[o].c[j]} is non-zero, its value must be positive, but it is otherwise
    ignored by the procedure.
    
    The parameter {ni} must be at least 1, and the number of directions
    {nd} must be at most {nv-2}, where {nv} is the number of coordinates
    {arad[i].c[j]} that are not zero. */

/* VARIABLE AND FIXED COODRINATES */

bool_t r2_align_coord_is_variable(r2_t arad[], int32_t i, int32_t j);
  /* {TRUE} if coordinate {j} of point {p[i]} is variable, {FALSE} if it is fixed.
    Namely, it returns {TRUE} if {arad[i].c[j]} is 1.0 pr more. */

int32_t r2_align_count_variable_coords(int32_t ni, r2_t arad[]);
  /* Returns the number of coordinates in the alignment vectors that are variable,
    as per {coord_is_variable(arad, i, j)} for {i} in {0..ni-1} and {j} in {0..1}. */
  
void r2_align_points_to_vars(int32_t ni, r2_t p[], r2_t arad[], r2_t (1,1)[], r2_t p0[], int32_t nv, double y[]);
  /* Stores the active displacements {p[0..ni-1][0..1]-p0[0..ni-1][0..1]} into {y[0..nv-1]}. */

void r2_align_vars_to_points(int32_t nv, double y[], int32_t ni, r2_t arad[], r2_t (1,1)[], r2_t p0[], r2_t p[]);
  /* Stores {y[0..nv-1]} into the active {p[0..ni-1][0..1]}, adding {p0[0..ni-1][0..1]}. */

/* MISMATCH FUNCTIONS */

typedef double r2_align_mismatch_t(int32_t ni, r2_t p[], i2_t iscale); 
  /* Type of a function that evaluates the mismatch of {ni} objects.
    in some neighborhood of points {p[0..ni-1]}, at a specified
    scale of resolution {iscale}.  The size of the neighborhood
    is client-defined, but is assumed to be proportional
    to {2^iscale.c[j]} along each axis {j}.
    
    Nevertheless, the function should ideally take about the same time
    whatever the {iscale}. That is, each object should be sampled with a
    fixed number of sampling points, spread over an interval
    proportional to {2^iscale.c[j]} along each axis {j}.
    
    The size of the neighborhood, measured on the scaled object, should
    be the same. That is, measured on the unscaled object, the
    neighborhod should be centered on {p[i]} and have size.
    
    The client must ensure that the sampling is reliable, e.g. by
    interpolating, smoothing, or shrinking the object according to
    Nyquist's criterion.
    
    Said another way, the mismatch should be computed as if every object
    is expanded or reduced by the scale factor {1/2^iscale.c[j]} along
    axis {j}, with a suitable smoothing filter, with each alignment
    point coordinate {p[i].c[j]} implicitly scaled by that same amount;
    and then sampled at a fixed number of points over a neighborhood of
    fixed size, independent of {iscale}.
    
    The function had better be C2 near the optimum alignment for
    quadratic optimization (see below). */

/* DEBUGGING/TESTING */

void r2_align_plot_mismatch
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    i2_t iscale,
    r2_align_mismatch_t *f2p
  );
  /* Writes to file {wr} a 2D plot data for the function {*f2p}, assumed
    defined over the ellipsoid of {(\RR^2)^ni} with center
    {ctr[0..ni-1].c[0..1]} and semi-axes {arad[0..ni-1].c[0..1]}.
    
    Namely, chooses two orthonormal balanced displacement vectors {u,v}
    in {(\RR^2)^ni} and enumerates a grid of points {p = ctr + ku*stu*u
    + kv*stv*v} with suitable range for th eintegers {ku,kv} and
    suitable steps {stu,stv}. For each sample point that lies inside the
    ellipsoid, writes to {wr} a line with {ku*stu}, {kv*stv}, and
    {f2p(ni,p,iscale)}. */
    
#endif
