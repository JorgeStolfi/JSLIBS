#ifndef r2_align_H
#define r2_align_H

/* General tools and concepts for translational alignment of 2D objects. */
/* Last edited on 2024-12-05 10:19:59 by stolfi */ 

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
  and {\BY[i] = (0,1)} for all {i} in {0..ni-1}. */

/* MISMATCH FUNCTIONS */

typedef double r2_align_mismatch_t (int32_t ni, r2_t p[]); 
  /* Type of a function that evaluates the mismatch of {ni} objects. in
    some neighborhood of points {p[0..ni-1]}. The size of the
    neighborhood is client-defined.
    
    The client must ensure that the sampling is reliable, e.g. by
    interpolating, smoothing, or shrinking the object according to
    Nyquist's criterion.
    
    The function had better be C2 near the optimum alignment for
    quadratic optimization (see below). */
  
/* 
  THE SEARCH DOMAIN
  
  The main functions in this library will search for an optimal
  alignment vector in a /search domain/ {\RF}, an ellipsoid in {\RC}
  defined by two alignment vectors {ctr[0..ni-1]} and {arad[0..ni-1]},
  and an optional balancing constraint. The domain consists of all alignment
  vectors {p} in {\RC} such that the difference {v = p - ctr} satisfies
  these properties:
  
    * {v} is /conformal/, meaning that {{v[i].c[j]} is zero whenever
    {arad[i].c[j]} is zero.
    
    * {v} is in the axis-aligned ellipsoid with center {ctr} and radii
    {arad}, meaning that  the squares of
    the coordinates {v[i].c[j]/arad[i].c[j]}, summed over all {i,j}
    for which the denominator is nonzero, add to at most 1.
    
  Optionally, the following /balancing constraint/ may be required,
  depending on a boolean parameter {bal}

    * {v} is /balanced/, meaning that for each {j} in
    {0..1}, the coordinates {v[i].c[j]} summed over all {i}, will add to
    zero.
    
  The balancing constraint for a given {j} is always dropped if there is
  only one {i} such that {arad[i].c[j]} is non-zero. */
    
/* 
  THE SEARCH BASIS AND RADII
  
  Let {nd} be the dimension of the seach ellipsoid {\RF}.
  The search for the optimum alignment vector in {\RF} 
  entails the construction of a /normal search basis/ {U}, 
  consisting of {nd} alignment vectors {u[0..nd-1]}, and a set
  of {nd} /normal radii/ {urad[0..nd-1]}, such that each vector {u[k]}:
  
    * is conformal
    
    * is balanced, if {bal} is true,
    
    * is /normalized/, meaning that the sum of squares of the coordinates
    {u[k][i].c[j]} over all {i,j} is 1,
    
    * is /orthogonal/ to every previous vector {u[r]}, meaning that 
    the sum {u[k][i].c[j]*u[r][i].c[j]} over all {i,j} is zero.
    
    * is a major direction of the ellipsoid {\RF}
    
  Furthermore, each {urad[k]} is the radius of the ellipsoid {\RF} in the 
  direction {u[k]}. 
  
  Thus the search ellipsoid {\RF} is the set of all alignment vectors
  {p = ctr + SUM{ k \in 0..nd-1 : s[k]*urad[k]*u[k]}} where {s}
  is a vector in the unit ball of {\RR^nd}. */

i2_t r2_align_count_variable_coords (int32_t ni, r2_t arad[]);
  /* Returns a integer pair {nv} such that {nv.c[j]} is the number of
    indices {i} such that {arad[i].c[j]} is not zero. */
  
int32_t r2_align_count_degrees_of_freedom(i2_t nv, bool_t bal);
  /* Returns the the dimension of the search domain {\RF}, given the number
    {nv.c[j]} of variable coords along each Cartesian axis {j} and the 
    balancing requirement {bal}.  
    
    Namely, if {bal} is false, it is just {nv.c[0]} plus {nv.c[1]}.  
    If {bal} is true, then each {nv.c[j]} is reduced by 1, unless
    it is 0 or 1. */

void r2_align_compute_search_ellipsoid
  ( int32_t ni,
    r2_t arad[],
    bool_t bal, 
    int32_t nd,
    r2_t U[], 
    double urad[]
  );
  /* Stores into the first {nd} rows of {U} a normal search basis 
    and into {urad[0..nd-1]} the corresponding normal radii 
    that that together describe the search ellipsoid {\RF},
    as defined by {arad[0..ni-1]} and {bal}.
      
    The radii {arad[i].c[j]} must be non-negative. The 
    parameter {nd} must be as computted by 
    {r2_align_count_degrees_of_freedom}. */

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

void r2_align_print_vector (FILE *wr, int32_t ni, char *name, int32_t ix, r2_t p[]);
  /* Prints the alignment vector {p[0..ni-1]} on {wr}, one point per line. 
     Prints "{name}[{ix}][{i}]" = " before each {p[i]}, or just "{name}[i] = " if {ix}
     is negative. */

void r2_align_plot_mismatch
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    bool_t bal,
    double step,
    r2_align_mismatch_t *F2
  );
  /* Writes to file {wr} a data for a 3D plot of the function {*F2}, assumed
    defined over the ellipsoid of {\RC = (\RR^2)^ni} with center {ctr[0..ni-1]} 
    and semi-axes {arad[0..ni-1]}.
    
    Namely, finds the directions {u0,u1} of the two longest axes of the
    search ellipsoid {\RF} determined by {ctr,arad,bal}, and enumerates
    a grid of points with step {tol} on the plane defined by them. For
    each sample point {(s0,s1)} of that grid that lies inside the
    ellipsoid, computes the corresponding alignment vector {p} and
    writes to {wr} a line with {s0}, {s1}, and {F2(ni,p)}.*/
            
void r2_align_throw_vector (int32_t ni, double rmin, double rmax, r2_t p[]);
  /* Fills {p[0..ni-1]} with an alignment vector uniformly distributed in the
    hollow ball of {\RC} centered at the origin with inner radius {rmin}
    and outer radius {rmax}. */
    
void r2_align_throw_ctr (int32_t ni, double cmax, r2_t ctr[], bool_t verbose);
  /* Stores into {arad[0..ni-1]} a random vector suitable for the center
    of a test. It will be uniformly distributed in the ball of {\RC}
    centered at the origin with radius {cmax}. If {verbose} is true also
    prints it to {stderr}. */
    
void r2_align_throw_arad (int32_t ni, double rmax, int32_t nvmin, r2_t arad[], bool_t verbose);
  /* Stores into {arad[0..ni-1]} a random radius vector. For each {j},
    some fraction of the coordinates {arad[0..ni-1].c[j]} will
    be zero; the other cordinates will be randomly chosen in the range
    {[0.5*rmax _ rmax)}.  However, at least {min(ni,nvmin)}
    will be nonzero.  If {verbose} is true, also prints it to {stderr}. */
            
double r2_align_dot (int32_t ni, r2_t p[], r2_t q[]);
  /* Returns the dot product of the alignment vectors {p} and {q} with {ni}
    components each. Namely, returns the sum of products {p[i].c[j]*q[i].c[j]}
    for {i} in {0..ni-1]} and {j} in {0..1}. */

double r2_align_dist_sqr (int32_t ni, r2_t p[], r2_t q[]);
  /* Given two alignment vectors {p,q} with {ni} components, returns the square of its total
    Euclidean distance. Namely, returns the sum of squares of {p[i].c[j] - q[i].c[j]}
    for {i} in {0..ni-1]} and {j} in {0..1}. */

double r2_align_rel_disp_sqr (int32_t ni, r2_t p[], r2_t q[], r2_t arad[]);
  /* Given two alignment vectors {p,q} with {ni} components, returns the
    total squared coordinate differences between them, relative to the radius
    vector {arad}. 
    
    If {q} is NULL, assumes it is the alignment vector with all zeros.
    
    That is, if {p-q} is not conformal to {arad}, returns {+INF}.
    Otherwise returns the sum of ((p[i].c[j] - q[i].c[j])/arad[i].c[j])^2 
    for all {i} in {0..ni-1} and {j} in {0..1} such that {arad[i].c[j]} is non-zero.
    
    Thus the basic domain ellipsoid {\CE} consists of all alignment
    vectors {p} such that {r2_align_rel_disp_sqr(ni, p, ctr,arad) <= 1}. */

#endif
