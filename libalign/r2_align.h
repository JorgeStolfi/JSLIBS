#ifndef r2_align_H
#define r2_align_H

/* General tools and concepts for translational alignment of 2D objects. */
/* Last edited on 2023-04-01 03:21:02 by stolfi */ 

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
  
/* 
  THE SEARCH DOMAIN
  
  The main functions in this library will search for an optimal
  alignment vector in a /seacrh domain/ {\RD}, an ellipsoid in {\RC}
  defined by two alignment vectors {ctr[0..ni-1]} and {arad[0..ni-1]}.
  The domain consists of all alignment vectors {p} in {\RC} such that the
  difference {v = p - ctr} satisfies these properties:
  
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
    
  The balancing constrait for a given {j} is always omitted if there is
  only one {i} such that {arad[i].c[j]} is non-zero. */
  
int32_t r2_align_count_degrees_of_freedom(int32_t ni, r2_t arad[], bool_t bal);
  /* Returns the the dimension of the search domain {\RD}. 
    This is the number of coordinates {arad[i].c[j]}
    that are nonzero, minus the number of applicable 
    balancing constraints.  */
    
/* 
  THE SEARCH BASIS AND RADII
  
  Let {nd} be the dimension of the seach ellipsoid {\RD}.
  The search for the optimum alignment vector in {\RD} 
  entails the consytuction of a /normal search basis/ {U}, 
  consisting of {nd} alignment vectors {u[0..nd-1]}, and a set
  of {nd} /normal radii/ {urad[0..nd-1]}, such that each vector {u[k]}:
  
    * is balanced and conformal,
    
    * is /normalized/, meaning that the sum of squares of the coordinates
    {u[k][i].c[j]} over all {i,j} is 1,
    
    * is /orthogonal/ to every previous vector {u[r]}, meaning that 
    the sum {u[k][i].c[j]*u[r][i].c[j]} over all {i,j} is zero.
    
    * is a major direction of the ellipsoid {\RD}
    
  Furthermore, each {urad[k]} is the radius of the ellipsoid {\RD} in the 
  direction {u[k]}. 
  
  Thus the search ellipsoid {\RD} is the set of all alignment vectors
  {p = ctr + SUM{ k \in 0..nd-1 : s[k]*urad[k]*u[k]}} where {s}
  is a vector in the unit ball of {\RR^nd}. */

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
    that that together describe the search ellipsoid {\RD},
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

void r2_align_print_vector (FILE *wr, int32_t ni, char *name, int32_t ix, r2_t p[], bool_t ctvars);
  /* Prints the alignment vector {p[0..ni-1]} on {wr}, one point per line. 
     Prints "{name}[{ix}][{i}]" = " before each {p[i]}, or just "{name}[i] = " if {ix}
     is snegative. 
     
     If {ctvars} is true, assumes that {p} is the {arad} ellipsoid size vector, 
     and also prints the number {nz} of zero coordinates and the number {nv=2*ni-2-nz} 
     of degrees of freedom of a conformal balanced delta vector. */

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
