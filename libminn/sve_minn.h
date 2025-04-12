#ifndef sve_minn_H
#define sve_minn_H

/* Quadratic minimzation by the simplex vertex-edge method. */
/* Last edited on 2025-04-01 09:44:30 by stolfi */

/* SIMPLICES

  An /{n}-simplex/ is a list {V} of {n+1} points of {R^n}, for some
  {n >= 0}. Those points are the /corners/ of the simplex, and any two
  distinct corners define an /edge/. The simplex is /degenerate/ if
  some corner can be written as an affine combination of the other
  corners; and is /proper/ otherwise.
  
  In this interface, {V(i)} denotes corner number {i} of the simplex {V} 
  
  SIMPLEX REPRESENTATION
  
  In this interface, an {n}-simplex {V} in {R^n} is represented
  as an array {v} with {(n+1)*n} elements, conceptually
  {n+1} rows and {n} columns; where row {i} contains the {n}
  coordinates of corner {i} of the simplex. 
  
  SIMPLEX NODES
  
  Let {V(i,j)} denote the midpoint of {V(i)} to {V(j)}. Note
  that {V(i,j)==V(j,i)} for all {i,j} in {0..n}, and also that
  {V(i,j)==V(i)==V(j)} when {i == j}.
  
  The main quadratic minimization routine ({sve_minn-step})
  requires the values of the goal function at all the corners
  and edge midpoints of some {n}-simplex {V}; that is, at 
  all points {V(i,j)} such that {0 <= j <= i <= n}.
  
  The number of such /nodes/ is {nf = (n+1)*(n+2)/2}. */

#include <stdint.h>

#include <bool.h>
#include <sign.h>

void sve_minn_quadratic_optimum(uint32_t n, double Fv[], double cm[], bool_t debug, bool_t debug_system);
  /* Given the values {Fv[0..nf-1]} of a quadratic function {F} at the 
    nodes of some {n}-simplex {V}, returns in {cm[0..n]} 
    the barycentric coordinates of the stationary point of {F} in the 
    affine subspace spanned by {V}.

    The procedure assumes that the array {Fv} has {nf} elements, and
    that {Fv[i*(i+1)/2+j] == F(V(i,j))} for all {i,j} such that 
    {0 <= j <= i <= n}. */

typedef double sve_goal_t(uint32_t n, const double x[]);
  /* The type of a procedure that can be provided as the {F} parameter to
    {sve_sample_function} and {sve_optimize} below. It should compute
    some function of the vector {x[0..n-1]}. The function had better be
    C2 near the optimum for fast convergence. */
    
void sve_sample_function(uint32_t n, sve_goal_t *F, double v[], double Fv[]);
  /* Evaluates the {n}-variate goal function {F} at the nodes of an
    {n}-simplex {V} in {R^n}, ad stores the values into
    {Fv[0..nf-1]}, in the order expected by {sve_minn_quadratic_optimum}.
    Assumes that {v} has {nf*n} elements and contains the 
    coordinates of the vertices of the simplex, stored by rows.
    
    More precisely, sets {Fv[i*(i+1)/2+j]} to {F(V(i,j))} for all
    {i,j} such that {0 <= j <= i <= n}. */

uint32_t sve_minn_single_step
  ( uint32_t n, 
    double x[], double *FxP,
    sve_goal_t *F, sign_t dir,
    double dCtr[], double dMax, bool_t dBox, 
    double *radiusP, double rMin,
    double v[], double Fv[],
    uint32_t *nEvalsP,
    double *dStepP,
    bool_t debug,
    bool_t debug_probes
  ); 
  /* Performs a single SVE quadratic optimization step of a goal
    function {F}, including clipping the guess to a specified
    box or ball domain, choosing the probe simplex, sampling the goal
    function at the simplex nodes, computing the quadratic optimum
    from those samples, clipping it to the domain, evaluating {F}
    at that point, and choosing the best of all
    values thus obtained.
    
    The search domain is assumed to have center {dCtr[0..n-1]} and radius
    {dMax}, being either a box or a ball depending on the {dBox} flag.
    
    On input, {x[0..n-1]} is assumed to be the guess for the optimum 
    (or previously computed optimum), and {*FxP} is the goal function
    value at that point.
    
    The procedure first chooses a probe simplex {v[p..nv*n-1]} with some
    radius {rad} and center {csi[0..n-1]}, where {nv} is {n+1}. If
    possible, the center {csi} will be {x} and the radius will be the
    input value of {*radiusP}, but both may have to be adjusted to
    ensure that the simplex lies within the search domain. The adjusted
    radius (which will not be less than {rMin}) is returned in
    {*radiusP}. The simplex will have {nv = n+1} vertices. Coordinate
    {j} of simplex vertex {i} will be stored in {v[i*n + j]}, for {i} in
    {0..nv-1} and {j} in {0..n-1}.
    
    The procedure then evaluates the function {F} at the
    {nf=(n+1)*(n+2)/2}} probe points, which are the vertices and edge
    midpoints of the simplex, using {sve_sample_function}. The values
    are stored in {Fv[0..nf-1]}.
    
    The procedure then computes the stationary point {y[0..n-1]}
    from those sampling points and sampled values.  If this 
    point lies outside the search domain, it is pulled towards
    the domain center {dCtr} until it lies on the boundary of the domain.
    Then the procedure computes the value {Fy} of the goal function {F}
    at this point.
    
    If {dir} is zero, the procedure just replaces the input vector {x}
    with the computed and clipped stationary point {y}, and sets {*FxP}
    to {Fy}. Otherwise it identifies either the maximum (if {dir} is
    {+1}) or minimum (if {dir} is {-1}) among the input {*FxP}, the
    value {Fy} at the computed stationary point, and the {nf} sampled
    values {Fv[0..nf-1]}.  It then sets {*FxP} to that maximum or minimum, and
    and sets {x} to the corresponding point ({x}, {y} one of the
    vertices in {v}, or one of the edge midpoints).
    
    The procedure increments {*nevals_P} with the number of calls to {F}
    that it made (normally {nf+1}, for the simplex sampling and 
    the stationary poont {y}).
    
    The procedure also sets {*dStepP} to the Euclidean distance
    between the input and output positions of the point {x}. 
    
    The return value will be 0 if the winner is the computed stationary
    point {y}, 1 if it is one of the probe points, and 2 if it is the
    input solution {x} (i.e. if the step failed altogether). */

/* LIMITS */

#define sve_minn_MIN_RADIUS (1.0e-100)
  /* The minimum acceptable radius. */

#endif
