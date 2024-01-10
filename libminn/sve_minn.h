#ifndef sve_minn_H
#define sve_minn_H

/* Quadratic minimzation by the simplex vertex-edge method. */
/* Last edited on 2024-01-10 13:30:53 by stolfi */

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
  
  The number of such /nodes/ is {K(n) = (n+1)*(n+2)/2}. */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <sign.h>

void sve_minn_step(int32_t n, double Fv[], double cm[], bool_t debug);
  /* Given the values {Fv[0..K(n)-1]} of a quadratic function {F} at the 
    nodes of some {n}-simplex {V}, returns in {cm[0..n]} 
    the barycentric coordinates of the stationary point of {F} in the 
    affine subspace spanned by {V}.

    The procedure assumes that the array {Fv} has {K(n)} elements, and
    that {Fv[i*(i+1)/2+j] == F(V(i,j))} for all {i,j} such that 
    {0 <= j <= i <= n}. */

typedef double sve_goal_t(int32_t n, double x[]);
  /* The type of a procedure that can be provided as argument to
    {sve_sample_function} and {sve_optimize} below. It should compute
    some function of the vector {x[0..n-1]}. The function had better be
    C2 near the optimum for fast convergence. */
    
void sve_sample_function(int32_t n, sve_goal_t *F, double v[], double Fv[]);
  /* Evaluates the {n}-variate goal function {F} at the nodes of an
    {n}-simplex {V} in {R^n}, ad stores the values into
    {Fv[0..K(n)-1]}, in the order expected by {sve_minn_step}.
    Assumes that {v} has {K(n)*n} elements and contains the 
    coordinates of the vertices of the simplex, stored by rows.
    
    More precisely, sets {Fv[i*(i+1)/2+j]} to {F(V(i,j))} for all
    {i,j} such that {0 <= j <= i <= n}. */

typedef bool_t sve_pred_t(int32_t n, double x[], double Fx);
  /* The type of a procedure that can be provided as argument to
    {sve_minn_iterate} below. It should check the current solution
    {x[0..n-1]} and the corresponding goal function value {Fx}, and
    return {TRUE} to stop the iteration, {FALSE} to continue it. */
    
void sve_minn_iterate
  ( int32_t n, 
    sve_goal_t *F, 
    sve_pred_t *OK,
    double x[],
    double *FxP,
    sign_t dir,
    double ctr[],
    double dMax,
    bool_t box,
    double rIni,
    double rMin, 
    double rMax,
    double stop,
    int32_t maxIters,
    bool_t debug
  );
  /*  Tries to find a stationary point {x[0..n-1]} of the {n}-argument
    function {F}, by repeated calls to {sve_minn_step}.
    
    Upon entry, the vector {x[0..n-1]} should contain the initial guess,
    and {*FxP} should contain its function value {F(n,x)}. Upon exit,
    {x} will contain the final guess, and {*FxP} with contain the value
    of {F(n,x)}. At each iteration, the procedure chooses a random probe
    simplex centered on the current guess {x[0..n-1]}. The radius {r} of
    the simplex is dynamically adjusted, starting with {rIni} but
    staying within the range {[rMin _ rMax]}.
    
    The parameter {dir} specifies the minimization drection. If {dir ==
    +1}, the procedure looks for a local maximum of {F}. If {dir == -1},
    it looks for a local minimum. If {dir == 0}, it looks for any
    stationary point of {F}, which may be a local minimum, a local
    maximum, or a saddle point. In the latter case, the value of {F} at
    the current guess {x} may increase or decrease during the search;
    especially if the function is neither concave nor convex, or has a
    significant non-quadratic behavior in the region searched. In that
    case, the final guess {x} may be neither the minimum nor the maximum
    of all sample points.
    
    If {dMax} is {+INF}, the search domain is the whole of {\RR}, and
    the {ctr} and {box} parameters are ignored. Otherwise, if {box} is {TRUE}, the
    search is limited to the signed unit cube {ctr + [-dMax _ +dMax]^n}. If {box}
    is {FALSE}, the search domain is the unit {n}-ball {{ x\in \RR^n :
    |x - ctr| <= dMax }}.  In any case, the initial guess {x} had better
    be inside the domain.
    
    The iterations will stop when (A) the current guess {x} satisfies
    the predicate {OK(n,x,F(n,x))}; or (B) the distance between
    successive guesses is less than {stop}; or (C) the quadratic
    minimization loop has been performed {maxIters} times.
    
    If {debug} is TRUE, the procedure prints various diagnostic messages. */

/* LIMITS */

#define sve_minn_MIN_RADIUS (1.0e-100)
  /* The minimum acceptable radius. */

#endif
