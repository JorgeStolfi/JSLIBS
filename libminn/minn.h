#ifndef minn_H
#define minn_H

/* Basic defs for n-dimensional minimization. */
/* Last edited on 2024-12-05 13:11:09 by stolfi */ 

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

/* 
  THE SEARCH DOMAIN
  
  The functions in this library search the minimum of a function in 
  a /search domain/ {\RD} that is either the
  {n}-dimensional signed unit cube {\RB = [-1 _ +1]^n} or the 
  {n}-dimensional unit ball {\RE} = {{v \in \RR^n : |v| <= 1}}. */

typedef double minn_goal_t (uint32_t n, double c[]); 
  /* Type of a function of {\RR^n} to {\RR} that is to be minimized. */

typedef enum
  { minn_method_ENUM, /* Exhaustive enumeration on grid. */
    minn_method_QUAD  /* Quadratic approx by the SVE method. */
  } minn_method_t;
  /* Minimization method to use. */

void minn_uniform
  ( uint32_t n,          /* Dimension of search space. */
    minn_goal_t *F,     /* Function to be minimized. */
    bool_t box,         /* True to search in the unit cube, false in the unit ball. */
    double atol[],      /* Desired precision along each coordinate. */
    minn_method_t meth, /* Minimizaton method do use.*/
    double v[],         /* (OUT) Minimum vector found. */
    double *Fval_P      /* (OUT) Goal function value at the minimum. */
  );
  /* Minimizes {F} over a search domain {\RD} of {\RR^n}d defined by {dMax} and {box},
    unig the mextod {meth}.
    
    If {dMax} is {+INF}, the search domain {\RD} is the whole of {\RR^n}, and {box}
    is ignored.  If {dMax} is finite, it must be positive; then the domain {\RD} is the 
    signed unit cube {[-dMax _ +dMax]^n} if {box} is true, or the ball of radius {dMax},
    {{v\in \RR^n : |v| <= dMax}} if {box} is false.
    
    The parameter {atol} must be an {n}-vector of positive numbers. 
    The procedure will attempt to achieve precision {atol[i]} along
    each coordinate axis {i}. 
    
    If {meth} is {minn_method_ENUM}, the procedure will enumerate a grid
    of sampling point in [\RD}, with steps {atol[i]} along each coordinate axis.
    In this case the radius {dMax} must be finite. Note that the number of
    probe points will be on the order of {PROD{dMax/atol[i] : i \in 0..n-1}}.
    
    If {meth} is {minn_method_QUAD}, the procedure uses quadratic minimization
    with minimum simplex radius {rMin} equal to the minimum of 
    {atol[k]} for {k} over {0..n-1}. */

void minn_subspace
  ( uint32_t n,          /* Dimension of search space. */
    minn_goal_t *F,     /* Function to be minimized. */
    uint32_t d,          /* Dimension of search domain. */
    double U[],         /* Main axis directions of the search domain. */
    double urad[],      /* Radii of the search domain. */
    bool_t box,         /* True the search domain is a box, false it is an ellipsoid. */
    double utol[],      /* Desired precision along each {U} row. */
    minn_method_t meth, /* Minimizaton method do use.*/
    double v[],         /* (OUT) Minimum vector found. */
    double *Fval_P      /* (OUT) Goal function value at the minimum. */
  );
  /* Minimizes {F} in the search domain {\RD} defined by the orthonormal 
    basis matrix {U}, the radius vector {urad}, and the {box} parameter.
    
    The matrix {U} must have {d} orthonormal rows and {n} columns, and
    {urad} must be a {d}-vector of finite positive numbers. If {U} is {NULL},
    {d} must be {n}, and {U} is assumed to be the identity matrix.
    
    If {box} is true, the search domain is the {d}-dimensional box {\RB(U,urad)}
    in {\RR^n} whose sides are parallel to the rows of {U}, and whose half-width 
    along each direction {U[k,*]} is {urad[k]}.  Namely, {\RB(U,urad)}
    is the set {{ s*U : s \in \RA(urad)}} where {\RA(urad)} is the axis-aligned
    box of {\RR^d} with hal-width {urad[k]} along each axis {k}.
    In particular, if {d} is {n} and {U} is the identity matrix of {NULL},
    then {\RB(U,urad)} is the axis-aligned box {\RA(urad)}.
   
    If {box} is false, the search domain is the {d}-dimensional ellipsoid {\RF(U,urad)}
    in {\RR^n} whose main directions are the rows of {U}, and whose radius
    along each direction {U[k,*]} is {urad[k]}. See {rmxn_ellipsoid.h}.
    IN particular, if {d} is {n} and {U} is the identity matrix of {NULL},
    then {\RF(U,urad)} is the axis-aligned ellipsoid {\RE(urad)}.
    
    The parameter {utol} must be a {d}-vector of positive numbers. The
    procedure will attempt to achieve precision {utol[k]} along each
    ellipsoid axis {U[k,*]}.
    
    If {meth} is {minn_method_ENUM}, the procedure will enumerate a
    {d}-dimensional grid of sampling vectors {v} aligned with the rows
    of {U}, with steps {utol[k]} along each row, and return the minimum
    found among those vectors.
    
    If {meth} is {minn_method_QUAD}, the procedure maps the ellipsoid
    {\RF{U,urad}} to the unit ball of {\RR^d}, then uses quadratic
    minimization on that ball with minimum simplex radius {rMin} equal
    to the minimum of {utol[k]/urad[k]} for {k} over {0..d-1}, and
    finally un-maps the minimum vector {u} of {\RR^d} found that way to
    the corresponding vector {v} of {\RR^n}. */

void minn_ellipsoid_constrained
  ( uint32_t n,          /* Dimension of search space. */
    minn_goal_t *F,     /* Function to be minimized. */
    double arad[],      /* Readii of the base ellipsoid. */
    uint32_t q,          /* Number of explicit constraints. */
    double A[],         /* Constraint matrix. */
    double tol,         /* Desired precision. */
    minn_method_t meth, /* Minimizaton method do use.*/
    double v[],         /* (OUT) Minimum vector found. */
    double *Fval_P      /* (OUT) Goal function value at the minimum. */
  );
  /* Minimizes {F} in the axis-aligned ellipsoid {\RE(arad)} defined by
    the non-negative radii {arad[0..n-1]} (see {rmxn_ellipsoid.h}),
    subject to the linear constraints in the matrix {A}, assumed to be
    {q} by {n} and stored by rows.
   
    Each row {A[k,*]} of {A} is intepreted as the coefficient vector
    of a linear constraint {v*A[k,*]' = 0}.  In addition, for each {i}
    such that {arad[i]} is zero, there is an implied constraint {v[i] = 0}.
    
    The procedure normalizes these constraints with
    {rmxn_ellipsoid_normalize_constraints}, then computes the {d} by {n}
    main directions matrix {U} and the radii {urad[0..d-1]} of the
    {d}-dimensional ellipsoid {\RF(U,urad)} which is the subset of the ellipsoid
    {\RE(arad)} that satisfies all those contraints.
    
    The parameter {tol} must be positive. The procedure will then
    perform {minn_ellipsoid(n,F,d,U,urad,utol,meth,v,Fval_P)} 
    where {utol[0..d-1]} is a vector whose elements are all
    equal to {tol}. This way, the procedure will attempt to achieve
    precision {tol} along each coordinate axis of {\RR^n}. */

#endif
