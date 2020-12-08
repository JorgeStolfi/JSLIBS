/* Univariate Bézier polinomials. */
/* Last edited on 2014-04-04 01:30:27 by stolfilocal */

#ifndef bz_basic_H
#define bz_basic_H

#include <stdint.h>

/* 
  BÉZIER REPRESENTATION OF A POLYNOMIAL
  
  A univariate polynomial of degree {g} can be specified by its
  /Bézier coefficients/ relative to some interval {B}. These are an
  array of {g+1} real values, conceptually located at a regularly
  spaced set points spanning {B}.
  
  BERNSTEIN-BÉZIER POLYNOMIALS
  
  For a fixed degree {g}, each Bézier coefficient is identified
  by a /Bézier index/ {e} in the range {0,.. g}. The
  coefficient with index {e} multiplies the
  /Bernstein-Bézier polynomial/ of degree {g} and index {j}
  
    { BB^g_e(z) = choose(g,e) z^{g-e} (1-z)^e } */
  
 
typedef int8_t bz_degree_t;
  /* Degree of a polynomial. */

typedef bz_degree_t bz_index_t;  
  /* Identifies a Bézier coeff along some axis. */

double bz_bernstein(bz_degree_t g, bz_index_t i, double z);
  /* Evaluates the Bernstein-Bézier polynomial of degree {g} and index
   {i} for the argument {z}. The index {i} must lie in {0..g}.  */

double bz_bernstein_max(bz_degree_t g, bz_index_t i);
  /* The maximum value of {bz_bernstein(g,i,z)} for {z} in [0_1]. */

void bz_split
  ( bz_degree_t g, /* Degree of curve. */
    double c[],    /* Bezier coeffs of a polynomial piece. */
    double wa,     /* Width of first half. */
    double a[],    /* OUT: Bézier coeffs of first half. */
    double wb,     /* Width of second half. */
    double b[]     /* OUT: Bézier coeffs of second half. */
  );
  /* Splits a univariate Bernstein-Bézier polynomial piece in two
    pieces by the DeCasteljau algorithm. 
    
    More precisely, given the Bézier coeffs {c[0..g]} of a polynomial
    {P} relative to some interval {I}, computes the Bézier coeffs
    {a[0..g]} and {b[0..g]} of {P} relative the two pieces {Ia,Ib} of {I}
    whose widths are {wa} and {wb}.
    
    The pointers {a} and/or {b} may be null, in which case the
    corresponding coeffs are not computed; otherwise {a[0..g]} and
    {b[0..g]} must be disjoint. The vector {c} may be identical to {a}
    or to {b}, in which case it is overwritten with the corresponding
    result; otherwise {c[0..g]} must be disjoint from both. */

void bz_diff
  ( bz_degree_t g,  /* Degree of curve. */
    double c[],     /* Bezier coeffs of a polynomial. */
    double w,       /* Actual width of interval. */
    double d[]      /* OUT: Bezier coeffs of derivative. */
  );
  /* Computes the derivative of a polynomial in Bézier form.
    
    More precisely, given the Bézier coeffs {c[0..g]} of a polynomial
    {P} relative to an interval of width {w}, returns in {d[0..g-1]} the
    coeffs of the derivative {P}, relative to the same interval. The vectors
    {c} and {d} need not be disjoint. */

void bz_integ
  ( bz_degree_t g,  /* Degree of curve. */
    double c[],     /* Bezier coeffs of a polynomial. */
    double w,       /* Actual width of interval. */
    double d[]      /* OUT: Bezier coeffs of integral. */
  );
  /* Computes the integral of a polynomial in Bézier form.
    
    More precisely, given the Bézier coeffs {c[0..g]} of a polynomial
    {P}, relative to an interval of width {w}, returns in {d[0..g+1]}
    the coeffs of {Q}, relative to the same interval; where {Q(x)} is
    the integral of {P} from the low end of the interval to {x}. The
    vectors {c} and {d} need not be disjoint. */

void bz_eval
  ( bz_degree_t g,   /* Degree of curve. */
    double c[],      /* Bezier coeffs. */
    double u,        /* Start of interval. */
    double v,        /* End of interval */
    double x,        /* Argument value. */
    bz_degree_t ord, /* Max derivative order desired. */
    double f[]       /* OUT: Image vector (size {ord+1}). */
  );
  /* Given the Bézier coeffs {c[0..g]} of a polynomial {P} relative to
    the interval {[u _ v]}, stores in {f[0..ord]} the value of {P} and
    all derivatives up to order {ord} at the argument {x}. The value
    is stored in {f[0]}, and the derivative of order {k} in {f[k]} */

#endif
