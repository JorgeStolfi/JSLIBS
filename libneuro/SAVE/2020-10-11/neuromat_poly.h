#ifndef neuromat_poly_H
#define neuromat_poly_H

/* NeuroMat filtering and spectral analysis tools. */
/* Last edited on 2013-11-19 02:23:43 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <bool.h>

/* A polynomial {P} of degree {g} is represented by its {g+1} coeffs {P[0..g]} such that
  { P(x) = SUM{x^k*P[k]:k\in 0..g} }. */

double neuromat_poly_eval(int g, double P[], double x);
  /* Evaluates the polynomial with coeffs {P[0..g]} at the given {x} value. */
  
void neuromat_poly_shift(int g, double P[], double a, double Q[]);
  /* Returns in {Q} the polynomial {P} of degree {g} shifted forward by {a}. 
  
    Namely, given the coeffs {P[0..g]} of a polynomial {P(·)}, returns
    the coeffs {Q[0..g]} of the polynomial {Q(·)} such that {Q(x) = P(x
    - a)} for all {x}. The vector {Q} must be either disjoint from {P}
    or the same as {P}. */

void neuromat_poly_stretch(int g, double P[], double h, double Q[]);
  /* Returns in {Q} the polynomial {P} of degree {g} stretched by the factor {h}. 
  
    Namely, given the coeffs {P[0..g]} of a polynomial {P(·)}, returns
    the coeffs {Q[0..g]} of the polynomial {Q(·)} such that {Q(x) =
    P(x/h)} for all {x}. The vector {Q} must be either disjoint from {P}
    or the same as {P}. */

void neuromat_poly_bezier(int g, int i, double a, double b, double P[]);
  /* Sets {P[0..g]} to the coefficients of the Bernstein-Bezier polynomial
    of degree {g} and index {i} for the domain interval {[a_b]}. */

void neuromat_poly_eval_multi(int g, double P[], int n, double x[], double s[]);
  /* Assumes that {x} and [s} have {n} elements. Sets {s[k]} to the
    value of the polynomial with coeffs {P[0..g]} at the argument
    {x[k]}, for {k} in {0..n-1]}. 
    
    If {x} is null, assumes {x[k] = 2*(k + 1/2)/n - 1} for {k} in {0..n-1}. */

void neuromat_poly_fit(int n, double x[], double y[], double w[], int g, double P[]);
  /* Finds the coeffs {P[0..g]} of the polynomial that best approximates the
    values {y[0..n-1]} at the argument values {x[0..n-1]}, in the least squares
    sense. 
    
    If {x} is null, assumes {x[k] = 2*(k + 1/2)/n - 1} for {k} in {0..n-1}.
    If {w} is not null, assumes {w[0..n-1]} are importance weights for
    the data pairs. If {w} is null, assumes the same weight for all data
    pairs. */

typedef void neuromat_poly_report_proc_t(int iter, double R[], double z[]);
  /* Type of procedure used by {neuromat_poly_fit_robust} to report its progress. */

void neuromat_poly_fit_robust
  ( int n, 
    double x[], 
    double y[], 
    double w[], 
    int maxiter, 
    int g, 
    double P[],
    neuromat_poly_report_proc_t *report
  );
  /* Similar to {neuromat_poly_fit}, but uses an iterative method that (hopefully)
    ignores the effect of outliers.  The parameter {maxiter} is the number of
    iterations; if zero, the procedure is equivalent to {neuromat_poly_fit}.
    
    If {report} is not null, at each iteration calls {report(iter,R,z)} where {iter}
    is the iteration count, {z[0..n-1]} is the current adjusted data values,
    and {R[0..g]} is the current fitted polynomial.*/

#endif
