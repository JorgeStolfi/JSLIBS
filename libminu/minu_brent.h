// Richard P. Brent's algorithm for minimizing a univariate function w/o derivatives.
// Converted from FORTRAN 77 to Modula-3 by J.Stolfi
// Last edited on 2009-02-08 21:13:19 by stolfi

#ifndef minu_brent_H
#define minu_brent_H

#include <minu_gen.h>

/* !!! Remove the {*parms} argument; can use nested procs instead. !!! */

void minu_brent_minimize
  ( void *parms,
    EvalProc eval,    
    double *x,        
    double *fx,       
    double *dfx,      
    double tol,       
    double dist,      
    double *a,        
    double *b,        
    CheckProc check,  
    bool_t debug       
  );
/*
  A reimplementation of Richard P. Brent's univariate function
  minimizer. See the `MinimizeProc' type in `minu_gen.h' for the
  meaning of the parameters. The `parms' argument is ignored.
  
  These comments were extracted from the FORTRAN 77 implementation:
  
  The `minimize' method does not use the derivative `dfx' returned by
  `eval'; it works only with the function values `fx'.  It uses a
  combination of golden section search and successive parabolic
  interpolation.  
  
  Convergence is never much slower than that for a Fibonacci
  search.  If `F' has a continuous second derivative which is
  positive at the minimum (which is not at `a' or `b'), then
  convergence is superlinear, and usually of the order of about
  1.324....
  
  The function `F' is never evaluated at two points closer together
  than `eps*abs(x)+(tol/3)', where `eps' is approximately the square
  root of the relative machine precision.  If `F' is a unimodal
  function and the computed values of `F' are always unimodal when
  separated by at least `eps*abs(x)+(tol/3)', then `x' approximates
  the abcissa of the global minimum of `F' on the interval `(a,b)' with
  an error less than `3*eps*abs(x)+tol'.  If `F' is not unimodal,
  then `x' may approximate a local, but perhaps non-global, minimum
  to the same accuracy.
  
  This algorithm is a somewhat rehashed version of the Algol 60
  procedure `localmin' given in Richard Brent, `Algorithms for
  minimization without derivatives', Prentice-Hall, Inc. (1973). */

#endif
