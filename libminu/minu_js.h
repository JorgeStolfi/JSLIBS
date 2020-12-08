// An ad-hoc algorithm for minimizing a univariate function without derivatives.
// M3 version created by J.Stolfi on May-Jul 1994; converted to C on Jul 2002.
// Last edited on 2002-07-15 13:18:17 by stolfi

#ifndef minu_js_H
#define minu_js_H

#include <minu_gen.h>

void minu_js_minimize
  (
    void *parms,
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
  A univariate function minimizer.  See the `MinimizeProc' type
  in `minu_gen.h' for the meaning of the parameters. The `parms'
  argument is ignored.
  
  This procedure does not use the derivative provided
  by `eval', and works only with the function values.  It looks
  for three points `u', `z', `v' in the given interval `[a__b]'
  such that

    1. The points are ordered and not too close:
        `u + delta <= z <= v - delta'
    3. They bracket at lest one local minimum:
        `F(z) <= MIN(F(u), F(v))'
    4. The gaps can't be split any further:
        `MAX(z-u, v-z) < 2*delta'

  where `delta = tol + sqrt(machinePrecision)*ABS(z)'

  When such a triple is found, or when either `eval' or `check'
  return `TRUE', `minimize' returns in `fx' the smallest of
  `F(u)' `F(z)', and `F(w)', and in `x' the corresponding
  argument. */

#endif 

