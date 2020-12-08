// An ad-hoc algorithm for minimizing a univariate function with derivatives.
// Created by J.Stolfi on Jul 2002.
// Last edited on 2007-11-01 20:19:12 by stolfi

#ifndef minu_herm_H
#define minu_herm_H

#include <minu_gen.h>

void minu_herm_minimize
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
  
  This procedure uses the derivative provided by `eval', and Hermite
  interpolation to estimate the minimum within an interval. It looks for
  two points `u,v' in `[a __ b]' such that

    1. The points are ordered:
        `u <= v'
    2. The points are not too close:
        `v - u >= delta'
    3. They bracket at least one local minimum:
        (`DF(u) <= 0' or `u == a') and (`DF(v) >= 0' or `v == b')
    4. The gap can't be split any further:
        `v - u < 2*delta'

  where `delta = tol + sqrt(machinePrecision)*ABS(z)'

  When such a triple is found, or when either `eval' or `check'
  return `TRUE', `minimize' returns in `fx' the smallest of
  `F(u),F(v)', and in `x' the corresponding argument. */

#endif 

