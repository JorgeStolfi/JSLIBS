/* Generic univariate minimizers */
/* Last edited on 2009-02-10 17:00:20 by stolfi */

#ifndef minu_gen_H
#define minu_gen_H

#include <bool.h>

typedef bool_t (*CheckProc)
  ( void *parms,
    double lo, double hi, 
    double x, double fx, double Dfx,double Dq, double DDfx,double DDqq
  );
  
typedef bool_t (*EvalProc)(void *parms, double x, double *fx, double *dfx);

typedef void (*MinimizeProc)(
    void *parms,     /* Parameters for the `eval' and `check' function. */
    EvalProc eval,   /* Function to minimize. */
    double *x,       /* In: starting guess; out: local min */
    double *fx,      /* In/out: function value at `x'. */
    double *dfx,     /* In/out: derivative at `x'. */
    double tol,      /* Desired uncertainty of result. */
    double dist,     /* Estimated signed distance from minimum */
    double *a,       /* Left endpoint of domain; may be `-INF'. */
    double *b,       /* Right endpoint of domain; may be `+INF'. */
    CheckProc check, /* Acceptance criterion. */
    bool_t debug       /* TRUE to print debugging info */
  );
/*
  Signature for a generic procedure that minimizes an arbitrary
  function of a single real variable, with or without using derivative
  information.  See `minu_brent.h' and `minu_js.h' for specific algorithms.
  
  Modula-3 version created by J.Stolfi on 07/1995; converted to C on 07/2002.

  A generic univariate minimizer looks for an approximate local
  minimum of a function `F' in an interval `[*a__*b]'.

  Upon entry `*x' must contain the initial guess for the minimum
  position (a finite value in `[*a__*b]'); `*fx' must contain the
  corresponding function value, and `*dfx' the corresponding
  derivative.

  The method will call `eval(parms, u, fu, dfu)' to evaluate the
  function and (if necessary) its derivative at additional points `u'
  in `[*a__*b]'.

  Upon exit, `*fx' will contain the minimum function value seen so
  far, `*x' will be the last point where that value occurred, and
  `*dfx' will be the corresponding derivative.

  The variables `*a' and `*b' will be updated to a ``min-range'' for
  this new `*x': a sub-interval `[a1__b1]' of the original `[*a__*b]'
  that contains `*x' and is guaranteed to contain at least one of the
  local minima of `F' in the original interval `[*a__*b]'.

  If `check' is not NULL, the minimizer will call `check(parms, a1,
  b1, *x, *fx, Dfx,Dq, DDfx,DDq)' whenever a new point `*x' is found
  where the function value is less than or equal to the previous
  values. At that point `[a1__b1]' will be a min-range for `*x'. In
  particular, the method will call `check(parms, *a, *b, *x, *fx,
  Dfx,Dq, DDfx,DDq)' right upon entry. For the meaning of the
  parameters `Dfx',`Dq',`DDfx', and `DDq', see below.

  The minimizer procedure will stop searching when it finds a min-range
  of size `tol' or less; or when either `eval' or `check' return
  `TRUE'. More precisely, if `eval' returns `TRUE', the point given to
  it will be discarded; if `check' returns `TRUE', the point given to
  it will be the returned as the best solution.

  The `dist' estimate is used to choose the initial probes around `*x';
  the details depend on the particular instance of the minimizer.

  Derivatives
  -----------

  Some specific types of minimizer ignore the derivative information,
  and work only with the functon values. When using such minimizers,
  the input value of the `*dfx' parameter is irrelevant, the `eval'
  function is not required to compute the derivative, and the output
  value of `*dfx' is meaningless.

  The `check' procedure is given also the first and second derivatives
  `Dfx/Dq' and `DDfx/DDq' of `F' at `*x'; the latter being estimated
  numerically, by divided differences. (If the minimizer does not use
  derivative information, then `Dfx/Dq' will be estimated, too.) The
  scaling factors `Dq' and `DDq' are never negative; but may be zero,
  meaning that the correponding derivative could not be estimated reliably. */

void ComputeDerivatives(
    double u, double fu, double v, double  fv, double x, double fx,
    double *Dfx, double *Dq, double *DDfx, double *DDq
  );
/*
  Computes estimates `Dfx/Dq' and `DDfx/DDq' for the first and second
  derivatives of the goal function `F' at `x', by fitting 
  a quadratic through the three points `(u, fu)', `(v, fv)',
  and `(x, fx)'.  

  The scaling factor `q' will be positive if all three arguments
  `u', `v', `w' are sufficiently distinct, and zero otherwise. */

#endif
