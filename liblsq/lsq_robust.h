#ifndef lsq_robust_H
#define lsq_robust_H

/* Fits a linear map of {R^nx} to {R} to samples by least squares, ignoring outliers. */
/* Last edited on 2019-12-18 16:06:46 by jstolfi */

#define lsq_robust_H_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <bool.h>

typedef void lsq_robust_report_t(int32_t iter, int32_t nt, int32_t nx, double Pc[], double Fc[], double Uc[]);
  /* Type of a procedure that is used by {lsq_robust_fit} to report the progress of the
    iteration.

    The argument {iter} will be the mumber of adjustment iterations
    already done. When {iter} is zero, {Fc[0..nt-1]} will be copy of the
    original function samples {F[0..nt-1]}, and {Pc} will be
    null. When {iter} is positive, {Fc[0..nt-1]} will be {F[0..nt-1]}
    adjusted so as to reduce the influence of outliers; and
    {Pc[0..nt-1]} will be the assumed probability of each sample in {F}
    being an inlier, based on the previous approximation. In either
    case, {Uc[0..nx-1]}} will be the coefficients of the current
    approximation, computed from {Fc}, as in {lin_fit}. */

void lsq_robust_fit
  ( int32_t nt,          /* Number of data points. */
    int32_t nx,          /* Number of independent variables (argument coordinates per data point). */
    double X[],      /* Argument coordiantes of data points ({nt} by {nx}). */
    double F[],      /* Function samples of data points ({nt} elements). */
    double W[],      /* Reliability weights ({nt} elements). */
    int32_t maxiter,     /* Max iteration count. */
    double U[],      /* (OUT) Fitted linear function coeffs ({nx} elements). */
    double P[],      /* (OUT) Assumed inlier probabilities, or NULL. */
    lsq_robust_report_t *report,
    bool_t verbose
  );
  /* Similar to {lsq_array_fit} from {lsq_array.h}, with {nf} set to 1;
    except that, after an initial least squares fit, it applies an
    iterative method that (hopefully) identifies exceptionally bad data
    points, and reduces their effect.

    The parameter {maxiter} is the number of adjustment iterations; if
    zero, the procedure is equivalent to {lsq_fit_arrays}.

    The procedure assumes that each data point {X[k*nx+{0..nx-1}],F[k]},
    may be one of two kinds, /inlier/ or /outlier/. 
    
    The function sample {F[k]} of an inlier point is assumed to be the true linear
    function of the arguments, plus an error with a Gaussian probability
    distribution that has zero mean and deviation {dev_gud/sqrt(W[k])},
    for some unknown parameter {dev_gud}.
    
    An outlier is assumed to have a broader Gaussian distribution, whose
    mean {avg_bad} is unknown but independent of the current approximation
    and of the index {k}, and whose deviation is {dev_bad/Wk}, for some
    unnown {dev_bad}.

    The procedure first computes an initial coefficient matrix {U} with
    plain weighted least squares, and with it computes an initial
    approximation value {S[k] = (X*U)[k]}} for each data point, and sets
    {P[k] = 0.5}. Then it performs up to {maxiter} corrective
    iterations. At each iteration, the procedure estimates the
    parameters {dev_gud,avg_bad,dev_bad} of the inlier and outlier
    distributions, based on the current values of {S,F,P}. With those
    parameters and the values {F[k]} and {S[k]}, it (re)computes {P[k]},
    the estimated probability that data point {k} is an inlier. Finally
    it recomputes the coefficient vector {U} by weighted least squares,
    but substituting the function sample {F[k]} by an adjusted sample
    {Fc[k]} that is {F[k]} if {P[k]} is 1, {S[k]} if {P[k]} is 0, and
    the linear interpolation of the two for other values of {P[k]}

    If {report} is not null, at each iteration the procedure calls
    {report(iter,P,Fc,U)} before each correcting iteration.
    
    Note that only the relative values of the weights are relevant. A
    data point with zero weight is ignored. Moreover, for any integer
    {N}, a data point that has weight {Wk} counts the same as {N} data
    points with same arguments and function samples, but with weight
    {Wk/N}.

    !!! Must prove that it converges. !!!

    If {P} is null, the procedure allocates it internally and discards
    it before returning. */

/* AUXILIARY PROCEDURES */

void lsq_robust_compute_stats
  ( int32_t nt,
    double Y[],   /* Values to analyze (residuals or function samples, {nt} elements). */
    double W[],
    double P[],
    bool_t inlier,
    double dev_ref,
    double *avgP,
    double *devP,
    double *priP,
    bool_t verbose
  );
  /* Computes the average {*avgP} and standard deviation {*devP} of {nt}
    values {Y[0..nt-1]}, for inliers ({inlier=TRUE})
    or outliers ({inlier=FALSE}). Also returns in {*priP} the estimated
    `a priori' probability of data point {k} being inlier or
    outlier, respectively.

    The vectors {W} and {P} must have {nt} elements. For each {k} in
    {0..nt-1}, the procedure assumes that {W[k]} is the original
    relevance weight of data point {k}, and {P[k]} is the current
    probability of that data point being an inlier.

    When {inlier} is true, each {Y[k]} should be the current residual of the
    data point {k} (the original function sample {F[k]} minus the
    current fitted model value {S[k]}). The effective weight of point {k} in the
    least squares fitting will be {W[k]*P[k]}. Although the mean {*avgP}
    is computed normally, the deviation {*devP} is computed assuming
    that the mean is zero (as it should if {Y} are the residuals from
    plain weighted least squares).

    ??? Check whether the mean should be considered. ???

    When {inlier} is false, each {Y[k]} must be the original sample
    function sample {F[k]}. The effective weight of data point {k} will
    be {W[k]*(1-P[k])}. In this case the computed deviation {*devP} will
    be relative to the computed mean {*avgP}.

    The argument {dev_ref} is used only when {inlier} is false, and must be the
    deviation of the inlier distribution.  It is used to ensure that
    the outlier distribution is wider than the inlier one.

    If {W} is null, assumes {w[k]==1.0} for all {k}.
    If {P} is null, assumes {p[k]==0.5} for all {k}. */

double lsq_robust_bayes
  ( double F,
    double pri_gud,
    double avg_gud,
    double dev_gud,
    double avg_bad,
    double dev_bad,
    bool_t verbose
  );
  /* Computes the `a posteriori' probability of a data point being a
    inlier, using Bayes's formula.

    Assumes that {F} is the function sample for some data point, and
    that {pri_gud} is the `a priori' probability of that being an
    inlier.

    The inliers and outliers are assumed to have Gaussian distributions
    with means {avg_gud,avg_bad} and deviations {dev_gud,dev_bad},
    respectively. */

void lsq_robust_debug_distr(FILE *wr, char *vname, char *tag, double avg, double dev, double pri);
  /* Prints to {wr} the parameters of the distribution identified by
    {vname} and {tag} -- either "samples" and "outlier", or "residuals" and "inlier".
    Namely, the average {avg} and the deviation {dev} of the Gaussian
    distribution, and the a priori probability {pri} of a data point being of that
    class. */

#endif
