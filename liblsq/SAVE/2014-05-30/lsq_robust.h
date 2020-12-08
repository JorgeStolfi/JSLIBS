#ifndef lsq_robust_H
#define lsq_robust_H

/* Fits a linear map of {R^nx} to {R} by least squares, given sample arrays. */
/* Last edited on 2014-05-30 12:09:49 by stolfilocal */

#define lsq_robust_H_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <bool.h>

typedef void lsq_robust_report_t(int iter, int nt, int nx, double Pc[], double Fc[], double Uc[]);
  /* Type of a procedure that is used by {lsq_robust_fit} to report the progress of the
    iteration.

    The argument {iter} will be the mumber of adjustment iterations
    already done.   When {iter} is zero, {Fc[0..nt-1]} will be copy of the original dependent variable samples {F[0..nt-1]},
    and {Pc} will be null.  When {iter} is positive,
    {Fc[0..nt-1]} will be {F[0..nt-1]} adjusted so as to reduce the
    influence of outliers; and {Pc[0..nt-1]} will be the assumed probability of each sample in {F} being an
    inlier, based on the previous approximation.  In either case, {Uc[0..nx-1]}} will be the coefficients of the
    current approximation, computed from {Fc}, as in {lin_fit}. */

void lsq_robust_fit
  ( int nt,       /* Number of cases. */
    int nx,       /* Number of independent variables. */
    double X[],   /* Values of independent variables ({nt} by {nx}). */
    double F[],   /* Values of dependent variable ({nt} elements). */
    double W[],   /* Corresponding weights ({nt} elements). */
    int maxiter,  /* Max iteration count. */
    double U[],   /* (OUT) Fitted linear function coeffs ({nx} elements). */
    double P[],   /* (OUT) Assumed inlier probabilities, or NULL. */
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
    may be one of two kinds, /inlier/ or /outlier/. The value {F[k]} of
    an inlier point is assumed to be the final approximation {(X*U)[k]},
    plus an error with a Gaussian probability distribution. An outlier
    is assumed to have a broader Gaussian distribution, independent of
    the current approximation.

    Before each corrective iteration, the procedure estimates the
    parameters (average and standard deviation) of the inlier and
    outlier distributions, based on the current approximation
    coefficients {U}, and (re)computes the vector {P[0..nt-1]}, setting
    {P[i]} to the estimated probability that data point {i} is an
    inlier. The coefficient vector {U} is then recomputed by weighted
    least squares, as if the actual value {Fc[i]} of each datum was
    {F[i]} if {P[i]} is 1, the current approximation value {Z[i]} if
    {P[i]} is 0, and the linear interpolation of the two for other
    values of {P[i]}

    ??? order of arguments of {report} should be order of computation. ???

    If {report} is not null, at each iteration the procedure calls
    {report(iter,U,Fc,P)} before each correcting iteration.

    !!! Must prove that it converges. !!!

    If {P} is null, the procedure allocates it internally and discards
    it before returning. */

/* AUXILIARY PROCEDURES */

void lsq_robust_compute_stats
  ( int nt,
    double Y[],   /* Values to analyze (deviations or function values, {nt} elements). */
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
    {0..nt-1}, the procedure assumes that {W[k]} is the relevance weight
    of {Y[k]}, and {P[k]} is the current probability of that data point being
    an inlier.

    When {inlier} is true, each {Y[k]} should be the residual of the
    dependent variable (the sampled value minus the fitted model value).
    The effective weight of point {k} in the least squares fitting will
    be {W[k]*P[k]}. Although the mean {*avgP} is computed normally, the
    deviation {*devP} is computed assuming that the mean is zero (as it
    should if {Y} are the deviations from plain least squares).

    ??? Check whether the mean should be considered. ???

    When {inlier} is false, each {Y[k]} must be a value of the dependent
    variable, as given. The effective weight of data point {k} will be
    {W[k]*(1-P[k])}. In this case the computed deviation {*devP} will be
    relative to the computed mean {*avgP}.

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

    Assumes that {F} is the measured value of the dependent variable
    for some data point, and that {pri_gud} is the `a priori'
    probability of that being an inlier.

    The inliers and outliers are assumed to have Gaussian distributions
    with means {avg_gud,avg_bad} and deviations {dev_gud,dev_bad},
    respectively.

    Note that for inlier values the expected value {avg_gud} must include
    the fitted formula value for the same independent variables as {F}. */

void lsq_robust_debug_distr(FILE *wr, char *tag, double avg, double dev, double pri);
  /* Prints to {wr} the parameters of the distribution identified by
    {tag} -- either "gud" (residuals of inliers) or "bad" (values of
    outliers). Namely, the average {avg} and the deviation {dev} of the Gaussian
    distribution, and the a priori probability {pri} of a data point being of that
    class. */

#endif
