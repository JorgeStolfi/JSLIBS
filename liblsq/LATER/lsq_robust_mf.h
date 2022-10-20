#ifndef lsq_robust_mf_H
#define lsq_robust_mf_H

/* Fits a linear map of {R^nx} to {R^nf} by least squares, given sample arrays. */
/* Last edited on 2022-10-20 06:31:27 by stolfi */

#define lsq_robust_mf_H_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <bool.h>

typedef void lsq_robust_mf_report_t(int32_t iter, double U[], double Fc[], double P[]);
  /* Type of a procedure that is used by {lsq_robust_mf_fit} to report the progress of the 
    iteration. 
    
    The argument {iter} will be the mumber of adjustment iterations
    done. the argument {U} will be the coefficient matrix of the
    current approximation (an array with {nx} rows and {nf} columns,
    linearized by rows, as in {lin_fit}). The argument {Fc} will be an
    {nt} by {nf} array containing a copy of {F}, adjusted so as to reduce the
    influence of outliers. If {it>0}, the argument {P} is a vector with {nt} elements,
    containing the estimate probability of each datum being an
    outlier. */

void lsq_robust_mf_fit
  ( int32_t nx,       /* Number of independent variables. */
    int32_t nf,       /* Number of dependent variables (functions to fit). */
    int32_t nt,       /* Number of cases to generate. */
    double X[],   /* Values of independent variables ({nt} by {nx}). */
    double F[],   /* Values of dependent variables ({nt} by {nf}). */
    double W[],   /* Corresponding weights ({nt} elements). */
    int32_t maxiter,  /* Max iteration count. */
    double U[],   /* (OUT) Fitted linear transformation matrix ({nx} by {nf]). */
    double P[],   /* (OUT) Assumed outlier probabilities, or NULL ({nt} elements). */
    lsq_robust_mf_report_t *report,
    bool_t verbose
  );
  /* Similar to {lsq_array_fit} from {lsq_array.h}, but after an initial
    least squares fit, applies an iterative method that (hopefully)
    identifies exceptionally bad data points, and reduces their effect.
    
    The parameter {maxiter} is the number of adjustment iterations; if
    zero, the procedure is equivalent to {lsq_fit_arrays}.
    
    The procedure assumes that each data point {xi,fi} (each pair of
    rows from {X} and {F}) may be one of two kinds, /inlier/ or
    /outlier/. An inlier point is assumed to be the current
    approximation {U * xi}, plus an error with a multivariate Gaussian
    probability distribution. An outlier is assumed to have a
    broader distribution, independent of the current approximation.
    
    Before each corrective iteration, the procedure computes a vector
    {P[0..nt-1]} such that {P[i]} is the estimated probability that
    data point {i} is an outlier. It also computes the average and
    covariance matrices of the inlier errors and of the outlier values,
    based on the currently fitted matrix {U}. The matrix {U} is then
    recomputed as if the actual weight of each datum {i} were
    {W[i]*(1-P[i])}.  !!! Check whether {P} is inlier or outlier. !!! 
    
    If {report} is not null, at each iteration the procedure calls
    {report(iter,U,Fc,P)} before each correcting iteration.
    
    !!! Must prove that it converges. !!!
    
    !!! Describe {K,L} !!!
    
    If {K}, {L}, and/or {P} are null, the procedure allocates them
    internally and discards them before returning. */

/* AUXILIARY PROCEDURES */

void lsq_robust_mf_compute_stats
  ( int32_t nf, 
    int32_t nt,
    double Y[],   /* Values to analyze (deviations or function values, {nt} by {nx}). */
    double W[],   /* Weights of data records ({nt} elements). */
    double P[],   /* Probabilities
    bool_t good, 
    double Vref,
    double E[], 
    double K[],
    double L[],
    double *priP,
    bool_t verbose
  );
  /* Computes the averages and covariances of {nt} vectors of 
    sample values of {nt} variables, for inliers ({good=TRUE}) or outliers
    ({good=FALSE}), and estimated `a priori' probability of a data record
    being good or bad, respectively.

    The array {Y} must have {nt} rows and {nf} columns, linearized by
    rows. The vectors {W} and {P} must have {nt} elements. For each {k}
    in {0..nt-1}, and each {j} in {0..nf-1}, the procedure assumes that
    {Y[nf*k + j]} is the value of variable {j} associated to data record {k},
    {W[k]} is the relevance weight of that data record, and {P[k]} is the
    current probability of that data record being a good one (`inlier').
    ??? Or is it `outlier'? ???
    
    When {good} is true, the {Y} matrix should contain the residuals of
    the independent variables (the sampled values minus the fitted model
    values). The effective weight of data record {k} will be
    {W[k]*P[k]}.
    
    When {good} is false, the {Y} matrix shoudl contain the
    sampled values of the independent variables, and the effective
    weight of data record {k} will be {W[k]*(1-P[k])}.
    
    The argument {Verf} is used only when {good} is false, and must be the largest 
    principal variance of the inlier distribution.  It is used to ensure that 
    the outlier distribution is wider than the inlier one. 
    
    If {W} is null, assumes {w[k]==1.0} for all {k}.
    If {P} is null, assumes {p[k]==0.5} for all {k}.
    
    The array {E} must have {nf} elements.  The procedure reeturns in {E[j]} the average 
    value of variable {j} (that is, {Y[k*nf + j]}).  
    
    The procedure returns also in {K}, which must be an {nf} by {nf} matrix, 
    the eigenvectors of the covariance matrix of the variables 
    specified in the {Y} matrix.  These eigenvectors are the principal directions of the joint distribution
    of those variables. The variances of those variables along those principal directions are returned in the 
    {L} vector}, which must be an {nf}-element vector.  If {K}, {L}, or {P} are null, they are allocated
    internally and then discarded.
    
    Finally, the procedure returns 
    in {*priP} (if not null) the estimated probability of 
    a data record belonging to the selectd class, which is
    the weighted mean of {P[0..nt-1]} or its complement, depending on {good}. */

double lsq_robust_mf_bayes
  ( int32_t nf, 
    double F[], 
    double P_gud,
    double E_gud[], 
    double K_gud[],
    double L_gud[],
    double E_bad[], 
    double K_bad[],
    double L_bad[],
    bool_t verbose
  );
  /* Computes the `a posteriori' probability of a data record being a good sample  (inlier),
    using Bayes's formula.
    
    Assumes that {F[0..nf-1]} are the measured values of the dependent
    variables for some data record, and that {P_gud} is the `a priori'
    probability of the data record being good.
    
    The good data records (inliers) and the bad samples (outliers) are
    assumed to have multivariate normal distributions. The arrays
    {E_gud,E_bad} must have {nf} elements, the expected values of the
    {F} variables for inliers and outliers, respectively; {K_gud,K_bad}
    must be {nf} by {nf} matrices, whose rows are the orthonormal
    eigenvector bases (principal directions) extracted from the
    covariance matrices; and {L_gud,L_bad} must be the corresponding
    eigenvalues (variances along the principal axes).
    
    Note that for good values the expected values {E_gud} must include
    the fitted formula values for the same independent variables as {F}.
    The procedure assumes that, a priori, the data record is good with
    probability {P_gud}. */

void lsq_robust_mf_fudge_covariance_matrix(int32_t nf, double alpha, double V[], double beta, double Q[]);
  /* Assumes that {V} and {Q} are {nf} by {nf} covariance matrices
    of {nf} random variables. Scales the covariance matrix {V} by the
    factor {alpha}, then adds a very small multiple of the identity to
    make sure that it is invertible. If {Q} is not null, also makes
    sure that the diagonal of {V} is at least {alpha} times the
    diagonal of {Q}.
    
    !!! Rethink this procedure. !!! */

void lsq_robust_mf_debug_distr(FILE *wr, char *tag, int32_t nf, double E[], double V[], double pri);
  /* Prints to {wr} the parameters of the distribution identified by
    {tag} -- either "gud" (residuals of inliers) or "bad" (values of
    outliers). The parameter {E} is a vector of {nf} numbers, and
    {V} is an {nf} by {nf} matrix stored by rows. 
    
    The printout has {nf} lines; the first number in each line {k} is
    the computed mean value {E[k]} of the distribution, and the next
    {nf} numbers are row {k} of {V}. The procedure also prints the
    computed a priori probability {pri} of a data record being of that
    class. */

#endif
