#ifndef lsq_array_H
#define lsq_array_H

/* Fits a linear map of {R^nx} to {R^nf} by least squares, given sample arrays. */
/* Last edited on 2014-05-25 16:08:06 by stolfilocal */

#define lsq_array_H_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#include <bool.h>

int lsq_array_fit
  ( int nt,     /* Number of cases to generate. */
    int nx,     /* Number of independent variables. */
    int nf,     /* Number of dependent variables (functions to fit). */
    double X[], /* Sample values of the independent variables ({nt} by {nx}). */
    double F[], /* Corresponding measured values of the dependent variables ({nt} by {nf}). */
    double W[], /* Corresponding weights ({nt} elements). */
    double U[], /* (OUT) Fitted linear transformation matrix ({nx} by {nt}). */
    bool_t verbose
  );
  /* Finds the linear function {s} from {R^nx} to {R^nf} that best
    approximates some function {f} from {R^nx} to {R^nf}, over a given
    set of cases (data records), by weighted least squares.
    
    Assumes that the data is given as two matrices {X,F} and a vector
    {W} (instead of a generating procedure {gen_case}, as in {lsq_fit}
    from {lsq.h}). The matrices {X} and {F} should have {nt} rows and
    {nx} or {nf} columns, respectively, linearized by rows; and the
    vector {W} should have {nt} elements. Each datum is therefore {xi[j]
    = X[i*nx + j]} for {j} in {0..nx-1}, {fi[k] = F[i*nf + k]} for {k}
    in {0..nf-1}, and {wi = W[i]}. If {W} is null, assumes unit weights 
    for all arguments.
     
    The matrix {U} should have {nx*nf} elements, interpreted as {nx}
    rows and {nf} columns. The approximation is defined as the
    vector-matrix product {S = X U}. Namely, column {k} of {U} is the
    coefficient vector of the linear combination of {xi[0..nx-1]} that
    best approximates coordinate {k} of the function {fi}, in the sense
    of minimizing the sum of squared errors {|X[i,*] U[*,k] - F[i,k]|^2}
    over all data records {i} in {0..nt-1}.
    
    Returns the rank {r} of the least-squares system. If {r} is less
    than {nx}, it means that the rows of {X} span a proper linear
    subspace of {R^{nx}}. In that case the best-fit map is not unique.
   
    This method uses internal working storage proportional {nx*(nx+nf)},
    like {lsq_fit}.
  */

/* AUXILIARY PROCEDURES */

void lsq_array_compute_matrix(int nt, int nx, double X[], double W[], double A[]);
  /* Assumes that {A} is an {nx} by {nx} matrix linearized by rows.
    Fills it with the linear system's matrix for the weighted least squares 
    problem with {nx} independent variables, {nt} data records, argument
    values matrix {X} ({nt} by {nx} elements), and weights vector {W} ({nt} elements).
    See {lsq_array_fit}. for details.
    
    Namely, sets {A} tp {X'*diag(W)*X} where {X} denotes matrix trasposition; that is,
    sets {A[i*nx+j]} to {SUM{ W[k]*X[k*nx+i]*X[k*nx+j] : k \in 0..nt-1 }}
    for all {i,j} in {0..nx-1}. */

void lsq_array_compute_rhs(int nt, int nx, int nf, double X[], double F[], double W[], double B[]);
  /* Assumes that {B} is an array with {nx} rows and {nf} columns.
    Fills it with the right-hand-side of the linear system for the weighted least squares 
    problem with {nx} independent variables, {nt} data records, argument
    values matrix {X} ({nt} by {nx} elements), function values matrix {F} ({nt} by {nf} elements),
    and weights vector {W} ({nt} elements).  See {lsq_array_fit}. for details.
    
    Namely, sets {B} to {X'*diag(W)*F}; that is, set {B[i*nf+j]} to
    {SUM{ W[k]*X[k*nx+i]*F[k*nf+j] : k \in 0..nt-1 }} for all {i} in {0..nx-1}
    and all {j} in {0..nf-1}. */

#endif
