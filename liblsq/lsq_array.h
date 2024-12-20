#ifndef lsq_array_H
#define lsq_array_H

/* Fits a linear map of {R^nx} to {R^nf} by least squares, given sample arrays. */
/* Last edited on 2024-12-05 12:55:00 by stolfi */

#define lsq_array_H_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#include <stdint.h>
#include <bool.h>

uint32_t lsq_array_fit
  ( uint32_t nt,     /* Number of data points. */
    uint32_t nx,     /* Number of independent variables (argument coordinates per data point). */
    uint32_t nf,     /* Number of dependent variables (function samples per data point). */
    double X[],     /* Argument coordinates for all data points ({nt} by {nx}). */
    double F[],     /* Corresponding function samples ({nt} by {nf}). */
    double W[],     /* Corresponding reliability weights ({nt} elements). */
    double U[],     /* (OUT) Fitted linear transformation matrix ({nx} by {nt}). */
    bool_t verbose
  );
  /* Finds the linear function {s} from {R^nx} to {R^nf} that best
    approximates some function {f} from {R^nx} to {R^nf}, over a given
    set of data points, by weighted least squares.
    
    Assumes that the data points are given as two matrices {X,F} and a
    vector {W} (instead of a generating procedure {gen_data_point}, as
    in {lsq_fit} from {lsq.h}). The matrices {X} and {F} should have
    {nt} rows and {nx} or {nf} columns, respectively, linearized by
    rows; and the vector {W} should have {nt} elements, interpreted as a
    column. Each data point is therefore {Xk[i] = X[k*nx + i]} for {i}
    in {0..nx-1}, {Fk[j] = F[k*nf + j]} for {j} in {0..nf-1}, and {Wk =
    W[k]}, for {k} in {0..nt-1}. If {W} is null, assumes unit weights
    for all arguments.
    
    The procedure stores into {U} the coefficient matrix of the 
    fitted linear approximation {s}.
    The matrix {U} should have {nx*nf} elements, interpreted as {nx}
    rows and {nf} columns. The approximation is defined as the
    vector-matrix product {S = X U}. Namely, column {j} of {U} is the
    coefficient vector of the linear combination of {Xk[0..nx-1]} that
    best approximates coordinate {j} of the function {Fk}, in the sense
    of minimizing the sum of squared errors {Wk*|X[k,*] U[*,j] - F[k,j]|^2}
    over all data points {k} in {0..nt-1}.
    
    Note that only the relative values of the weights are relevant. A
    data point with zero weight is ignored. Moreover, for any integer
    {N}, a data point that has weight {Wk} counts the same as {N} data
    points with same arguments and function samples, but with weight
    {Wk/N}.
    
    Returns the rank {r} of the least-squares system. If {r} is less
    than {nx}, it means that the rows of {X} span a proper linear
    subspace of {R^{nx}}. In that case, the best-fit map is not unique.
   
    This method uses internal working storage proportional {nx*(nx+nf)},
    like {lsq_fit}.
  */

/* AUXILIARY PROCEDURES */

void lsq_array_compute_matrix(uint32_t nt, uint32_t nx, double X[], double W[], double A[]);
  /* Assumes that {A} is an {nx} by {nx} matrix linearized by rows.
    Fills it with the linear system's matrix for the weighted least squares 
    problem with {nx} argument coordinates, {nt} data points, argument
    coordinates matrix {X} ({nt} by {nx} elements), and weights vector {W} ({nt} elements).
    See {lsq_array_fit}. for details.
    
    Namely, sets {A} tp {X'*diag(W)*X} where {X} denotes matrix trasposition; that is,
    sets {A[i*nx+j]} to {SUM{ W[k]*X[k*nx+i]*X[k*nx+j] : k \in 0..nt-1 }}
    for all {i,j} in {0..nx-1}. */

void lsq_array_compute_rhs(uint32_t nt, uint32_t nx, uint32_t nf, double X[], double F[], double W[], double B[]);
  /* Assumes that {B} is an array with {nx} rows and {nf} columns.
    Fills it with the right-hand-side of the linear system for the weighted least squares 
    problem with {nx} argument coordinates, {nt} data points, argument
    coordinates matrix {X} ({nt} by {nx} elements), function samples matrix {F} ({nt} by {nf} elements),
    and weights vector {W} ({nt} elements).  See {lsq_array_fit}. for details.
    
    Namely, sets {B} to {X'*diag(W)*F}; that is, set {B[i*nf+j]} to
    {SUM{ W[k]*X[k*nx+i]*F[k*nf+j] : k \in 0..nt-1 }} for all {i} in {0..nx-1}
    and all {j} in {0..nf-1}. */

#endif
