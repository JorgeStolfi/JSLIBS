/* pst_interpolate.h -- procedures for interpolating weighted images. */
#ifndef pst_interpolate_H
#define pst_interpolate_H

/* Created on 2011-06-17 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-03-16 00:01:07 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <float_image.h>
 
void pst_interpolate_values(uint32_t n, double vS[], double wS[], uint32_t m, double *vR_P, double *wR_P);
  /* Given a list {vS[0..n-1]} of values of a function {f} at {n}
    consecutive equidistant arguments {x} in {0..n-1}, and an integer
    {m} in {0..n}, the procedure estimates the value {vR = f(m - 0.5)}.
    It uses extrapolation (if {m} is 0 or {n}) or intrepolation of
    {vS[0..n-1]} with a polynomial of degree {n-1} (constant, affine,
    quadratic, or cubic).  
    
    The procedure also assumes that {wS[k]}, for {k} in {0..n-1}, is a
    positive reliability weight for {vS[k]}, proportional to the
    reciprocal of the variance of the noise in {vS[k]}. The procedure
    computes a reliability weight {wR} for the result {vR}, from those
    weights and from the extrapolation/interpolation formula used. It
    returns the results in {*vR_P} and {*wR_P}.
    
    In particular, if {n} is zero, the result {vR} will be {NAN} and its
    weight {wR} will be zero. Idem if overflow or underflow occurs.
    
    All the values {vS[0..n-1]} must be finite, and all the weights {wS[0..n-1]}
    must be finite and positive. */
    
typedef void pst_iterpolate_get_data_func_t (int32_t j, double *vR_P, double *wR_P);
  /* Type of a procedure that fetches a data value {vR} and its reliability weight {wR}
    given an integer index {j}, possibly negative, and returns them in {*vr_P}
    and {*wR_P}.  The returned results may be {vR=NAN} and {wR=0}.  Otherwise
    {vR} must be finite and {wR} must be finite and positive. 
    The procedure must eventually return {vR=NAN} and {wR=0} if {j} is too large 
    or too small. */
    
void pst_interpolate_select_data
  ( int32_t j0,
    pst_iterpolate_get_data_func_t func,
    int32_t *ja_P,
    int32_t *jb_P,
    int32_t *n_P,
    int32_t *m_P
  );
  /* Evaluates {func(j)} trying to obtain up to four data pairs with
    positive weights, at consecutive {j} indices {ja..jb} that include
    {j0} and/or {j1=j0+1} and are as balanced as possible around the
    midpoint {ctr=j0+0.5} of the two. Returns the indices {ja,jb} in
    {*ja_P} and {*jb_P}.
    
    Also returns in {*n_P} the number of pairs {jb-ja+1}
    and in {*m_P} the number {m=j1-ja} of indices not greater than {j0}.
    
    This procedure is useful to obtain data for {pst_interpolate_values}. */

#endif
