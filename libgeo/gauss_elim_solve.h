/* gauss_elim_solve.h - Solving linear eqs by Gaussian elimination. */
/* Last edited on 2024-11-25 04:13:08 by stolfi */

#ifndef gauss_elim_solve_H
#define gauss_elim_solve_H

#include <stdint.h>
   
uint32_t gauss_elim_solve(uint32_t m, uint32_t n, double A[], uint32_t p, double B[], double X[], double tiny);
  /* Solves the linear system {A X = B}, where {A} is a known matrix
    of size {m × n}, {B} is a known matrix of size {m × p}, and {X} is
    an unknown matrix of size {n × p}. The arrays {A} and {B} are not
    modified.
    
    Returns the number {r} of equations used in the computation of
    {X}, in the range {0..min(m,n)}. If {r < m}, the matrix {A} has linearly
    dependent rows, and the result {X} is meaningless. If {r == m} but
    {m < n}, then the system has multiple solutions, and {n-m} rows of
    {X} have been set to zero.
    
    During triangularization, any entry whose absolute value gets
    reduced to {tiny} or less times its previous value is set to zero.
    If {tiny} is zero or negative, this cleanup is supressed. */

void gauss_elim_solve_packed
  ( uint32_t m,
    uint32_t n,
    uint32_t p,
    double AB[],
    double X[],
    double tiny,
    uint32_t *rank_P,
    double *det_P
  );
  /* Like {gauss_elim_solve}, except that 
  
      (1) the arrrays {A} and {B} are taken to be the first {n} and the
      last {p} columns of the {m×(n+p)} array {AB[0..m*(n+p)-1]}.
    
      (2) the contents of the array {AB} is lost (modified by the Gauss
      elimination method).
      
      (3) if {rank_P} is not {NULL}, the number {r} of equations used in the computation of
      {X} is returned in {*rank_P};
      
      (4) if {det_P} is not {NULL}, the determinant of the original matrix {A}
      is returned in {*det_P}.  It will be zero if {m != n}. */

uint32_t gauss_elim_extract_solution
  ( uint32_t m,
    uint32_t n,
    double M[],
    uint32_t p,
    double X[]
  );
  /* Assumes that {M} is an {m×n} matrix that has been
    triangularized with {total = TRUE}, diagonalized, and normalized.
    
    The procedure interprets {M} as two blocks {M = (A B)} side by
    side, where {A} is {m×q} and {B} is {m×p}, with {q = n-p}. The
    procedure tries to compute a {q × p} matrix {X} such that {A X =
    B}. That is, it tries to solve the {p} linear systems {A x = u}
    where {x} is an unknown vector of size {q} and {u} is replaced in
    turn by each of the columns of {B}.
    
    Returns the number {r} of equations used in the computation of
    {X}, in the range {0..m}. If {r < m}, the matrix {Y = A X - B} may
    have {m-r} nonzero rows. If {r < q}, then the system has multiple
    solutions, and {q-m} rows of {X} have been set to zero.
    
    Specifically,
    
      * The returned result {r} is the number of rows {i < min(m,q)} such that
      {Lead(M,i) < q}.
    
      * If {j} in {0..q-1} is not {Lead(M,i)} for any
      {i}, then {X[j,0..p-1]} will be set to 0.
      
      * For every {i} in {0..m-1} such that {M[i,0..q-1] == 0},
      {Y[i,0..p-1]} will be equal to {-B[i,0..p-1]}.
    
    The matrix {X[i,k]} is stored in the given vector {X},
    linearized by rows. */



#endif
