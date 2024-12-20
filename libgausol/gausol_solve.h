/* gausol_solve.h - Solving linear eqs by Gaussian elimination. */
/* Last edited on 2024-11-30 23:34:04 by stolfi */

#ifndef gausol_solve_H
#define gausol_solve_H

#include <stdint.h>

#include <bool.h>

/* NOMENCLATURE

  The procedures in this and related modules are concerned with solving
  the system of "linear" (more properly, affine) equations {A X = B},
  where {A}, {X}, and {B} are matrices with sizes respectively {m×n},
  {n×p}, and {m×p}.
    
  Many procedures in this library will take or generate two tables
  {prow[0..m-1]}, a permutation of {0..m-1}, and {pcol[0..n-1]}, a
  permutation of {0..n-1}.  These tables define the /permuted versions/
  {P}, {Y}, and {Q} of {A}, {X} and {B}, respectively, where
        
        {P[i,j] = A[prow[i],pcol[j]} 
        {Y[j,k] = X[pcol[j],k]}
        {Q[i,k] = B[prow[i],k]}
        
    for {i} in {0..m-1}, {j} in {0..n-1}, and {k} in {0..p-1}.
    
    Note that the matrix product {P Y} is a permuted version of {A X},
    specificaly {(P Y)[i,k] = (A X)[prow[i],k]} for all relevant {i,j}.
    Thus the system {P Y = Q} is perfectly equivalent to {A X = B}
    
    Furthermore, many procedures in this module determine or accept an
    integer {rank} in {0..min(m,n)}, that defines a partition of the
    permuted matrices into blocks. The permutation {P} of {A} is split
    into into a {2×2} arrangement of 4 submatrices {P00,P01,P10,P11}
    respectively of size {rank×rank}, {rank×(n-rank)}, {(m-rank}×rank},
    and {(m-rank)×(n-rank)}. The permutation {Y} of {X} is split by rows
    into {Y0} of size {rank×p} and {Y1} of size {(n-rank)×p}. And the
    permutation {Q} of {B} is partitioned by rows into {Q0} of size
    {rank×p} and {Q1} {(m-rank)×p}.
    
    Note that {P00,P01,P10,Y0,Q0} are empty if {rank=0}; that {P10,P11,Q1} are
    empty if {rank=m}; and that {P01,P11,Y1} are empty if {rank=n}. All four
    submatrices of {P} are non-empty if and only if {rank} is in
    the range {1..min(m,n)-1}.
    
    The equation system {P Y = Q} then can be split into two 
    systems, {P00 Y0 + P01 Y1 = Q0} and {P10 Y0 + P11 Y1 = Q1}.
    
    All matrix arguments are assumed to be stored by rows into
    one-dimensional vectors of the proper size.  Thus the element
    {M[i,j]} of an {u×v} matrix argument {M} is actually {M[i*u + v]}.
    
    If the {prow} table argument is {NULL}, it is generally assumed that
    it is the identity permutation, {prow[i] = i} for all {i}.
    Likewise, if {pcol} table is {NULL}, it is assumed that it is the
    identity permutation, {pcol[j] = j} for all {j}. */
   
void gausol_solve
  ( uint32_t m,
    uint32_t n, double A[],
    uint32_t p, double B[], double X[],
    bool_t pivot_rows, bool_t pivot_cols,
    double tiny,
    double *det_P,
    uint32_t *rank_P
  );
  /* Uses the Gaussian elimination method to solve the system of
    "linear" (affine) equations {A X = B}, where {A} is a known matrix
    of size {m × n}, {B} is a known matrix of size {m × p}, and {X} is
    an unknown matrix of size {n × p}. The contents of the arrays {A}
    and {B} are not modified.
    
    If {pivot_rows} and {pvot_cols} are both true, uses full pivoting at
    each step. If only one of them is true, uses simple pivoting of rows
    or columns, respectively. All three options produce equivalent
    results for most purposes; the full-pivoting option is slower but
    may be less affected by roundoff errors.
    
    If both {pivot_rows} and {pvot_cols} are false, the procedure does
    no rearrangement of rows or columns at all. The pivot at each step
    will be the current disgonal element, if nonzero. This may fail
    unless the matrix {A} is strongly dominated by its diagonal.
    
    If {rank_P} is not {NULL}, the procedure returns in {*rank_P} the
    number {rank} of equations used in the computation of {X}, in the
    range {0..min(m,n)}. 
    
    If at least one of {pivot_rows} and {pivot_cols} is true, or {rank =
    min(m,n)}, or {rank = m-1 = n-1}, this will be the rank of the
    matrix {A}, apart from roundoff accidents. Otherwise the
    actual rank of {A} may anything in {rank..min(m,n)}.
    
    In any case, if {rank < m}, only {rank} of the equations (or
    combinations thereof) were used to compute {X}, and {m-rank}
    redundant or incompatible ones were ignored; so {X} probably does
    not satisfy the system. If {rank < n} then the procedure
    did not find enough independent equations to determine {X} uniquely,
    so {n-rank} rows of {X} have been set arbitrarily to zero.
      
    If {*det_P} is not {NULL}, the procedure will set {*detA} to its
    estimate of the determinant of the matrix {A}.  However, if the 
    matrix is not square, or {rank} may not be the actual rank of 
    {A}, as noted above, the returned value will be {NAN}.

    During the Gauss elimination method, the procedure will set to zero
    any entry in the internal matrices whose absolute value gets reduced to
    {tiny} or less. If {tiny} is zero or negative, this cleanup is
    supressed. */
   
void gausol_solve_in_place
  ( uint32_t m,
    uint32_t n, double A[],
    uint32_t p, double B[], double X[],
    bool_t pivot_rows, bool_t pivot_cols,
    double tiny,
    double *det_P,
    uint32_t *rank_P
  ); 
  /* Like {gausol_solve}, but the contents of the arrays 
    {A} and {B} are destroyed in the process. It 
    avoids allocating internal copies of the same. */
    

#endif
