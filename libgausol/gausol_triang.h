/* gausol_triang_reduce.h - Gaussian triangulation with full pivoting. */
/* Last edited on 2024-11-29 22:55:24 by stolfi */

#ifndef gausol_triang_reduce_H
#define gausol_triang_reduce_H

#include <stdint.h>

#include <bool.h>

void gausol_triang_reduce
  ( uint32_t m, uint32_t prow[],
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    double tiny,
    uint32_t *rank_P,
    double *det_P
  );
  /* Modifies the {m×n} matrix {A} and the {m×p} matrix {B} of a linear
    system {A X = B} according the Gaussian elimination method, leaving
    {A} upper triangular modulo a permutation of the rows and columns.
    
    Specifically, the matrices {A} abd {B} should be stored by rows.
    That is, {A} and {B} should have {m*n} and {m*p} elements,
    respectively, {A[i,j]} should be actually {A[i*n + j]}, and {B[i,k]}
    should be {B[i*p + k]}.

    Effect of the procedure:
    
    The procedure modifies {A} and {B} by a sequence of row operations
    {A <- R A} and {B <- R B} where {R} is an {m×m} matrix with unit
    determinant.
    
    During these row operations, any element whose absolute value was or
    becomes smaller than {tiny} in absolute value gets set to zero. If
    {tiny} is zero or negative, this cleanup is supressed.
    
    The procedure also returns in {*rank_P} an integer {rank} in {0..min(m,n)} that
    the max number of linearly independent rows (or columns) of {A}.
    
    The procedure also fills tables {prow[0..m-1]} and {pcol[0..n-1]}
    with permutations of the indices {0..m-1} and {0..n-1},
    respectively.
    
    Also, if {*det_P} is not {NULL}, the matrix {A} is square, and
    {rank} is {m}, the procedure stores into {*det_P} the determinant of
    {A}, computed as the product of the diagonal elements of the
    triangulated and permuted matrix {P = P00}, with sign reversed if
    the two permutations imply an odd total number of pair swaps. If the
    matrix is square, {rank} is less than {m}, and at least one of
    {prow} or {pcol} is not {NULL} (meaning that at least partial
    pivoting was used), sets {*det_P} to zero. In all other cases, sets
    {*detP} to {NAN}.
    
    The tables {prow} and {pcol} define permuted versions {P} of {A}
    and {Q} of {B}, as explained in {gausol_solve.h}.    
    On output of the procedure,
    
      (1) {P[i,i]} is nonzero for {i} in {0..rank-1}.
      
      (2) {P[k,j]} is zero for {j} in {0..rank-1} and {k} in {j+1..m-1}.
      
    These conditions say that {P00} is upper triangular with nonzero
    diagonal, and that {P10} is all zeros. Moreover, condition (3) must
    hold, which is, either {rank == min(m,n)}, in which case {P11} is
    empty, or one of the following is
    
      (3.00) The first element {P11[0,0] = P[rank,rank]} is zero.

      (3.10) If {prow} is not {NULL}, the first column of {P11} (namely
      {P[i,rank]} for {i} in {rank..m-1}) is zero.
    
      (3.01) If {pcol} is not {NULL}, the first row of {P11} (namely
      {P[rank,j]} for {j} in {rank..n-1}) is zero..
      
      (3.11) if both {pcol} and {prow} are non-{NULL}, then the whole of {P11} 
      (namely {P[i,j]} for {i} in {rank..m-1} and {j} in {rank..n-1}) is zero.
    
    Note that (3.10) and (3.01) imply (3.00), and (3.11) implies all
    three. In cases (3.10), (3.01), and (3.11), the rank of {A} was
    {rank} (apart from roundoff errors). If only (3.00) holds (that is,
    {prow} and {pcol} were both {NULL}), the rank of {A} was at least
    {rank}.
    
    Uses for {prow=pcol=NULL}:
    
    This combination can be used when it is known that the coefficient
    matrix {A} has a strongly dominating diagonal, so that pivoting will
    not be needed, not even the simple version. It could be used, for
    example, to refine the inverse of a matrix {A} after computing a
    first approximation {M0}. For that one may compute the inverse {D1} of
    {C1 = A M0} and set the new inverse as {M1 = M0 D1}. The matrix {C1}
    should be very close to the identity so no pivoting should be
    necessary.
    
    Uses when exactly one of {prow} and {pcol} is {NULL}:
    
    A typical application of this case is finding the rank of a
    matrix {A}. This is the same as that of {P}, hence it is the
    returned result {rank}.
    
    Another typical application for this case is to compute the
    determinant of a square {m×m} matrix {A}. If the returned {rank} is
    less than {m} then the determinant is zero, otherwise it is {sgn}
    times the determinant of {P00}, and the latter is the product of
    {P[i,i]} for {i} in {0..rank-1}.
    
    Another typical application is solving the {p} linear equation systems 
    {A X = B} where {X} is {n×p} {p}. This equation is equivalent to 
    {P Y = Q} where {Y[j,k] = X[pcol[j],k]} for {j} in {0..n-1} and {k} in
    {0..p-1}. This in turn is equivalent to {P00 Y0 + P01 Y1 = Q0} and 
    {P11 Y1 = Q1} where {Y0} and {Y1} are the first
    {rank} and last {n-rank} rows of {Y}, respectively.  But
    if {rank<min(m,n)} then {P11} has at least one row and one column
    and the first column is all zeros; meaning that the whole equation
    {P Y = Q} has either no solution (if {Q1} is nonzero) 
    or infintely many solutions (if {Q1} is all zeros).
    
    A special case of this application is when {m=n=p} and {B} is the identity,
    in which case the required solution {X} is the inverse of {A}.  
    If {rank<m} (and so {rank<n}) the input matrix {A} is singular, and the 
    strict inverse does not exist. 
    
    Uses with both {prow} and {pcol} non-null:
    
    This use case has the same applciations as the previous one,
    but it may be more resistant to roundoff errors, although it takes 
    slightly more time.
    
    Another application of this use case, specifically, is finding SOME
    solution to a system of {p} linear equation systems {A X = B} which
    may have many solutions. If {rank<min(m,n)} but the system has a solution, then,
    {Q1} shoul be all zeros too, so that the equation {P11 Y1 = Q1} is 
    trivially satisfied by ANY {Y1}.  Then we can set {Y1} to all zeros
    (or any other matrix} and use the equation {P00 Y0 = Q0 - P01 Y1}
    to compute {Y0}.  Since {P00} is upper triangular, this system
    can be solved one row at a time, from {rank-1} down to {0}.
    
    Alternatively, if {Q1} is not zero, the procedure can be used to find a 
    matrix {X} that satisfies the equations {pcol[0..rank-1]}, ignoring the
    other equations because they are either redundant or impossible
    given those {rank} are satisfied.
    
    Another application is finding a basis for the linear subspace
    spanned by the rows of a matrix {A}. For this purpose, one would
    call this procedure with {p=0}. Upon return the dimension of the row
    space will be {rank}, and the desired basis will be in rows
    {prow[0..rank-1]} of {A}.
    
    Method used by this procedure:
    
    The procedure starts with {rank=0}, {pcol[j] = j} for {j} in {0..n-1},
    and {prow[i] = i} for {i} in {m-1}. At this point conditions (1) and
    (2) are trivially true. It then tries to increase {rank} by one at each
    iteration, while maintaining conditions (1) and (2). When it fails to do so,
    one of conditions (3.00), (3.10), (3.01), or (3.11) will be true. 
    
    At each iteration, the procedure looks for a nonzero /pivot/ element
    {P[i,j] = A[prow[i],pcol[j]]} in some subset of the subarray {P11},
    that is, with {i} in {rank..m-1} and {j} in {rank..n-1}. If it finds
    such pivot, it adds a multiple of row {i} of {P} and {Q} to every
    other row {k} in {rank..m-1} except {i}, with the multiplier chosen
    so as to make {P[k,j]} zero. Then it modifies {pcol[rank..n-1]} and
    {prow[rank..m-1]} so that the pivot becomes {P[rank,rank]}. Then it
    increments {rank}. The procedure stops when it fails to find a pivot
    -- either because {rank} became equal to {m} or {n} (so that {P11}
    is empty), or because the searched part of {P11} is all zeros.
    
    When {prow=pcol=NULL}, the procedure tries to use {P11[0,0]=P[rank,rank]} as the pivot.  
    If that element is zero, condition (3.00) is satisfied.
    
    When {prow!=NULL} but {pcol=NULL}, it looks for the pivot only among
    the first column {P11}, that is, column {rank} of {P}, rows
    {rank..m-1}. If all those elements are zero, condition (3.10) holds.
    
    Similarly, when {prow=NULL} but {pcol!=NULL}, the procedure looks
    for the pivot only among the first row {P11}, that is, row {rank} of
    {P}, columns {rank..n-1}. If all those elements are zero, condition
    (3.01) holds.
    
    When both {prow} and {pcol} are non-{NULL}, the procedure looks for
    a pivot in the whole submatrix {P11}.  If that whole matrix is zero,
    condition (3.11) is satisfied.
    
    In any case, if there are two or more nozero candidates for the
    pivot, the procedure chooses the one which has largest absolute
    value. */

void gausol_triang_diagonalize
  ( uint32_t m, uint32_t prow[],
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    uint32_t rank,
    double tiny
  );
  /* To be called after {gausol_triang_reduce}, to turn the 
    main submatrix of {A} (permuted) from upper triangular to diagonal,
    updating {B} accordingly.
  
    Specifically, assumes that the {m×n} matrix {A} and the {m×p} matrix
    {B} have been processed by {gausol_triang_reduce} with parameters
    {m,n,p,tiny} resulting in apparent rank {rank} and row and column
    permutation tables {prow[0..m-1]} (if not {NULL}) and {pcol[0..n-1]}
    (if not {NULL}).
    
    See {gausol_solve,h} for the definition of the permuted versions
    {P} and {Q} of {A} abd {B}, and their submatrices {P00,P01,Q0,P10,P11,Q1}.
    
    This procedure assumes that {P00} (the first {rank} rows and colums of
    {P}) is upper triangular, that is, {P[i,j]} is zero for {i} in
    {0..rank-1} and {j} in {0..i-1}; and that the diagonal elements are
    nonzero. It applies row operations to {A} and {B}, such that the submatrix
    {P00} becomes diagonal, without changing the diagonal elements.  Only the
    submatrices {P01} and {Q0} and the part of {P00} above the main
    diagonal are affected. */

void gausol_triang_normalize
  ( uint32_t m, uint32_t prow[],
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    uint32_t rank,
    double tiny
  );
  /* To be called after {gausol_triang_reduce} and
    {gausol_triang_diagonalize}, to turn the main submatrix of {A}
    (permuted) from diagonal to the identity, while updating {B}
    accordingly.
  
    Specifically, assumes that the {m×n} matrix {A} and the {m×p} matrx
    {B} have been processed by {gausol_triang_reduce}, resulting in
    apparent rank {rank} and row and column permutation tables
    {prow[0..m-1]} (if not {NULL}) and {pcol[0..n-1]} (if not {NULL});
    and then by {gausol_triang_diagonalize} with parameters
    {m,n,p,tiny}.
    
    See {gausol_solve,h} for the definition of the permuted versions
    {P} and {Q} of {A} abd {B}, and their submatrices {P00,P01,Q0,P10,P11,Q1}.
    
    This procedure assumes that {P00} is diagonal with nonzero diagonal
    elements. It divides each row {i} of {P01,Q0} by the diagonal
    element {P[i,i]}, and sets {P[i,i]=1}, for {i} in {0..rank-1}. The rest
    of {P} is not affected. */

#endif
