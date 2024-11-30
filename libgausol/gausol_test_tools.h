/* gausol_test_tools.h - basic tools for testing {gausol_XXX}.h. */
/* Last edited on 2024-11-30 04:53:17 by stolfi */

#ifndef gausol_test_tools_H
#define gausol_test_tools_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

void gausol_test_tools_choose_system
  ( uint32_t trial,
    uint32_t m_max, uint32_t n_max, uint32_t p_max,
    uint32_t *m_P, 
    uint32_t *n_P, double **A_P,
    uint32_t *p_P, double **B_P,
    double **X_P,
    double tiny,
    bool_t verbose
  );
  /* Choses the problem size {m,n,p}, and allocates the arrays
    {A[0..m-1,0..n-1]}, {B[0..m-1,0..p-1]}, and {X[0..n-1,0..p-1}. Then
    fills {A} and {X} with random numbers, then computes {B = A*X}. If
    {verbose}, also prints the system to {stderr}. Computed {B}
    entries whose absolute value is less than {tiny} are set to zero. */

void gausol_test_tools_check_triang_reduce
  ( uint32_t m, uint32_t prow[], 
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    double tiny,
    uint32_t rank, 
    double det_ref, double det_cmp,
    double rms_A, 
    bool_t verbose
  );
  /* Checks whether {A} and {B} have the shape expected after
    {gausol_triang_reduce} returned {prow}, {pcol}, {rank} and {sgn}.
    Namely, conditions (1), (2), and one of (3.00), (3.01), (3.10), (3.11) should be
    true. In fact, either (3a), (3b), (3c), or (3d) should be true,
    depending on whether {prow} and/or {pcol} are {NULL} or not
    
    Also checks whether {prow} and {pcol}, if not {NULL}, are permutations of {1..m}
    and {1..n}, respectively, and {sgn} is equal to the number of swaps
    performed by them. 
    
    Also checks whether the determinant {det_ref} of {A} computed by some 
    other means matches the determinant {det_cmp} returned by {gausol_triang_reduce}.
    This test takes place only if {det_ref} is not {NAN}, which should happen
    only if {A} is square.
    
    The procedure assumes that the RMS value of the entries of {A}, before
    triangulation, was {rms_A}, 

    Should check whether the linear span of the rows is preserved, but
    currently does not. */

void gausol_test_tools_check_diagonalize
  ( uint32_t m, uint32_t prow[], 
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    uint32_t rank,
    bool_t verbose
  );
  /* Checks whether {A} and {B} have the shape expected after
    {gausol_triang_diagonalize}, asuming that they passed
    {check_triang_reduce} with the given parameters.

    Namely, checks that {P[i,j]} is zero for 
    {i} in {0..rank-1} and {j} in {i+1..rank-1}. */

void gausol_test_tools_check_normalize
  ( uint32_t m, uint32_t prow[],
    uint32_t n, uint32_t pcol[], double A[],
    uint32_t p, double B[],
    uint32_t rank,
    bool_t verbose
  );
  /* Checks that {A} and {B} have the shape expected after
    {gausol_triang_normalize}, assuming that they passed 
    {check_triang_reduce} and {check_diagonalize}. */

void gausol_test_tools_check_satisfaction
  ( uint32_t m,  
    uint32_t n,  double A[],
    uint32_t p, double B[],
    double X_ref[],
    bool_t verbose
  );
  /* Expects {A} and {B} to be {m×n} and {m×p} matrces, stored by rows. 
    Checks whether the {n×p} matrix {X_ref} is still a
    solution of the system {A X = B}. */

double gausol_test_tools_check_solve
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p, double B[],
    double X_ref[], double X_cmp[],
    uint32_t rank,
    bool_t verbose
  );
  /* Assumes that {A}, {B} are the original {m×n} and {m×p} matrices of
    a system of equations {A X = B} (before {gausol_triang_reduce}), and
    {X_ref} is the {n×p} nominal solution of the system -- that is, {B}
    was computed as {A Xref}; that {X_cmp} is the {n×p} solution of {A X
    = B} computed by {gausol_solve}; and that {rank} is the rank of {A}
    as determined by {gausol_solve}.
    
    The rank {rank} must be in {0..min(m,n)}. If {m=n=rank}, assumes that the
    system {A X + B} has a unique solution, so it compare the solutions
    {Xref} and {X_cmp}.  IN this case, the RMS difference between the elements 
    of these arrays is returned as result.  If {rank != m} or {rank != n},
    returns {NAN}.
    
    If {rank=n} (for any {m}) assumes that the system {A X = B} has at
    least one solition, so it also compares the product {A X_cmp} with the
    RHS matrix {B}. */

double gausol_test_tools_compare_solutions
  ( uint32_t n, uint32_t p,
    double X_ref[], double X_cmp[],
    double rel_tol, double abs_tol,
    bool_t verbose
  );
  /* Checks a putative solution {X_cmp} of a system {A X = B} against
    the `true' solution {X_ref}, with relative tolerance {rel_tol} and
    absolute toleance {abs_tol}. Assumes that both {X_cmp} and {X_ref}
    are {n × p} matrices.  Returns the RMS difference between the elements
    of the two arrays. */

void gausol_test_tools_check_residual
  ( uint32_t m, uint32_t n, double A[],
    uint32_t p, double B[],
    double X[],
    double rel_tol, double abs_tol,
    bool_t verbose
  );
  /* Assumes that {A}, {B}, and {X} are matrices with sizes
    {m×n}, {m×p}, and {n×p}, respectively. Compares the product
    {A X} with {B}. */

void gausol_test_tools_compare_determinants
  ( uint32_t m, uint32_t n,
    uint32_t rank,
    double rms_A,
    double det_ref, double det_cmp,
    bool_t verbose
  ); 
  /* Assumes {det-ref} is the determinant of the original
    matrix {A} computed by some other means, and {det_cmp} 
    is the determinant returned by {gausol_triang_reduce}.
    Compares the two, and bombs out if they differ by too
    much. 
    
    Assumes that the RMS value of the entries of {A}, before
    triangulation, was {rms_A}, and that entries less than {tiny} in
    absolute value were replaced by zeros. 
    
    The check is skipped if {det_ref} and/or {det_cmp} are
    {NAN}.  This is the case when the matrix {A} is not square,
    or pivoting was suppressed and the apparent {rank} was less
    than {m}. */

void gausol_test_tools_make_row_dependent(uint32_t i, uint32_t m, uint32_t n, double A[]);
  /* The matrix {A} must have {m} rows and {n} columns, and {i} must be in {0..m-1}.
    Fills row {i} of {A} with a random linear combination of the other rows.
    In particular, if {m} is 1, then {i} must be 0, and the procedure
    fills that only row with zeros. */

void gausol_test_tools_multiply
  ( uint32_t m, uint32_t n, uint32_t p,
    double A[], double X[], double B[],
    double tiny
  );
  /* Multiplies the {m×n} matrix {A} by the {n×p} matrix {X} 
    and stores the result into the {m×p} matrix {B}. 
    Any computed element that is less than {tiny}
    is set to zero. */

#define gausol_test_tools_det_by_enum_SIZE_MAX 9
  /* Max value of {q} for {rmxn_det_by_enum}. Note that {9! = 362'880}. */
    
double gausol_test_tools_det_by_enum(uint32_t m, uint32_t n, double A[], uint32_t q);
  /* Determinant of the first {q} rows and columns of the {m×n} matrix {A}, computed by
    the elementary definition (sum of {q!} products of elements of {A})
    Returns zero if {q > m} or {q > n}.  Otherwise {q}
    must not exceed {gausol_test_tools_det_by_enum_SIZE_MAX}. */
     
double gausol_test_tools_elem_RMS(uint32_t m, uint32_t n, double A[]);
  /* The root mean square value of all elements of the {m×n} matrix {A}. */
 
double *gausol_test_tools_throw_matrix
  ( uint32_t m, uint32_t n,
    double pzero,
    double mag_min, double mag_max,
    double tiny,
    char *head,
    char *name,
    bool_t verbose
  );
  /* Generates a random {m × n} matrix {A}. 
  
    Some elements will be set to zero, with probability
    {pzero/min(m,n)}. Thus, the average number of zeros in each row or
    column (whichever is shorter) will be {pzero}.
    
    Each elements that is not set to zero by that lottery will be a
    random number uniform in {[-mag _ +mag]}, where {mag} is a random
    number with exponential distribution in {[mag_min _ mag_max]}.
    
    However, elements that would be less than {tiny} in absolute value
    are set to zero. 
    
    If {verbose} is true, prints the array to {stderr} with the 
    given {head} and {name}.  */

void gausol_test_tools_throw_system
  ( uint32_t m,
    uint32_t n, double **A_P,
    uint32_t p, double **X_P,
    bool_t verbose
  );
  /* Stores into {*A_P} a random {m × n} coefficient matrix {A}, and 
    into {*X_P} a random {m × p} solution matrix {X}. */
    
#endif

