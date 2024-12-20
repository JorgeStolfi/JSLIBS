/* Tools for testing the {libeigen} routines. */
/* Last edited on 2024-12-05 17:33:45 by stolfi */

#ifndef sym_eigen_test_tools_H
#define sym_eigen_test_tools_H

#include <stdio.h>
#include <stdint.h>

/* In the following procedures, a matrix {M} is assumed to be stored 
  into a vector by rows of size {rowsz}, so that element {M[i,j]}
  is actually {M[i*rowsz + j]}. */
  
#define sym_eigen_test_tools_fill_matrix_NUM_TYPES 5  

void sym_eigen_test_tools_fill_matrix(uint32_t rowsz, uint32_t n, double A[], uint32_t tnum);
  /* Fills matrix {A[0..n-1,0..n-1]} with the 
    symmetric test matrix of type {tnum}. */

void sym_eigen_test_tools_print_matrix
  ( FILE *wr, char *fmt, uint32_t rowsz,
    uint32_t m, uint32_t n,  double A[],
    double tiny
  );
  /* Prints the (sub)matrix {A[0..m-1, 0..n-1]}, replacing 
    values less than {tiny} by blanks.  (Pass {tiny} negative to 
    print all elements, incluiding zeros.)  */

void sym_eigen_test_tools_print_tridiag(FILE *wr, char *fmt, uint32_t n, double d[], double e[], double tiny);
  /* Prints the symmetric tridiagonal matrix whose diagonal is
    {d[0..n-1]} and whose sub-diagonal is {e[1..n-1]}. Replaces values
    in the latter which less than {tiny} by blanks. (Pass {tiny}
    negative to print all elements, incluiding zeros.) */

void sym_eigen_test_tools_print_eigenvalues(FILE *wr, char *fmt, uint32_t nev, double d[]);
  /* Prints the eigenvalue list {d[1..nev-1]}. */

void sym_eigen_test_tools_transform(uint32_t rowsz, uint32_t m, uint32_t n, double R[], double A[], double v[]);     
  /* Conjugates the {n × n} matrix {A} by the orthogonal {m × n} matrix
    {R}, that is, computes the {m × m} matrix {B=R*A*(R^t)}, using
    {v[0..n-1]} as temp storage.
    The result is stored back into the first {m} rows and columns of
    {A}. */

void sym_eigen_test_tools_check_orthogonal(uint32_t rowsz, uint32_t m, uint32_t n, double R[]);
  /* Checks whether the {m × n} matrix {R} is orthogonal.
   That is, whether {S = R*R^t} is the {m × m} identity
   apart from roundoff. */
   
void sym_eigen_test_tools_check_tridiagonal(uint32_t rowsz, uint32_t n, double A[]);
  /* Checks whether the {n × n} matrix {A} is indeed
    tridiagonal; that is, whether elements {A[i,j]}
    with {|i-j|>1} are zero, apart from roundoff. */
   
void sym_eigen_test_tools_check_diagonal(uint32_t rowsz, uint32_t n, double A[], double d[]);
  /* Checks whether the {n × n} matrix {A} is indeed
    diagonal, with diagonal elements {d[0..n-1]}; that is, 
    whether {A[i,i]=d[i]} and {A[i,j]=0} for all {i}
    and all {j!=i}, apart from roundoff. */

#endif
