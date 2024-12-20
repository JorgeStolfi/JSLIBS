/* Last edited on 2024-12-05 21:44:47 by stolfi */
/* test of sym_eigen_triang.h */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>

#include <sym_eigen_test_tools.h>

#include <sym_eigen_tridiag.h>

void tsetr_do_tests(FILE *wr1, FILE *wr2, FILE *wr3, uint32_t n, uint32_t tnum);
  /* Tests the three cases of {sym_eigen_tridiag} on an {n x n} array of
    type {tnum}. Writes the results to {wr1}, {wr2}, and {wr3},
    respectively. */

int32_t main(int32_t argn, char **argc);

#define NT (sym_eigen_test_tools_fill_matrix_NUM_TYPES) 
  /* Number of test matrix types. */

int32_t main(int32_t argn, char **argc)
  {
   
    FILE *wr1 = open_write("out/case1.txt", TRUE);
    FILE *wr2 = open_write("out/case2.txt", TRUE);
    FILE *wr3 = open_write("out/case3.txt", TRUE);
    
    for (uint32_t tnum = 0;  tnum < NT; tnum++)
      { tsetr_do_tests(wr1, wr2, wr3, 1, tnum);
        tsetr_do_tests(wr1, wr2, wr3, 2, tnum);
        tsetr_do_tests(wr1, wr2, wr3, 7, tnum);
      }
      
    fclose(wr1);
    fclose(wr2);
    fclose(wr3);

    return 0;
  }

void tsetr_do_tests(FILE *wr1, FILE *wr2, FILE *wr3, uint32_t n, uint32_t tnum)
  { 
    double A[n*n], R[n*n], d[n], e[n], v[n];
    double tiny = 0.5e-11;
    
    auto void ttr1(FILE *wr);
      /* Tests {sym_eigen_tridiag} with {R=NULL} */
      
    auto void ttr2(FILE *wr);
      /* Tests {sym_eigen_tridiag} with independent {R} */
      
    auto void ttr3(FILE *wr);
      /* Tests {sym_eigen_tridiag} with {R=A} */

    ttr1(wr1);
    ttr2(wr2);
    ttr3(wr3);

    return;
       
    void ttr1(FILE *wr)
      { 
        fprintf(wr, "\n");
        fprintf(wr, "=== C OUT - TRIDIAG (NO R) N=%d TYPE=%d ===\n", n, tnum);
        fprintf(wr, "\n");
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        fprintf(wr, "C - input matrix {A}\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,A, tiny);
        
        sym_eigen_tridiag(n,A,d,e,NULL);
        fprintf(wr, "C - computed tridiagonal form {d,e}\n");
        sym_eigen_test_tools_print_tridiag(wr, "%14.10f", n,d,e, tiny);
      }
       
    void ttr2(FILE *wr)
      { 
        fprintf(wr, "\n");
        fprintf(wr, "=== C OUT - TRIDIAG (NEW R) N=%d TYPE=%d ===\n", n, tnum);
        fprintf(wr, "\n");
        
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        fprintf(wr, "C - input matrix {A}\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,A, tiny);
        
        sym_eigen_tridiag(n,A,d,e,R);
        fprintf(wr, "C - computed tridiagonal form {d,e}\n");
        sym_eigen_test_tools_print_tridiag(wr, "%14.10f", n,d,e, tiny);
        fprintf(wr, "C - transformation matrix {R} (should be orthogonal)\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,R, tiny);
        sym_eigen_test_tools_check_orthogonal(n, n,n,R);

        /* Check whether the matrix {R} converts {A} to tridiagonal form: */
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        sym_eigen_test_tools_transform(n, n,n,R,A, v);
        fprintf(wr, "C - transformed matrix {R*A*R^t} (should be tridiagonal)\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,A, tiny);
        sym_eigen_test_tools_check_tridiagonal(n, n,A);
      }
       
    void ttr3(FILE *wr)
      { 
        fprintf(wr, "\n");
        fprintf(wr, "=== C OUT - TRIDIAG (R = A) N=%d TYPE=%d ===\n", n, tnum);
        fprintf(wr, "\n");
        
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        fprintf(wr, "C - input matrix {A}\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,A, tiny);
        
        /* Copy {A} to {R}: */
        for (int32_t ij = 0; ij < n*n; ij++) { R[ij] = A[ij]; }
        
        /* Destructive tridiagonalizarion of {R}: */
        sym_eigen_tridiag(n,R,d,e,R);
        fprintf(wr, "C - computed tridiagonal form {d,e}\n");
        sym_eigen_test_tools_print_tridiag(wr, "%14.10f", n,d,e, tiny);
        fprintf(wr, "C - transformation matrix {R} (should be orthogonal)\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,R, tiny);
        sym_eigen_test_tools_check_orthogonal(n, n,n,R);

        /* Check whether matrix {R} converts {A} to tridiagonal form: */
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        sym_eigen_test_tools_transform(n, n,n,R,A, v);
        fprintf(wr, "C - transformed matrix {R*A*R^t} (should be tridiagonal)\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,A, tiny);
        sym_eigen_test_tools_check_tridiagonal(n, n,A);
      }
  }
