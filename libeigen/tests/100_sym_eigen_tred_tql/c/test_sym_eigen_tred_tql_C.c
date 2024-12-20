/* Last edited on 2024-12-05 18:31:04 by stolfi */
/* test of sym_eigen_tred1.h, sym_eigen_tred2.h, sym_eigen_tql1.h, sym_eigen_tql2.h,  */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>

#include <sym_eigen_test_tools.h>

#include <sym_eigen_tred1.h>
#include <sym_eigen_tred2.h>
#include <sym_eigen_tql1.h>
#include <sym_eigen_tql2.h>

void tsett_do_tests(FILE *wr1, FILE *wr2, uint32_t n, uint32_t tnum);
  /* Tests the two versions of the tridiagonalization and QL 
    procedures on an {n x n} array of type {tnum}.  Writes
    the results to {wr1} and {wr2}, respectively. */
    
int32_t main(int32_t argn, char **argc);

#define NT (sym_eigen_test_tools_fill_matrix_NUM_TYPES) 
  /* Number of test matrix types. */

int32_t main(int32_t argn, char **argc)
  {
  
    FILE *wr1 = open_write("out/version1.txt", TRUE);
    FILE *wr2 = open_write("out/version2.txt", TRUE);
    
    for (uint32_t tnum = 0;  tnum < NT; tnum++)
      { tsett_do_tests(wr1, wr2, 1, tnum);
        tsett_do_tests(wr1, wr2, 2, tnum);
        tsett_do_tests(wr1, wr2, 7, tnum);
      }
      
    fclose(wr1);
    fclose(wr2);

    return 0;
  }
    
void tsett_do_tests(FILE *wr1, FILE *wr2, uint32_t n, uint32_t tnum)
  { 
    double A[n*n], R[n*n], d[n], e[n], v[n];
    double tiny = 0.5e-11;

    auto void ttr1(FILE *wr);
      /* Tests {sym_eigen_tred1}, {sym_eigen_tql1}. */
    
    auto void ttr2(FILE *wr);
      /* Tests {sym_eigen_tred2}, {sym_eigen_tql2}. */
    
    ttr1(wr1);
    ttr2(wr2);
    
    return;
       
    void ttr1(FILE *wr)
      { 
        fprintf(wr, "\n");
        fprintf(wr, "=== C OUT - TRED1/TQL1 (NO Z) N=%d TYPE=%d ===\n", n, tnum);
        fprintf(wr, "\n");
        
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        fprintf(wr, "C OUT - input matrix\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,A, tiny);
        
        sym_eigen_tred1(n,A,d,e,R);
        fprintf(wr, "C OUT - tridiagonal form\n");
        sym_eigen_test_tools_print_tridiag(wr, "%14.10f", n,d,e, tiny);

        /* Compute eigenvalues: */
        uint32_t nev;
        sym_eigen_tql1(n,d,e,&nev);
        fprintf(stderr, "  n = %d nev = %d\n", n, nev);
        assert(nev <= n);
        fprintf(wr, "C OUT - eigenvalues\n");
        sym_eigen_test_tools_print_eigenvalues(wr, "%14.10f", nev,d);
      }
       
    void ttr2(FILE *wr)
      { 
        fprintf(wr, "\n");
        fprintf(wr, "=== C OUT - TRED2/TQL2 (Z.NE.A) N=%d TYPE=%d ===\n", n, tnum);
        fprintf(wr, "\n");
        
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        fprintf(wr, "C OUT - input matrix\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,A, tiny);
        
        sym_eigen_tred2(n,A,d,e,R);
        fprintf(wr, "C OUT - tridiagonal form\n");
        sym_eigen_test_tools_print_tridiag(wr, "%14.10f", n,d,e, tiny);
        fprintf(wr, "C OUT - orthogonal transformation\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,R, tiny);

        /* Check whether matrix {R} makes {A} tridiagonal: */
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        sym_eigen_test_tools_transform(n, n,n,R,A, v);
        fprintf(wr, "C OUT - transformed matrix (tridiagonal)\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, n,n,A, tiny);

        /* Compute eigenvalues and eigenvectors: */
        uint32_t nev;
        sym_eigen_tql2(n,d,e,R,&nev);
        fprintf(stderr, "  n = %d nev = %d\n", n, nev);
        assert(nev <= n);
        fprintf(wr, "C OUT - eigenvalues\n");
        sym_eigen_test_tools_print_eigenvalues(wr, "%14.10f", nev,d);
        fprintf(wr, "C OUT - eigenvectors\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, nev,n,R, tiny);

        /* Check whether matrix {R} turns {A} diagonal: */
        sym_eigen_test_tools_fill_matrix(n, n,A, tnum);
        sym_eigen_test_tools_transform(n, nev,n,R,A,v);
        fprintf(wr, "C OUT - transformed matrix (diagonal)\n");
        sym_eigen_test_tools_print_matrix(wr, "%14.10f", n, nev,nev,A, tiny);
        fprintf(wr, "\n");
      }
  }
