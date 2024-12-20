
/* See {sym_eigen_test_tools.h} */
/* Last edited on 2024-12-05 11:23:12 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsprintf.h>

#include <sym_eigen_test_tools.h>
  
char *sym_eigen_test_tools_blank_zero(char *fmt);
  /* Returns a stack-allocated "blank" string to replace
    tiny {float} or {double} elements that would be printed 
    with {fmt};  The {fmt} is expected to generate a '.'
    somewhere, and the result will have a '.' in the same
    position. */

void sym_eigen_test_tools_fill_matrix(uint32_t rowsz, uint32_t n, double A[], uint32_t tnum)
  { 
    auto double elem_diagonal(int32_t i, int32_t j);
      /* Diagonal with {A[i,i]=i+1}. */
      
    auto double elem_triblock(int32_t i, int32_t j);
      /* Tridiagonal with {2 × 2} blocks. */
      
    auto double elem_randomic(int32_t i, int32_t j);
      /* Pseudorandom. */
      
    auto double elem_rankdeux(int32_t i, int32_t j);
      /* A rank 2 matrix. */

    auto double elem_hartleys(int32_t i, int32_t j);
      /* A matrix with Hartley basis eigenvectors 
        and distinct eigenvalues. */
 
    for (int32_t i = 0;  i < n; i++) 
      { for (int32_t j = 0;  j <= i; j++) 
          { double Aij;
            switch(tnum)
              { case 0: Aij = elem_diagonal(i,j); break;
                case 1: Aij = elem_triblock(i,j); break;
                case 2: Aij = elem_randomic(i,j); break;
                case 3: Aij = elem_rankdeux(i,j); break;
                case 4: Aij = elem_hartleys(i,j); break;
                default: demand(FALSE, "invalid test number");
              }
            A[i*(int32_t)rowsz + j] = Aij;
            A[j*(int32_t)rowsz + i] = Aij;
          }
      }
    return;
    
   double elem_diagonal(int32_t i, int32_t j)
      { return (i == j ? i+1 : 0); }
      
    double elem_triblock(int32_t i, int32_t j)
      { int32_t k = (i/2); /* Block index. */
        if (i == j)
          { return 1.5 + 2*k; }
        else if ((j+1 == i) && ((i+j) % 4 == 1))
          { return 0.5; }
        else
          { return 0.0; }
      }
      
    double elem_randomic(int32_t i, int32_t j)
      { double s = 1.0 + i + j - n;
        double t = (i+1.0)*(j+1.0);
        return 0.0001/((fabs(s)+1.0)*sqrt(t));
      }
      
    double elem_rankdeux(int32_t i, int32_t j)
      { /* Rank 2: */
        double ui0 = 1.0/(1+i);
        double ui1 = 1.0/(1+i*i);
        double uj0 = 1.0/(1+j);
        double uj1 = 1.0/(1+j*j);
        return ui0*uj0 + ui1*uj1;
       }
       
    double elem_hartleys(int32_t i, int32_t j)
      { double Aij = 0;
        for (int32_t r = 0; r < n; r++)
          { double Rir = cos((i*r*2.0/n + 0.25)*M_PI);
            double Rjr = cos((j*r*2.0/n + 0.25)*M_PI);
            double dr = r + 1;
            Aij += Rir*dr*Rjr/(2*n);
          }
        return Aij;
      }
  }

char *sym_eigen_test_tools_blank_zero(char *fmt)
  { 
   /* Figure out the char count {iwd} up to the '.' and {fwd} after it:*/ 
    char *xel = jsprintf(fmt, -M_PI);
    int32_t wd = (int32_t)strlen(xel);
    char *pt = strchr(xel, '.');
    int32_t iwd = (pt == NULL ? wd : (int32_t)(pt - xel + 1));
    int32_t fwd = wd - iwd;
    assert((iwd > 0) && (fwd >= 0));
    free(xel);
    
    /* Create the blank string to replace tiny elements: */
    char *xzero = jsprintf("%*s%*s", iwd, ".", fwd, "");
    return xzero;
  }

void sym_eigen_test_tools_print_matrix
  ( FILE *wr, char *fmt, uint32_t rowsz,
    uint32_t m, uint32_t n,  double A[],
    double tiny
  )
  { 
    char *xzero = sym_eigen_test_tools_blank_zero(fmt);
    /* Print array: */
    fprintf(wr, "\n");
    for (int32_t i = 0;  i < m; i++) 
      { for (int32_t j = 0;  j < n; j++) 
          { double Aij = A[i*(int32_t)rowsz + j];
            if (fabs(Aij) < tiny)
              { fputs(xzero, wr); }
            else
              { fprintf(wr, fmt, Aij); }
          }
        fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
    free(xzero);
  }

void sym_eigen_test_tools_print_tridiag(FILE *wr, char *fmt, uint32_t n, double d[], double e[], double tiny)
  { 
    char *xzero = sym_eigen_test_tools_blank_zero(fmt);
    fprintf(wr, "\n");
    for (int32_t i = 0;  i < n; i++) 
      { for (int32_t j = 0;  j < n; j++)
        if (i == j)
          { fprintf(wr, fmt, d[i]); }
        else  
          { double eij;
            if (i == j+1)
              { eij = e[i]; }
            else if (j == i+1)
              { eij = e[j]; }
            else
              { eij = NAN; }
            if (isnan(eij) || fabs(eij) < tiny)
              { fputs(xzero, wr); }
            else
              { fprintf(wr, fmt, eij); }
          }
        fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
  }

void sym_eigen_test_tools_print_eigenvalues(FILE *wr, char *fmt, uint32_t nev, double d[])
  { fprintf(wr, "\n");
    for (uint32_t i = 0;  i < nev; i++) 
      { fprintf(wr, fmt, d[i]); fputc('\n', wr); }
    fprintf(wr, "\n");
  }

void sym_eigen_test_tools_transform(uint32_t rowsz, uint32_t m, uint32_t n, double R[], double A[], double v[]) 
  {
    /* Set {A[0..m-1,0..n-1]} to {R} times {A}: */
    for (int32_t j = 0;  j < n; j++)
      { /*  Set {v[0..m-1]} to {R} times the column {j} of {A}: */
        for (int32_t i = 0;  i < m; i++)
          { double s = 0.0;
            for (int32_t k = 0;  k < n; k++) { s += R[i*(int32_t)rowsz + k]*A[k*(int32_t)rowsz + j]; }
            v[i] = s;
          }
        /* Now store {v[0..m-1]} into colum {j} of {A}: */
        for (int32_t i = 0;  i < n; i++) { A[i*(int32_t)rowsz + j] = (i < m ? v[i] : 0); }
      }

    /* Set {A[0..m-1,0..m-1]} to {A} times {R} transposed: */
    for (int32_t i = 0;  i < m; i++) /* do 500 */
      { /*  Set {v[0..m-1]} to row {i} of {A} times {R} transposed: */
        for (int32_t j = 0; j < m; j++)
          { double s = 0.0;
            for (int32_t k = 0;  k < n; k++) { s += A[i*(int32_t)rowsz + k]*R[j*(int32_t)rowsz + k]; }
            v[j] = s;
          }
        /* Now store {v[0..m-1]} into row {i} of {A}: */
        for (int32_t j = 0; j < n; j++) { A[i*(int32_t)rowsz + j] = (j < m ? v[j] : 0); }
      }
  }

void sym_eigen_test_tools_check_tridiagonal(uint32_t rowsz, uint32_t n, double A[])
  { double emax = 1.0e-14;
    for (int32_t i = 0;  i < n; i++) 
      { for (int32_t j = 0;  j < n; j++)
          { double Aij = A[i*(int32_t)rowsz + j];
            if ((abs(i-j) > 1) && (fabs(Aij) > emax))
              { fprintf(stderr, "** A[%d,%d] = %24.16e (max = %24.16e)\n", i, j, Aij, emax);
                demand(FALSE, "matrix is not tridiagonal");
              }
          }
      }
  }

void sym_eigen_test_tools_check_diagonal(uint32_t rowsz, uint32_t n, double A[], double d[])
  { double emax = 1.0e-14;
    for (int32_t i = 0;  i < n; i++) 
      { for (int32_t j = 0;  j < n; j++)
          { double Aij = A[i*(int32_t)rowsz + j];
            double Dij = (i == j ? d[i] : 0.0);
            if (fabs(Aij - Dij) > emax)
              { fprintf(stderr, "** A[%d,%d] = %24.16e be %24.16e", i, j, Aij, Dij);
                fprintf(stderr, " (dif = %24.16e, max = %24.16e)\n", Aij - Dij, emax);
                demand(FALSE, "matrix is not {d}-diagonal");
              }
          }
      }
  }

void sym_eigen_test_tools_check_orthogonal(uint32_t rowsz, uint32_t m, uint32_t n, double R[])
  { 
    double emax = 1.0e-13;
    /* Compute {S=R&R^t} and compare with identity: */
    for (int32_t i = 0;  i < m; i++) 
      { for (int32_t j = 0;  j < m; j++)
          { double Sij_exp = (i == j ? 1.0 : 0.0); /* Expected value of {S[i,j]}. */
            double Sij_cmp = 0;
            for (int32_t k = 0; k < n; k++)
              { double Rik = R[i*(int32_t)rowsz + k];
                double Rjk = R[j*(int32_t)rowsz + k];
                Sij_cmp += Rik*Rjk;
              }
            if (fabs(Sij_cmp - Sij_exp) > emax)
              { fprintf(stderr, "** (R*R^t)[%d,%d] = %24.16e (max = %24.16e)\n", i, j, Sij_cmp, emax);
                demand(FALSE, "matrix is not orthogonal");
              }
          }
      }
  }
  
             
