/* Last edited on 2021-06-09 20:41:05 by jstolfi */
/* test of TriRed.c, TriQL.c */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>

#include <sym_eigen.h>

int32_t main(int32_t argn, char **argc);
void filla(int32_t n, double *A, int32_t it);
void prta(int32_t n, double *A);
void prtri(int32_t n, double *d, double *e);
void preval(int32_t n, double *d, int32_t p);
void prevec(int32_t n, double *R, int32_t p);
void appsim(int32_t n, double *R, double *A, double *v);

int32_t main(int32_t argn, char **argc)
  {
    int32_t n = 7;
    double A[n*n], R[n*n], d[n], e[n], v[n];
    int32_t p;
    int32_t test1, test2, absrt;
    test1 = 1;
    test2 = 1;
    absrt = 0;

    int32_t nt = 3; /* Number of test matrices. */
    int32_t it;
    for (it = 0; it < nt; it++)
      { 
        if (test1)
          { 
            fprintf(stderr, "\n");
            fprintf(stderr, "=== test %d ===\n\n", it);
            fprintf(stderr, "\n");
            
            fprintf(stderr, "\n");
            fprintf(stderr, "C - WITHOUT R\n\n");
            fprintf(stderr, "\n");
            filla(n,A,it);
            fprintf(stderr, "input matrix\n");
            prta(n,A);
            syei_tridiagonalize(n,A,d,e,NULL);
            fprintf(stderr, "tridiagonal form\n");
            prtri(n,d,e);
            syei_trid_eigen(n,d,e,NULL,&p,absrt);
            fprintf(stderr, "eigenvalues\n");
            preval(n,d,p);
            fprintf(stderr, "\n");
          }

        if (test2)
          {
            fprintf(stderr, "\n");
            fprintf(stderr, "C - WITH R\n");
            fprintf(stderr, "\n");
            filla(n,A,it);
            fprintf(stderr, "input matrix\n");
            prta(n,A);
            syei_tridiagonalize(n,A,d,e,R);
            fprintf(stderr, "tridiagonal form\n");
            prtri(n,d,e);
            fprintf(stderr, "orthogonal transf (tridiagonal)\n");
            prevec(n,R,n);

            /* Check whether matrix {R} returns tridiagonal form: */
            filla(n,A,it);
            appsim(n,R,A,v);
            fprintf(stderr, "transformed matrix (tridiagonal)\n");
            prta(n,A);

            syei_trid_eigen(n,d,e,R,&p,absrt);
            fprintf(stderr, "eigenvalues\n");
            preval(n,d,p);
            fprintf(stderr, "eigenvectors\n");
            prevec(n,R,p);

            /* Check whether matrix {R} returns diagonal form: */
            filla(n,A,it);
            appsim(n,R,A,v);
            fprintf(stderr, "transformed matrix (diagonal)\n");
            prta(n,A);
            fprintf(stderr, "\n");
          }
      }

    return 0;
  }

void filla(int32_t n, double *A, int32_t it)
  /* Generate a symmetric matrix A(n,n) */
  { 
    int32_t i,j;
    double s, t;
    for (i = 0; i < n; i++) 
      { for (j = 0; j <= i; j++) 
          { double Aij;
            if (it == 0)
              { /* Randomish: */
                s = (i+j-n+1);
                t = (i+1)*(j+1);
                Aij = 1.0/((fabs(s)+1.0)*sqrt(t));
              }
            else if (it == 1)
              { /* Tridiagonal, with 2x2 blocks: */
                Aij = ((i == j) ? 2.0 : (abs(i - j) == 1 ? 1 : 0));
              }
            else if (it == 2)
              { /* Rank 2: */
                double ui0 = 1.0/(1+i);
                double ui1 = 1.0/(1+i*i);
                double uj0 = 1.0/(1+j);
                double uj1 = 1.0/(1+j*j);
                Aij = ui0*uj0 + ui1*uj1;
              }
            else
              { assert(FALSE); }
            A[n*i+j] = Aij;
            A[n*j+i] = A[n*i+j];
          }
      }
  }

void prta(int32_t n, double *A)
  /* prints matrix A(n,n) */
  { 
    int32_t i,j;
    fprintf(stderr, "\n");
    for (i = 0; i < n; i++) 
      { for (j = 0; j < n; j++) 
          { fprintf(stderr, "%10.6f", A[n*i+j]); }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }

void prtri(int32_t n, double *d, double *e)
  /* prints symmetric tridiagonal matrix */
  /* d[0..n-1] is diagonal, e[1..n-1] is sub-diagonal */
  { int32_t i, j;
    fprintf(stderr, "\n");
    for (i = 0; i < n; i++) 
      { for (j = 0; j < i-1; j++)  { fprintf(stderr, "%8s", ""); }
        if (i > 0) fprintf(stderr, "%10.6f", e[i]);
        fprintf(stderr, "%10.6f", d[i]);
        if (i < n-1) fprintf(stderr, "%10.6f", e[i+1]);
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }

void preval(int32_t n, double *d, int32_t p)
  /* prints eigenvalue list d[1..n-1] */
  { int32_t i;
    fprintf(stderr, "\n");
    for (i = 0; i < p; i++) 
      { fprintf(stderr, "%10.6f\n", d[i]); }
    fprintf(stderr, "\n");
  }

void prevec(int32_t n, double *R, int32_t p)
  /* prints eigenvectors R[0..p-1,*] */
  { int32_t i, j;
    fprintf(stderr, "\n");
    for (i = 0; i < p; i++)  
      { for (j = 0; j < n; j++)
          { fprintf(stderr, "%10.6f", R[n*i+j]); }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }    

void appsim(int32_t n, double *R, double *A, double *v)      
  /* applies similarity transformation "R" to "A", */
  /* i.e. computes "R*A*(R^t)", using "v[0..n-1]" as temp storage */
  {
    double s;
    int32_t i, j, k;

    /*set "A" to "R" times "A": */
    for (j = 0; j < n; j++)
      { /*  set "v" to "R" times the column "j" of "A": */
        for (i = 0; i < n; i++)
          { s = 0.0;
            for (k = 0; k < n; k++) { s += R[n*i+k]*A[n*k+j]; }
            v[i] = s;
          }
        /*  now store "v" into colum "j" of "A": */
        for (i = 0; i < n; i++) { A[n*i+j] = v[i]; }
      }

    /*set "A" to "A" times "R" transposed: */
    for (i = 0; i < n; i++) /* do 500 */
      { /*  set "v" to row "i" of "A" times "R" transposed: */
        for (j = 0; j < n; j++)
          { s = 0.0;
            for (k = 0; k < n; k++) { s += A[n*i+k]*R[n*j+k]; }
            v[j] = s;
          }
        /*  now store "v" into row "i" of "A": */
        for (j = 0; j < n; j++) { A[n*i+j] = v[j]; }
      }
  }

  
