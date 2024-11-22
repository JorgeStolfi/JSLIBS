/* Last edited on 2024-11-21 21:36:10 by stolfi */
/* test of TriRed.c, TriQL.c */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>

#include <sym_eigen.h>

int32_t main(int32_t argn, char **argc);
void filla(uint32_t n, double *A, uint32_t it);
void prta(uint32_t n, double *A);
void prtri(uint32_t n, double *d, double *e);
void preval(uint32_t n, double *d, uint32_t p);
void prevec(uint32_t n, double *R, uint32_t p);
void appsim(uint32_t n, double *R, double *A, double *v);

int32_t main(int32_t argn, char **argc)
  {
    uint32_t n = 7;
    double A[n*n], R[n*n], d[n], e[n], v[n];
    uint32_t absrt = 0;

    uint32_t nt = 3; /* Number of test matrices. */
    for (int32_t it =  0; it < nt; it++)
      { 
        fprintf(stderr, "\n");
        fprintf(stderr, "C - WITH R\n");
        fprintf(stderr, "\n");
        filla(n,A,it);
        fprintf(stderr, "C - input matrix\n");
        prta(n,A);
        sym_eigen_tridiagonalize(n,A,d,e,R);
        fprintf(stderr, "C - tridiagonal form\n");
        prtri(n,d,e);
        fprintf(stderr, "C - orthogonal transf (tridiagonal)\n");
        prevec(n,R,n);

        /* Check whether matrix {R} returns tridiagonal form: */
        filla(n,A,it);
        appsim(n,R,A,v);
        fprintf(stderr, "C - transformed matrix (tridiagonal)\n");
        prta(n,A);

        uint32_t p;
        sym_eigen_trid_eigen(n,d,e,R,&p,absrt);
        fprintf(stderr, "C - eigenvalues\n");
        preval(n,d,p);
        fprintf(stderr, "C - eigenvectors\n");
        prevec(n,R,p);

        /* Check whether matrix {R} returns diagonal form: */
        filla(n,A,it);
        appsim(n,R,A,v);
        fprintf(stderr, "C - transformed matrix (diagonal)\n");
        prta(n,A);
        fprintf(stderr, "\n");
      }

    return 0;
  }

void filla(uint32_t n, double *A, uint32_t it)
  /* Generate a symmetric matrix A(n,n) */
  { 
    for (int32_t i =  0; i < n; i++) 
      { for (int32_t j =  0; j <= i; j++) 
          { double Aij;
            if (it == 0)
              { /* Randomish: */
                double s = 1.0 + i + j - n;
                double t = (i+1.0)*(j+1.0);
                Aij = 0.0001/((fabs(s)+1.0)*sqrt(t));
              }
            else if (it == 1)
              { /* Tridiagonal, with 2x2 blocks: */
                Aij = ((i == j) ? 2.0 : (abs((int32_t)i - (int32_t)j) == 1 ? 1 : 0));
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

void prta(uint32_t n, double *A)
  /* prints matrix A(n,n) */
  { 
    fprintf(stderr, "\n");
    for (int32_t i =  0; i < n; i++) 
      { for (int32_t j =  0; j < n; j++) 
          { fprintf(stderr, "%14.10f", A[n*i+j]); }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }

void prtri(uint32_t n, double *d, double *e)
  /* prints symmetric tridiagonal matrix */
  /* d[0..n-1] is diagonal, e[1..n-1] is sub-diagonal */
  { fprintf(stderr, "\n");
    for (int32_t i =  0; i < n; i++) 
      { for (int32_t j =  0; j < i-1; j++)  { fprintf(stderr, "%14s", ""); }
        if (i > 0) fprintf(stderr, "%14.10f", e[i]);
        fprintf(stderr, "%14.10f", d[i]);
        if (i < n-1) fprintf(stderr, "%14.10f", e[i+1]);
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }

void preval(uint32_t n, double *d, uint32_t p)
  /* prints eigenvalue list d[1..n-1] */
  { fprintf(stderr, "\n");
    for (int32_t i =  0; i < p; i++) 
      { fprintf(stderr, "%14.10f\n", d[i]); }
    fprintf(stderr, "\n");
  }

void prevec(uint32_t n, double *R, uint32_t p)
  /* prints eigenvectors R[0..p-1,*] */
  { fprintf(stderr, "\n");
    for (int32_t i =  0; i < p; i++)  
      { for (int32_t j =  0; j < n; j++)
          { fprintf(stderr, "%14.10f", R[n*i+j]); }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }    

void appsim(uint32_t n, double *R, double *A, double *v)      
  /* applies similarity transformation "R" to "A", */
  /* i.e. computes "R*A*(R^t)", using "v[0..n-1]" as temp storage */
  {
    /*set "A" to "R" times "A": */
    for (int32_t j =  0; j < n; j++)
      { /*  set "v" to "R" times the column "j" of "A": */
        for (int32_t i =  0; i < n; i++)
          { double s = 0.0;
            for (int32_t k =  0; k < n; k++) { s += R[n*i+k]*A[n*k+j]; }
            v[i] = s;
          }
        /*  now store "v" into colum "j" of "A": */
        for (int32_t i =  0; i < n; i++) { A[n*i+j] = v[i]; }
      }

    /*set "A" to "A" times "R" transposed: */
    for (int32_t i =  0; i < n; i++) /* do 500 */
      { /*  set "v" to row "i" of "A" times "R" transposed: */
        for (int32_t j =  0; j < n; j++)
          { double s = 0.0;
            for (int32_t k =  0; k < n; k++) { s += A[n*i+k]*R[n*j+k]; }
            v[j] = s;
          }
        /*  now store "v" into row "i" of "A": */
        for (int32_t j = 0; j < n; j++) { A[n*i+j] = v[j]; }
      }
  }

  
