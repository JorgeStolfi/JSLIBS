/* See {lsq_array.h} */
/* Last edited on 2014-05-25 21:13:53 by stolfilocal */

#define lsq_array_C_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <rmxn.h>
#include <jsmath.h>
#include <gauss_elim.h>

#include <lsq.h>
#include <lsq_array.h>

int lsq_array_fit
  ( int nt,     /* Number of cases to generate. */
    int nx,     /* Number of independent variables. */
    int nf,     /* Number of dependent variables (functions to fit). */
    double X[], /* Sample values of the independent variables ({nt} by {nx}). */
    double F[], /* Corresponding measured values of the dependent variables ({nt} by {nf}. */
    double W[], /* Corresponding weights ({nt} elements). */
    double U[], /* (OUT) Fitted linear transformation matrix ({nx} by {nt}). */
    bool_t verbose
  )
  {
    double *A = rmxn_alloc(nx,nx);
    lsq_array_compute_matrix(nt, nx, X, W, A);
    double *B = rmxn_alloc(nx,nf);
    lsq_array_compute_rhs(nt, nx, nf, X, F, W, B);
    int rank = lsq_solve_system(nx, nf, A, B, U, verbose);
    free(B);
    free(A);
    return rank;
  }


void lsq_array_compute_matrix(int nt, int nx, double X[], double W[], double A[])
  {
    rmxn_zero(nx, nx, A);
    /* Fill the lower triangular half of {A}: */
    int i, j, k;
    for (k = 0; k < nt; k++)
      { double* Xk = &(X[k*nx]);
        double Wk = (W != NULL ? W[k] : 1);
        for (i = 0; i < nx; i++)
          { for (j = 0; j <= i; j++)
              { double Xkij = Xk[i]*Xk[j];
                A[i*nx + j] += Wk*Xkij;
              }
          }
      }
    /* Replicate the lower half of {A} into the upper half: */
    for (i = 1; i < nx; i++) { for (j = 0; j < i; j++) { A[j*nx + i] = A[i*nx + j]; } }
  }

void lsq_array_compute_rhs(int nt, int nx, int nf, double X[], double F[], double W[], double B[])
  {
    rmxn_zero(nx, nf, B);
    int i, j, k;
    for (k = 0; k < nt; k++)
      { double* Xk = &(X[k*nx]);
        double* Fk = &(F[k*nf]);
        double Wk = (W != NULL ? W[k] : 1);
        for (i = 0; i < nx; i++)
          { for (j = 0; j < nf; j++) 
              { B[i*nf + j] += Wk*Xk[i]*Fk[j]; }
          }
      }
  }

