/* See {lsq_array.h} */
/* Last edited on 2019-12-18 17:12:11 by jstolfi */

#define lsq_array_C_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
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

int32_t lsq_array_fit
  ( int32_t nt,     /* Number of data points to generate. */
    int32_t nx,     /* Number of independent variables (argument coordinates per data point). */
    int32_t nf,     /* Number of dependent variables (function samples per data point). */
    double X[], /* Argument coordinates of data points ({nt} by {nx}). */
    double F[], /* Corresponding function samples ({nt} by {nf}. */
    double W[], /* Corresponding weights ({nt} elements). */
    double U[], /* (OUT) Fitted linear transformation matrix ({nx} by {nt}). */
    bool_t verbose
  )
  {
    double *A = rmxn_alloc(nx,nx);
    lsq_array_compute_matrix(nt, nx, X, W, A);
    double *B = rmxn_alloc(nx,nf);
    lsq_array_compute_rhs(nt, nx, nf, X, F, W, B);
    int32_t rank = lsq_solve_system(nx, nf, A, B, 0,NULL,NULL, U,NULL, verbose);
    free(B);
    free(A);
    return rank;
  }

void lsq_array_compute_matrix(int32_t nt, int32_t nx, double X[], double W[], double A[])
  {
    rmxn_zero(nx, nx, A);
    /* Fill the lower triangular half of {A}: */ 
    for (uint32_t k = 0;  k < nt; k++)
      { double* Xk = &(X[k*nx]);
        double Wk = (W != NULL ? W[k] : 1);
        for (uint32_t i = 0;  i < nx; i++)
          { for (uint32_t j = 0;  j <= i; j++)
              { double Xkij = Xk[i]*Xk[j];
                A[i*nx + j] += Wk*Xkij;
              }
          }
      }
    /* Replicate the lower half of {A} into the upper half: */
    for (uint32_t i = 1;  i < nx; i++) 
      { for (uint32_t j = 0;  j < i; j++) 
         { A[j*nx + i] = A[i*nx + j]; } 
      }
  }

void lsq_array_compute_rhs(int32_t nt, int32_t nx, int32_t nf, double X[], double F[], double W[], double B[])
  {
    rmxn_zero(nx, nf, B);
    for (uint32_t k = 0;  k < nt; k++)
      { double* Xk = &(X[k*nx]);
        double* Fk = &(F[k*nf]);
        double Wk = (W != NULL ? W[k] : 1);
        for (uint32_t i = 0;  i < nx; i++)
          { for (uint32_t j = 0;  j < nf; j++) 
              { B[i*nf + j] += Wk*Xk[i]*Fk[j]; }
          }
      }
  }

