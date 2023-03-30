/* See {dspmat_linsys_ALT_solve.h}. */

#define dspmat_linsys_ALT_solve_C_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 10:49:40 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <dspmat.h>
#include <jsmath.h>

#include <dspmat_linsys_ALT.h>

#define debug_level (1)
  /* Define it as {1,2,...} to get increasingly detailed debugging printouts. */

double dspmat_linsys_ALT_compute_residuals
  ( double b[],
    dspmat_t *A, 
    double x[],
    double d[],
    double abs_tol, 
    double rel_tol
  )
  {
    double max_error = 0.0; /* Max residual relative to {abs_tol,rel_tol*b[i]}. */

    /* Compute {d[i] = SUM{A[i,j]*x[j] : j}} for all {i}: */
    dspmat_index_t row;
    dspmat_pos_t pos;
    for (row = 0; row < A->rows; row++) { d[row] = 0.0; }
    for (pos = 0; pos < A->ents; pos++)
      { dspmat_entry_t *aP = &(A->e[pos]);
        d[aP->row] += aP->val * x[aP->col];
      }
    for (row = 0; row < A->rows; row++) 
      { /* Compute the error indicator for {b[row]}: */
        double error = abs_rel_diff(b[row], d[row], abs_tol, rel_tol);

        /* Record the maximum error: */
        if (fabs(error) > max_error) { max_error = fabs(error); }

        /* Set {d[i]} to the residual {b[i] - SUM{A[i,j]*x[j] : j}}: */
        d[row] = b[row] - d[row];
      }

    if (debug_level >= 1)
      { fprintf(stderr, "  max error = %14.7f\n", max_error);
        if (debug_level >= 2)
          { for (row = 0; row < A->rows; row++)
              { fprintf(stderr, "  d[%4d] = %15.7f = %24.16e\n", row, d[row], d[row]); }
          }
      }
                
    return max_error;
  }

void dspmat_linsys_ALT_recompute_solution
  ( double b[],
    dspmat_t *A,
    double x[], 
    double d[],
    double c[],
    double r[],
    double M[]
  )
  {
    dspmat_index_t col;
    dspmat_pos_t pos;
    
    /* Compute {M[pos] = M[pos] + A[pos]*d[row]/r[row]}: */
    for (pos = 0; pos < A->ents; pos++) 
      { dspmat_entry_t *aP = &(A->e[pos]);
        M[pos] = aP->val * (x[aP->col] + d[aP->row]/r[aP->row]);
      }
      
    /* Compute {x[col] = SUM{M[pos] : row}/c[col]}: */
    for (col = 0; col < A->cols; col++) { x[col] = 0.0; }
    for (pos = 0; pos < A->ents; pos++) { x[A->e[pos].col] += M[pos]; }
    for (col = 0; col < A->cols; col++) { x[col] /= c[col]; }
    
    if (debug_level >= 2)
      { for (col = 0; col < A->cols; col++)
          { fprintf(stderr, "  x[%4d] = %15.7f = %24.16e\n", col, x[col], x[col]); }
      }
  }
  
#define TINY_SUM (1.0e-300)

void dspmat_linsys_ALT_solve
  ( double b[],
    dspmat_size_t nb,
    dspmat_t *A, 
    double x[], 
    dspmat_size_t nx,
    int32_t max_iter, 
    double omega,
    double abs_tol, 
    double rel_tol
  )
  {
    if (debug_level >= 1) { fprintf(stderr, "entering %s ...\n", __FUNCTION__); }
    if (debug_level >= 1)
      { fprintf(stderr, "{A.rows} = %d  {A.cols} = %d  {A.ents} = %d\n", A->rows, A->cols, A->ents); }

    demand(A->rows == A->cols, "matrix is not square");
    demand(A->cols == nx, "{x} has the wrong length");
    demand(A->rows == nb, "{b} has the wrong length");
    dspmat_size_t n = A->rows;

    double *c = notnull(malloc(A->cols * sizeof(double)), "no mem");  /* {A} col sums. */
    double *r = notnull(malloc(A->rows * sizeof(double)), "no mem");  /* {A} row sums. */
    double *d = notnull(malloc(A->rows * sizeof(double)), "no mem");  /* Row discrepancies. */
    double *M = notnull(malloc(A->ents * sizeof(double)), "no mem");  /* Terms of {A*x}. */
    
    /* Compute the col sums {c[0..n-1]} and the row sums {r[0..n-1]}: */
    dspmat_index_t row, col;
    dspmat_pos_t pos;
    for (col = 0; col < n; col++) { c[col] = 0.0; }
    for (row = 0; row < n; row++) { r[row] = 0.0; }
    for (pos = 0; pos < A->ents; pos++) 
      { dspmat_entry_t *aP = &(A->e[pos]);
        r[aP->row] += aP->val;
        c[aP->col] += aP->val;
      }

    /* If any {c[col]} is zero, set it to {TINY_SUM} so that {x[col]} will be 0 instead of {NAN}: */
    for (col = 0; col < n; col++) { if (c[col] == 0.0) { c[col] = TINY_SUM; } }

    /* If any {r[row]} is zero, set it to {TINY_SUM} so that {x[col]} will be 0 instead of {NAN}: */
    for (row = 0; row < n; row++) 
      { if (r[row] == 0.0) 
          { demand(fabs(b[row]) < abs_tol, "system has no solution"); 
            r[row] = TINY_SUM;
          }
      }
    
    /* Clear the vector {x}: */
    for (col = 0; col < n; col++) { x[col] = 0.0; }
    /* Iteration loop: */
    int32_t iter = 0;
    while (TRUE)
      { /* Check number of iterations: */
        if (iter >= max_iter) 
          { if (debug_level >= 1) { fprintf(stderr, "no convergence in %d iterations", iter); }
            break;
          }
        /* Compute the residuals {d[0..n-1]}: */
        if (debug_level >= 2) { fprintf(stderr, "---------------------------\n"); }
        if (debug_level >= 1) { fprintf(stderr, "iteration %3d", iter); }
        double err = dspmat_linsys_ALT_compute_residuals(b, A, x, d, abs_tol, rel_tol);
        if (debug_level >= 2) { fprintf(stderr, "---------------------------\n"); }
        if (err <= 1.0) { break; }
        /* Recompute the solution {x[0..n-1]}: */
        if (debug_level >= 2) { fprintf(stderr, "---------------------------\n"); }
        dspmat_linsys_ALT_recompute_solution(b, A, x, d, c, r, M);
        if (debug_level >= 2) { fprintf(stderr, "---------------------------\n\n"); }
        iter++;
      }
    /* Reclaim auxiliary storage: */
    free(M);
    free(d);
    free(r);
    free(c);
  }
