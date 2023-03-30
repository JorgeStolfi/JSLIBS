/* See {dspmat_linsys_GS.h}. */

#define dspmat_linsys_GS_C_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 10:49:20 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>
#include <dspmat.h>

#include <dspmat_linsys_GS.h>

#define debug_level (1)
  /* Define it as {1,2,...} to get increasingly detailed debugging printouts. */

double dspmat_linsys_GS_solve_iteration
  ( double b[],
    dspmat_t *A, 
    double x[],
    double omega,
    double abs_tol, 
    double rel_tol
  )
  {
    if (debug_level >= 1) { fprintf(stderr, "entering %s ...\n", __FUNCTION__); }
    if (debug_level >= 1)
      { fprintf(stderr, "{A.rows} = %d  {A.cols} = %d  {A.ents} = %d\n", A->rows, A->cols, A->ents); }

    demand(A->rows == A->cols, "matrix must be square");
    dspmat_size_t n = A->rows;
    double max_error = 0.000; /* Max error in any coord, rel. to tolerance. */
    dspmat_pos_t posA = 0; /* Scans the entries of {A}. */
    dspmat_index_t row;
    for (row = 0; row < n; row++)
      { /* Compute {sum = SUM{ A[row,j]*x[j] : j \neq row}}, save {Aii=A{row,row]}: */
        double sum = 0.0;
        double Aii = 0.0;
        while (posA < A->ents)
          { dspmat_entry_t *aP = &(A->e[posA]);
            demand(aP->row >= row, "matrix {A} is not sorted");
            assert(aP->row < n);
            assert(aP->col < A->cols);
            if (aP->row > row) { break; }
            if (aP->col == row) 
              { Aii = aP->val; }
            else
              { dspmat_index_t col = aP->col;
                sum += aP->val * x[col];
              }
            posA++;
          }

        /* Compute the error indicator for this row: */
        double error = abs_rel_diff(b[row], sum + Aii*x[row], abs_tol, rel_tol);
 
        /* Record the maximum error: */
        if (fabs(error) > max_error) { max_error = fabs(error); }

        /* Solve for the new value {xNew} of {x[row]}: */
        demand(Aii != 0.0, "matrix {A} has zero in diagonal");
        double xRaw = (b[row] - sum) / Aii;
        demand(isfinite(xRaw), "overflow");
        
        /* Mix the old and new solutions: */
        double xNew = omega*xRaw + (1-omega)*x[row];

        /* Update {x[row]}: */
        x[row] = xNew;
      }
    assert(posA == A->ents);
    if (debug_level >= 1)
      { fprintf(stderr, "  max error = %14.7f\n", max_error);
        if (debug_level >= 2)
          { for (row = 0; row < n; row++)
              { fprintf(stderr, "  x[%4d] = %15.7f = %24.16e\n", row, x[row], x[row]); }
          }
      }
            
    return max_error;
  }

void dspmat_linsys_GS_solve
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

    demand(A->rows == A->cols, "matrix {A} is not square");
    demand(A->cols == nx, "{x} has the wrong length");
    demand(A->rows == nb, "{b} has the wrong length");
    dspmat_size_t nA = A->rows;
    dspmat_index_t col;
    /* Clear the vector {x}: */
    for (col = 0; col < nA; col++) { x[col] = 0.0; }
    /* Gauss-Seidel loop: */
    int32_t iter = 0;
    while (iter < max_iter)
      { /* Gauss-Seidel pass, recomputes {x[0..nA-1]}: */
        if (debug_level >= 2) { fprintf(stderr, "---------------------------\n"); }
        if (debug_level >= 1) { fprintf(stderr, "iteration %3d", iter); }
        double error = dspmat_linsys_GS_solve_iteration(b, A, x, omega, abs_tol, rel_tol);
        if (debug_level >= 2) { fprintf(stderr, "---------------------------\n\n"); }
        iter++;
        if (error <= 1.0) { break; }
      }
    if (debug_level >= 1) { fprintf(stderr, "no convergence in %d iterations", iter); }
  }

void dspmat_inv_mul_GS
  ( dspmat_t *A, 
    dspmat_t *B, 
    dspmat_t *X,
    int32_t max_iter, 
    double omega,
    double abs_tol, 
    double rel_tol
  )
  {
    fprintf(stderr, "entering %s ...\n", __FUNCTION__);
    
    demand(A->rows == A->cols, "matrix is not square");
    demand(A->cols == B->rows, "{B} inconsistent row count");

    X->rows = A->rows;
    X->cols = B->cols;
    
    int32_t nb = B->rows; 
    double *b = notnull(malloc(nb * sizeof(double)), "no mem"); /* Will hold a col of {B}. */
    
    int32_t nx = X->rows;
    double *x = notnull(malloc(nx * sizeof(double)), "no mem"); /* Will hold a col of {X}. */
    
    /* Sort the entries of {B} by column: */
    dspmat_sort_entries(B, 0, +1); 
    if (debug_level > 0)
      { fprintf(stderr, "  matrix A: %d cols %d rows %d entries\n", A->cols, A->rows, A->ents);
        fprintf(stderr, "  matrix B: %d cols %d rows %d entries\n", B->cols, B->rows, B->ents);
      }
    
    /* Process column by column: */
    dspmat_pos_t posB = 0;
    dspmat_pos_t posX = 0;
    int32_t j;
    if (debug_level > 0) { fprintf(stderr, "  processing columns o {B}:"); }
    for (j = 0; j < B->cols; j++)
      { dspmat_pos_t oposB = posB;
        posB = dspmat_extract_col(B, posB, j, b, nb); 
        dspmat_linsys_GS_solve(b, nb, A, x, nx, max_iter, omega, abs_tol, rel_tol);
        dspmat_pos_t oposX = posX; 
        posX = dspmat_add_col(X, posX, j, x, nx);
        if (debug_level > 0)
          { fprintf(stderr, " %u:%u", posB - oposB, posX - oposX); }
      }
    assert(posB == B->ents); 
    dspmat_trim(X, posX);
    if (debug_level > 0) 
      { fprintf(stderr, "\n");
        fprintf(stderr, "  matrix X: %d cols %d rows %d entries\n", X->cols, X->rows, X->ents);
      }
      
    free(x);
    free(b);
  }
