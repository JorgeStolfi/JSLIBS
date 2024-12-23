/* See {dspmat_linsys_SUBGE.h}. */

#define dspmat_linsys_SUBGE_C_COPYRIGHT "Copyright � 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2024-12-05 10:40:22 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <rn.h>
#include <affirm.h>
#include <dspmat.h>
#include <dspmat_extra.h>

#include <dspmat_linsys_SUBGE.h>

#define debug_level (1)
  /* Define it as {1,2,...} to get increasingly detailed debugging printouts. */

void dspmat_linsys_SUBGE_precond(dspmat_t *A, dspmat_t *B)
  {
    assert(A->cols == A->rows);
    dspmat_size_t n = A->cols;
    
    assert(B->rows == A->rows);
    dspmat_size_t p = B->cols;
    
    /* Sort the entries of {A,B} by row: */
    dspmat_sort_entries(A, +2, +1); 
    dspmat_sort_entries(B, +2, +1); 
    
    /* Allocate temporary work matrices: */
    dspmat_t AT = dspmat_new(n, n, A->ents);
    dspmat_t BT = dspmat_new(n, p, B->ents);
    
    /* Effective entry counts: */
    dspmat_count_t entsA = A->ents; /* Counts used entries of {A}. */
    dspmat_count_t entsB = B->ents; /* Counts used entries of {B}. */
    dspmat_count_t entsAT = 0; /* Counts used entries of {AT}. */
    dspmat_count_t entsBT = 0; /* Counts used entries of {BT}. */
    
    /* Phase I: Normalize each row so that the largest {A} element is 1. */
    dspmat_linsys_SUBGE_norm_rows(A, entsA, B, entsB, &AT, *entsAT, &BT, *entsBT);
    
    /* Phase II: Make {A} almost upper triangular. */
    /* Elements in the diagonal are 1. */
    /* All elements below the diagonal are less than {1/n}.*/
    /* All elements above the diagonal are at most 1. */
    
    for (i = 0; i < n; i++)
      { 
        /* Choose a row {k} among rows {i..n-1} of {AT}: */
        for (k = i; k < n; k++)
          { 
        /* Move row {k} to {i}. */
        /* Eliminate row {i} from rows {i+1..n-1} if needed. */
      }
  }

void dspmat_linsys_SUBGE_solve
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
      { fprintf(stderr, 
          "{A.rows} = %d  {A.cols} = %d  {A.ents} = %d\n", 
          A->rows, A->cols, A->ents);
      }
    
    demand(A->rows == A->cols, "matrix {A} is not square");
    demand(A->cols == nx, "{x} has the wrong length");
    demand(A->rows == nb, "{b} has the wrong length");
    dspmat_size_t n = A->rows;

    dspmat_sort_entries(A, +2, +1); /* Sort {A} entries by row, then col. */
    
    /* Initialize the output {x} with the {b} matrix: */
    rn_copy(n, b, x);
    
    /* Scratch arrays and vectors: */
    double *y = notnull(malloc(n * sizeof(double)), "no mem"); 
    dspmat_t Q = dspmat_new(0,0,0);
            
    /* Create a matrix {M = I-A}: */
    dspmat_t M = dspmat_new(0,0,0);
    dspmat_identity(&Q, n);
    dspmat_mix(1.0, A, -1.0, &Q, &M);
    
    /* Bootstrap loop: */
    int32_t iter = 0;
    while (TRUE)
      { if (debug_level >= 2) { fprintf(stderr, "---------------------------\n"); }
        if (debug_level >= 1) { fprintf(stderr, "iteration %3d", iter); }
        if (debug_level >= 1) { fprintf(stderr, "  {M.ents} = %d", M.ents); }
        /* Check the residual {B - A*X}: */
        dspmat_map_col(A, x, n, y, n);
        double error = rn_abs_rel_diff(n, b, y, abs_tol, rel_tol);
        if (debug_level >= 1) { fprintf(stderr, "  max error = %14.7f\n", error); }
        if (error <= 1.0) { break; }
      
        /* Bootstrap pass: */
        dspmat_map_col(&M, x, n, y, n);
        rn_add(n, x, y, x);
        
        dspmat_mul(&M, &M, &Q);
        dspmat_copy(&Q, &M);
        if (debug_level >= 2) { fprintf(stderr, "---------------------------\n\n"); }
        iter++;
        if (iter >= max_iter)
          { if (debug_level >= 1) 
              { fprintf(stderr, "no convergence in %d iterations", iter); }
            break;
          }
      }
      
    /* Reclaim storage: */      
    free(y);
    dspmat_trim(&Q, 0);
    dspmat_trim(&M, 0);
  }

void dspmat_inv_mul_SUBGE
  ( dspmat_t *A,
    dspmat_t *B,
    dspmat_t *X,
    int32_t max_iter, 
    double abs_tol, 
    double rel_tol
  )
  {
    if (debug_level >= 1) { fprintf(stderr, "entering %s ...\n", __FUNCTION__); }
    if (debug_level >= 1)
      { fprintf(stderr, "{A.rows} = %d  {A.cols} = %d  {A.ents} = %d\n", A->rows, A->cols, A->ents); }

    demand(A->rows == A->cols, "matrix {A} is not square");
    demand(A->rows == B->rows, "{A.cols} does not match {B.rows}");
    dspmat_size_t n = A->rows;

    dspmat_sort_entries(A, +2, +1); /* Sort {A} entries by row, then col. */
    dspmat_sort_entries(B, +2, +1); /* Sort {B} entries by row, then col. */
    
    /* Initialize the output {X} with the {B} matrix: */
    dspmat_copy(B,X);
    
    /* Scratch arrays: */
    dspmat_t Y = dspmat_new(0,0,0); 
    dspmat_t Z = dspmat_new(0,0,0);
    dspmat_t Q = dspmat_new(0,0,0);
            
    /* Create a matrix {M = I-A}: */
    dspmat_t M = dspmat_new(0,0,0);
    dspmat_identity(&Q, n);
    dspmat_mix(1.0, A, -1.0, &Q, &M);
    
    /* Bootstrap loop: */
    int32_t iter = 0;
    while (TRUE)
      { if (debug_level >= 2) { fprintf(stderr, "---------------------------\n"); }
        if (debug_level >= 1) { fprintf(stderr, "iteration %3d", iter); }
        /* Check the residual {B - A*X}: */
        dspmat_mul(A, X, &Y);
        double error = dspmat_abs_rel_diff(B, &Y, abs_tol, rel_tol);
        if (debug_level >= 1) { fprintf(stderr, "  max error = %14.7f", error); }
        if (error <= 1.0) { break; }
      
        /* Bootstrap pass: */
        dspmat_mul(&M, X, &Y);
        dspmat_mix(1.0, X, 1.0, &Y, &Z);
        dspmat_copy(&Z, X); 
        
        dspmat_mul(&M, &M, &Q);
        dspmat_copy(&Q, &M);
        if (debug_level >= 1) { fprintf(stderr, "  {X.ents} = %d", X->ents); }
        if (debug_level >= 1) { fprintf(stderr, "  {M.ents} = %d\n", M.ents); }
        if (debug_level >= 2) { fprintf(stderr, "---------------------------\n\n"); }
        iter++;
        if (iter >= max_iter)
          { if (debug_level >= 1) 
              { fprintf(stderr, "no convergence in %d iterations", iter); }
            break;
          }
      }
      
    /* Reclaim storage: */      
    dspmat_trim(&Y, 0);
    dspmat_trim(&Z, 0);
    dspmat_trim(&Q, 0);
    dspmat_trim(&M, 0);
  }
  
