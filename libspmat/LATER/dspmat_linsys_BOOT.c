/* See {dspmat_linsys_BOOT.h}. */

#define dspmat_linsys_BOOT_C_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2009-08-31 21:44:36 by stolfi */

#define _GNU_SOURCE
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

#include <dspmat_linsys_BOOT.h>

#define debug_level (1)
  /* Define it as {1,2,...} to get increasingly detailed debugging printouts. */

void dspmat_linsys_BOOT_solve
  ( double b[],
    dspmat_size_t nb,
    dspmat_t *A, 
    double x[], 
    dspmat_size_t nx,
    int max_iter, 
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
    int iter = 0;
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

void dspmat_inv_mul_BOOT
  ( dspmat_t *A,
    dspmat_t *B,
    dspmat_t *X,
    int max_iter, 
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
    int iter = 0;
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
  
