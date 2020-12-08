/* See {dspmat_extra.h}. */

#define dspmat_extra_C_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2009-08-31 21:42:54 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <jsmath.h>
#include <dspmat.h>

#include <dspmat_linsys_GS.h>
#include <dspmat_linsys_ALT.h>

#include <dspmat_extra.h>

#define debug_level (1)
  /* Define it as {1,2,...} to get increasingly detailed debugging printouts. */

double dspmat_abs_rel_diff(dspmat_t *A, dspmat_t *B, double abs_tol, double rel_tol)
  {
    double max_error = 0.0;
    
    auto void abs_rel_diff_entries
      ( dspmat_index_t row, 
        dspmat_index_t col, 
        double Aij, 
        double Bij
      );
    
    void abs_rel_diff_entries
      ( dspmat_index_t row, 
        dspmat_index_t col, 
        double Aij, 
        double Bij
      )
      { double error = abs_rel_diff(Aij, Bij, abs_tol, rel_tol);
        if (fabs(error) > max_error) { max_error = fabs(error); }
      }
    
    dspmat_merge(A, B, &abs_rel_diff_entries);
    
    return max_error;
  }

void dspmat_normalize_rows(dspmat_t *A, dspmat_t *R)
  {
    dspmat_size_t rows = A->rows;
    dspmat_size_t cols = A->cols;
    
    R->rows = rows;
    R->cols = cols;
    
    dspmat_pos_t posA = 0; /* Scans the entries of {A}. */
    dspmat_pos_t posR = 0; /* Scans the entries of {A}. */
    dspmat_index_t row;
    for (row = 0; row < rows; row++)
      { double amax = dspmat_max_element_in_row(A, row, posA);
        double scale = (fabs(amax) == 0.0 ? 1.0 : 1.0/amax);
        /* Scale the elements: */
        dspmat_scale_row(A, row, &posA, scale, R, &posR);
      }
    demand(posA == A->ents, "matrix {A} was not sorted by rows");
    dspmat_trim(R, posR);
  }
  
double dspmat_max_element_in_row
  ( dspmat_t *A, 
    dspmat_index_t row, 
    dspmat_pos_t pos
  )
  {
    double vmax = 0.0; /* Element with max abs value. */
    
    auto void get_max(dspmat_index_t row, dspmat_index_t col, double val);
    void get_max(dspmat_index_t row, dspmat_index_t col, double val)
      { if (fabs(val) > fabs(vmax)) { vmax = val; } }
      
    dspmat_scan_row(A, pos, row, get_max);
    return vmax;
  }
 
void dspmat_scale_row
  ( dspmat_t *A, 
    dspmat_index_t row, 
    dspmat_pos_t *posAp, 
    double scale, 
    dspmat_t *R, 
    dspmat_pos_t *posRp
  )
  {
    dspmat_pos_t posA = (*posAp);
    dspmat_pos_t posR = (*posRp);
    
    auto void scale_elem(dspmat_index_t rowA, dspmat_index_t colA, double valA);
    void scale_elem(dspmat_index_t rowA, dspmat_index_t colA, double valA)
      { assert(rowA == row);
        double valR = scale * valA;
        posR = dspmat_add_element(R, posR, rowA, colA, valR);
      }
      
    posA = dspmat_scan_row(A, posA, row, scale_elem);
    
    (*posAp) = posA;
    (*posRp) = posR;
  }
