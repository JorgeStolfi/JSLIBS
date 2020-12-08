#ifndef dspmat_extra_H
#define dspmat_extra_H
/* Extra tools for sparse matrices with {double} entries */

#define dspmat_extra_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-08-17 by J.Stolfi, UNICAMP */
/* Last edited on 2009-01-17 19:29:49 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <dspmat.h>

double dspmat_abs_rel_diff(dspmat_t *A, dspmat_t *B, double abs_tol, double rel_tol);
  /* Computes the maximum difference between pairs of corresponding
    elements of arrays {A,B}, divided by {abs_tol} or {rel_tol} times
    the largest of the two elements. See {abs_rel_diff} in {jsmath.h}
    for details. Both arrays must be sorted by rows, and by columns
    within each row. */
    
void dspmat_normalize_rows(dspmat_t *A, dspmat_t *R);
  /* Scales each non-zero row of {A} by a positive factor
    so that the element with maximum absolute value
    in the row is {±1}. The matrix {A} must be sorted 
    by increasing row index. The result is stored in 
    matrix {R}, whose dimensions are set to those of {A}. */ 

double dspmat_max_element_in_row
  ( dspmat_t *A, 
    dspmat_index_t row, 
    dspmat_pos_t pos
  );
  /* !!! document !!! */

void dspmat_scale_row
  ( dspmat_t *A, 
    dspmat_index_t row, 
    dspmat_pos_t *posAp, 
    double scale, 
    dspmat_t *R, 
    dspmat_pos_t *posRp
  );
  /* !!! document !!! */
  
#endif
