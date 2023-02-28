/* See {tf_matrix.h}. */
/* Last edited on 2023-02-25 16:13:18 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <rmxn.h>
#include <rn.h>

#include <tf_matrix.h>

void tf_solve_system_mxn (mat_rm_t M, rm_t U, rm_t b, double wgts[])
{
    mat_rm_t Mt, MtM, iMtM, Mdag;
    
    Mt = tf_transpose_mat_rm (M);
    
    tf_apply_weigths (Mt, wgts);

    MtM = tf_alloc_mat_rm (Mt->nrows, M->ncols);

    rmxn_mul (Mt->nrows, Mt->ncols, M->ncols, Mt->c, M->c, MtM->c);

    iMtM = tf_alloc_mat_rm (MtM->nrows, MtM->ncols);

    (void) rmxn_inv (MtM->nrows, MtM->c, iMtM->c);
    
    Mdag = tf_alloc_mat_rm (iMtM->nrows, Mt->ncols);
    
    rmxn_mul (iMtM->nrows, iMtM->ncols, Mt->ncols, iMtM->c, Mt->c, Mdag->c);

    rmxn_map_col (Mdag->nrows, Mdag->ncols, Mdag->c, b->c, U->c);

    tf_free_mat_rm_structure (Mt);
    tf_free_mat_rm_structure (MtM);
    tf_free_mat_rm_structure (iMtM);
    tf_free_mat_rm_structure (Mdag);
}

void tf_solve_system_mxn_mxp (mat_rm_t M, mat_rm_t U, mat_rm_t B, double wgts[])
{
  assert(M->ncols == U->nrows);
  assert(M->nrows == B->nrows);
  assert(U->ncols == B->ncols);

    mat_rm_t Mt, MtM, iMtM, Mdag;
    
    Mt = tf_transpose_mat_rm (M);
    
    tf_apply_weigths (Mt, wgts);

    MtM = tf_alloc_mat_rm (Mt->nrows, M->ncols);

    rmxn_mul (Mt->nrows, Mt->ncols, M->ncols, Mt->c, M->c, MtM->c);

    iMtM = tf_alloc_mat_rm (MtM->nrows, MtM->ncols);

    (void) rmxn_inv (MtM->nrows, MtM->c, iMtM->c);
    
    Mdag = tf_alloc_mat_rm (iMtM->nrows, Mt->ncols);
    
    rmxn_mul (iMtM->nrows, iMtM->ncols, Mt->ncols, iMtM->c, Mt->c, Mdag->c);

    rmxn_mul (Mdag->nrows, Mdag->ncols, B->ncols, Mdag->c, B->c, U->c);

    tf_free_mat_rm_structure (Mt);
    tf_free_mat_rm_structure (MtM);
    tf_free_mat_rm_structure (iMtM);
    tf_free_mat_rm_structure (Mdag);
}

mat_rm_t tf_transpose_mat_rm (mat_rm_t m)
{
    int32_t i, j;

    mat_rm_t r = tf_alloc_mat_rm (m->ncols, m->nrows);

    for (i = 0; i < m->nrows; i++) {
        for (j = 0; j < m->ncols; j++) {
            r->c[j*m->nrows + i] = m->c[i*m->ncols + j];
        }
    }

    return r;
}

void tf_apply_weigths (mat_rm_t m, double wgts[])
{
  if (wgts == NULL) { return; }
    int32_t i, j;

    for (i = 0; i < m->nrows; i++) {
        for (j = 0; j < m->ncols; j++) {
            m->c[i*m->ncols + j] *= wgts[j];
        }
    }
}

void tf_print_mat_rm (mat_rm_t m)
{
    int32_t i, j;

    for (i = 0; i < m->nrows; i++) {
        for (j = 0; j < m->ncols; j++) {
            fprintf (stdout, "  %5.15lf\n", m->c[i*m->ncols + j]);
        }
        fprintf (stdout, "\n");
    }
}


rmxn_t tf_alloc_rmxn (int32_t nrows, int32_t ncols) 
{
    int32_t row; /* Array of rows */

    rmxn_t m = (rmxn_t)malloc(sizeof(struct _rmxn_t));

    /* Allocates rows */
    m->nrows = nrows;
    m->ncols = ncols;
    m->c  = (double **)malloc(nrows * sizeof (double *));

    if (m->c == NULL ){
        fprintf( stderr, "error: allocating double matrix" );
        return NULL;
    }

    /* Allocates columns */
    for (row = 0; row < nrows; row++){
        m->c[row] = rn_alloc(ncols);
        if (m->c[row] == NULL ){
            fprintf( stderr, "error: allocating array" );
            return NULL;
        }
    }
    return m; /* Returns bidimensional array */
}

void tf_free_rmxn_structure (rmxn_t m)
{
    int32_t row;

    for (row = 0; row < m->nrows; row++){
        free(m->c[row]);
    }

    free(m->c);    
    free(m);    
}


rm_t tf_alloc_rm (int32_t size)
{
    rm_t v = (rm_t)malloc(sizeof(struct _rm_t));

    v->size = size;
    v->c  = rn_alloc(size);

    return v;
}

void tf_free_rm_structure (rm_t v)
{
    free(v->c);
    free(v); 
}

void tf_copy_vector_rm (rm_t source, rm_t target)
{
    int32_t i;

    target->size = source->size;

    for (i = 0; i < source->size; i++)
        target->c[i] = source->c[i];
}

mat_rm_t tf_alloc_mat_rm (int32_t nrows, int32_t ncols)
{
   mat_rm_t m = (mat_rm_t)malloc(sizeof(struct _mat_rm_t));

   m->nrows = nrows;
   m->ncols = ncols;

   m->c  = rmxn_alloc(nrows,ncols);

   return m;
}

void tf_free_mat_rm_structure (mat_rm_t m)
{
    free(m->c);
    free(m);
}
