/* Functions for linear algebra and matrix manipulation. */
/* Last edited on 2011-05-15 00:39:26 by stolfi */

#ifndef tf_matrix_H
#define tf_matrix_H

#include <stdio.h>
#include <ctype.h>
#include <r3.h> 
#include <r2.h> 
#include <r4x4.h> 

typedef struct _rmxn_t {
    int  nrows;
    int  ncols;
    double **c;             
} *rmxn_t;

typedef struct _mat_rm_t {
    int  nrows;
    int  ncols;
    double *c;             
} *mat_rm_t;

typedef struct _rm_t {
    int size; 
    double *c;             
} *rm_t;

void                     tf_copy_vector_rm (rm_t source, rm_t target);

mat_rm_t                 tf_alloc_mat_rm (int nrows, int ncols);

void                     tf_free_mat_rm_structure (mat_rm_t m);

rm_t                     tf_alloc_rm (int size);

void                     tf_free_rm_structure (rm_t v);

rmxn_t                   tf_alloc_rmxn (int nrows, int ncols);

void                     tf_free_rmxn_structure (rmxn_t m);


void                     tf_solve_system_mxn (mat_rm_t M, rm_t a, rm_t b, double wgts[]);

void tf_solve_system_mxn_mxp (mat_rm_t M, mat_rm_t U, mat_rm_t B, double wgts[]);
/* Solves {M U = B} for {U} in least squares sense. */

void                     tf_apply_weigths (mat_rm_t m, double wgts[]);

mat_rm_t                 tf_transpose_mat_rm (mat_rm_t m);

void                     tf_print_mat_rm (mat_rm_t m);

#endif
