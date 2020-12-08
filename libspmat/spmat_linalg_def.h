#ifndef spmat_linalg_def_H
#define spmat_linalg_def_H
/* The ugly entrails of {spmat_linalg.h}. */

#define spmat_linalg_def_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2017-01-02 21:39:29 by jstolfi */

/* These inclusions are necessary if this file is included or compiled on its own: */
#include <stdint.h>
#include <stdlib.h>
#include <spmat.h>
#include <spmat_linalg.h>
  
#define spmat_DECLARE_identity(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_identity(MATRIX_TYPE *M, spmat_size_t size)

#define spmat_DECLARE_mix(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_mix \
    ( ELEM_TYPE a, \
      MATRIX_TYPE *A, \
      ELEM_TYPE b, \
      MATRIX_TYPE *B, \
      MATRIX_TYPE *C \
    )

#define spmat_DECLARE_mul(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_mul(MATRIX_TYPE *A, MATRIX_TYPE *B, MATRIX_TYPE *C)

#define spmat_DECLARE_map_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_map_col \
    ( MATRIX_TYPE *M, \
      ELEM_TYPE a[], \
      spmat_size_t na, \
      ELEM_TYPE b[], \
      spmat_size_t nb \
    )

#define spmat_DECLARE_map_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_map_row \
    ( ELEM_TYPE a[], \
      spmat_size_t na, \
      MATRIX_TYPE *M, \
      ELEM_TYPE b[], \
      spmat_size_t nb \
    )

#define spmat_LINALG_DEF(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_identity(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_mix(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_mul(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_map_col(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_map_row(MATRIX_TYPE,PREFIX,ELEM_TYPE)

/* IMPLEMENTATION MACROS */

#define spmat_IMPLEMENT_identity(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_identity(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { M->rows = M->cols = size; \
      PREFIX##_trim(M, size); \
      spmat_pos_t pos = \
        PREFIX##_fill_diagonal(M, 0, 0,0, PREFIX##_elem_one, size); \
      assert(pos == size); \
      assert(M->ents == size);  \
    }

#define spmat_IMPLEMENT_mix(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_mix(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { spmat_pos_t posC = 0; \
       \
      auto void mix_entries \
        ( spmat_index_t row, \
          spmat_index_t col, \
          const ELEM_TYPE Aij, \
          const ELEM_TYPE Bij \
        ); \
       \
      void mix_entries \
        ( spmat_index_t row, \
          spmat_index_t col, \
          const ELEM_TYPE Aij, \
          const ELEM_TYPE Bij \
        ) \
        { ELEM_TYPE val; \
          if (PREFIX##_elem_is_trivial(Aij)) \
            { val = PREFIX##_elem_mul(b, Bij); } \
          else if (PREFIX##_elem_is_trivial(Bij)) \
            { val = PREFIX##_elem_mul(a, Aij); } \
          else \
            { ELEM_TYPE va = PREFIX##_elem_mul(a, Aij); \
              ELEM_TYPE vb = PREFIX##_elem_mul(b, Bij); \
              val = PREFIX##_elem_add(va, vb); \
            } \
          if (! PREFIX##_elem_is_trivial(val)) \
            { PREFIX##_expand(C, posC); \
              PREFIX##_entry_t *eC = &(C->e[posC]); \
              eC->row = row; eC->col = col; eC->val = val; \
              posC++; \
            } \
        } \
       \
      C->rows = A->rows; \
      C->cols = A->cols; \
      PREFIX##_merge(A, B, &mix_entries); \
      PREFIX##_trim(C, posC); \
    }

#define spmat_IMPLEMENT_mul(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_mul(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(A->cols == B->rows, "incompatible row/column counts"); \
      C->rows = A->rows; C->cols = B->cols; \
      ELEM_TYPE *av = malloc(A->cols * sizeof(ELEM_TYPE)); /* A row of {A} */ \
      ELEM_TYPE *cv = malloc(B->cols * sizeof(ELEM_TYPE)); /* A row of {B} */ \
      spmat_pos_t posA = 0; /* Start of next row of {A}. */ \
      spmat_pos_t posC = 0; /* Next free eentry in {C}. */ \
      while (posA < A->ents) \
        { spmat_index_t row = A->e[posA].row; \
          posA = PREFIX##_extract_row(A, posA, row, av, A->cols); \
          PREFIX##_map_row(av, A->cols, B, cv, B->cols); \
          posC = PREFIX##_add_row(C, posC, row, cv, B->cols); \
        } \
      PREFIX##_trim(C, posC); \
      free(av); free(cv); \
    }

#define spmat_IMPLEMENT_map_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_map_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(M->cols == na, "incompatible column counts"); \
      demand(M->rows == nb, "incompatible row counts"); \
      spmat_index_t i; \
      for (i = 0; i < nb; i++) { b[i] = PREFIX##_trivial_elem; } \
      spmat_pos_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          demand(eP->row < M->rows, "invalid row index"); \
          demand(eP->col < M->cols, "invalid col index"); \
          ELEM_TYPE *bP = &(b[eP->row]); \
          ELEM_TYPE p = PREFIX##_elem_mul(eP->val, a[eP->col]); \
          ELEM_TYPE s = PREFIX##_elem_add(p, *bP); \
          (*bP) = s; \
        } \
    }

#define spmat_IMPLEMENT_map_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_map_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(M->rows == na, "incompatible row counts"); \
      demand(M->cols == nb, "incompatible column counts"); \
      spmat_index_t j; \
      for (j = 0; j < nb; j++) { b[j] = PREFIX##_trivial_elem; } \
      spmat_pos_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          demand(eP->row < M->rows, "invalid row index"); \
          demand(eP->col < M->cols, "invalid col index"); \
          ELEM_TYPE *bP = &(b[eP->col]); \
          ELEM_TYPE p = PREFIX##_elem_mul(a[eP->row], eP->val); \
          ELEM_TYPE s = PREFIX##_elem_add(p, *bP); \
          (*bP) = s; \
        } \
    }

#define spmat_LINALG_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_identity(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_mix(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_mul(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_map_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_map_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  extern void PREFIX##_linalg_bOgUs /* To eat the semicolon. */

#endif

