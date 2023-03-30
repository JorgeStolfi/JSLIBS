#ifndef spmat_def_H
#define spmat_def_H
/* The ugly entrails of {spmat.h}. */

#define spmat_def_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 10:46:34 by stolfi */

/* These inclusions are necessary if this file is included or compiled on its own: */
#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <spmat.h>

/* DEFINITION MACROS */

#define spmat_DEFINE_BASIC_TYPES(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  typedef spmat_size_t PREFIX##_size_t; \
  typedef spmat_index_t PREFIX##_index_t; \
  typedef spmat_count_t PREFIX##_count_t; \
  typedef spmat_pos_t PREFIX##_pos_t; \
  typedef ELEM_TYPE PREFIX##_elem_t

#define spmat_DEFINE_ENTRY_TYPE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  typedef struct PREFIX##_entry_t \
    { spmat_index_t row; \
      spmat_index_t col; \
      ELEM_TYPE val; \
    } PREFIX##_entry_t

#define spmat_DEFINE_MATRIX_TYPE(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  typedef struct MATRIX_TYPE \
    { spmat_size_t rows; \
      spmat_size_t cols; \
      PREFIX##_entry_t *e; \
      spmat_count_t ents; \
    } MATRIX_TYPE

#define spmat_DECLARE_new(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  MATRIX_TYPE PREFIX##_new \
    ( spmat_size_t rows, \
      spmat_size_t cols, \
      spmat_count_t ents \
    )

#define spmat_DECLARE_expand(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_expand(MATRIX_TYPE *M, spmat_pos_t pos)

#define spmat_DECLARE_trim(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_trim(MATRIX_TYPE *M, spmat_count_t ents)

#define spmat_DECLARE_make_desc(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  MATRIX_TYPE PREFIX##_make_desc \
    ( spmat_size_t rows, \
      spmat_size_t cols, \
      PREFIX##_entry_t *e, \
      spmat_count_t ents \
    )

#define spmat_DECLARE_copy(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_copy(MATRIX_TYPE *M, MATRIX_TYPE *N)

#define spmat_DECLARE_find_element(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_find_element \
    ( MATRIX_TYPE *M, \
      spmat_pos_t posIni, \
      spmat_pos_t posLim, \
      spmat_index_t row, \
      spmat_index_t col \
    )

#define spmat_DECLARE_add_element(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_add_element \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t row, \
      spmat_index_t col, \
      ELEM_TYPE val \
    )

#define spmat_DECLARE_add_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_add_row \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t row, \
      ELEM_TYPE val[], \
      spmat_size_t nv \
    )

#define spmat_DECLARE_add_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_add_col \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t col, \
      ELEM_TYPE val[], \
      spmat_size_t nv \
    )

#define spmat_DECLARE_add_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_add_diagonal \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t row, \
      spmat_index_t col, \
      ELEM_TYPE val[], \
      spmat_size_t nv \
    )

#define spmat_DECLARE_fill_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_fill_diagonal \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t row, \
      spmat_index_t col, \
      ELEM_TYPE val, \
      spmat_size_t nv \
    )

#define spmat_DECLARE_extract_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_extract_row \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t row, \
      ELEM_TYPE val[], \
      spmat_size_t nv \
    )

#define spmat_DECLARE_extract_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_extract_col \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t col, \
      ELEM_TYPE val[], \
      spmat_size_t nv \
    )

#define spmat_DECLARE_entry_scan_proc_t(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  typedef void PREFIX##_entry_scan_proc_t \
    ( spmat_index_t row, \
      spmat_index_t col, \
      ELEM_TYPE val \
    )

#define spmat_DECLARE_scan_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_scan_row \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t row, \
      PREFIX##_entry_scan_proc_t proc \
    )

#define spmat_DECLARE_scan_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_pos_t PREFIX##_scan_col \
    ( MATRIX_TYPE *M, \
      spmat_pos_t pos, \
      spmat_index_t col, \
      PREFIX##_entry_scan_proc_t proc \
    )

#define spmat_DECLARE_sort_entries(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_sort_entries(MATRIX_TYPE *M, int32_t orow, int32_t ocol)

#define spmat_DECLARE_sort_entries_ins(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_sort_entries_ins \
    ( MATRIX_TYPE *M, \
      int32_t orow, \
      int32_t ocol, \
      spmat_pos_t posIni, \
      spmat_pos_t posLim \
    )

#define spmat_DECLARE_transpose(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_transpose(MATRIX_TYPE *M)

#define spmat_DECLARE_entry_merge_proc_t(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  typedef void PREFIX##_entry_merge_proc_t \
    ( spmat_index_t row, \
      spmat_index_t col, \
      ELEM_TYPE valA, \
      ELEM_TYPE valB \
    )

#define spmat_DECLARE_merge(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_merge(MATRIX_TYPE *A, MATRIX_TYPE *B, PREFIX##_entry_merge_proc_t *proc)

#define spmat_DECLARE_entry_condense_proc_t(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  typedef ELEM_TYPE PREFIX##_entry_condense_proc_t \
    ( spmat_index_t row, \
      spmat_index_t col, \
      ELEM_TYPE v0, \
      ELEM_TYPE v1 \
    )

#define spmat_DECLARE_condense(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_condense(MATRIX_TYPE *A, PREFIX##_entry_condense_proc_t proc)

#define spmat_TYPEDEF(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DEFINE_BASIC_TYPES(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DEFINE_ENTRY_TYPE(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DEFINE_MATRIX_TYPE(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_new(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_expand(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_trim(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_make_desc(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_copy(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_find_element(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_add_element(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_add_row(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_add_col(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_add_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_fill_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_extract_row(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_extract_col(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_entry_scan_proc_t(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_scan_row(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_scan_col(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_sort_entries(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_sort_entries_ins(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_transpose(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_entry_merge_proc_t(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_merge(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_entry_condense_proc_t(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_condense(MATRIX_TYPE,PREFIX,ELEM_TYPE)

/* IMPLEMENTATION MACROS */

#define spmat_IMPLEMENT_new(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_new(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(rows <= spmat_MAX_ROWS, "too many rows"); \
      demand(cols <= spmat_MAX_COLS, "too many cols"); \
      void *e = spmat_alloc(ents, sizeof(PREFIX##_entry_t)); \
      return (MATRIX_TYPE) \
        { .rows = rows, .cols = cols, .e = (PREFIX##_entry_t *)e, .ents = ents }; \
    }

#define spmat_IMPLEMENT_expand(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_expand(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { if (pos >= M->ents) \
        { spmat_expand \
            ( (void**)&(M->e), \
              &(M->ents), \
              pos, \
              sizeof(PREFIX##_entry_t) \
            ); \
        } \
    }

#define spmat_IMPLEMENT_trim(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_trim(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { if (ents != M->ents) \
        { spmat_trim \
            ( (void**)&(M->e), \
              &(M->ents), \
              ents, \
              sizeof(PREFIX##_entry_t) \
            ); \
        } \
    }

#define spmat_IMPLEMENT_make_desc(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_make_desc(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(rows <= spmat_MAX_ROWS, "too many rows"); \
      demand(cols <= spmat_MAX_COLS, "too many cols"); \
      demand(ents <= spmat_MAX_ENTS, "too many entries"); \
      return (MATRIX_TYPE) \
        { .rows = rows, .cols = cols, .e = e, .ents = ents }; \
    }

#define spmat_IMPLEMENT_copy(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_copy(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { N->rows = M->rows; N->cols = M->cols; \
      PREFIX##_trim(N, M->ents); \
      int32_t k; \
      for (k = 0; k < M->ents; k++) { N->e[k] = M->e[k]; } \
    }

#define spmat_BIN_SEARCH_THRESHOLD 5
  /* Min number of entries for which binary search is worth the
    trouble. If we count comparisons only, the number should be 3.
    However, binary search calls {spmat_compare_indices} for
    each test, while sequential search just compares the indices with
    {==}.  Hence we'd better set the threshold higher than 3. */

#define spmat_IMPLEMENT_find_element(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_find_element(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { if (posLim > M->ents) { posLim = M->ents; } \
      if (posIni >= posLim) { return posLim; } \
      spmat_pos_t pos = posIni; \
      PREFIX##_entry_t *eP = &(M->e[pos]); \
      while (pos < posLim) \
        { if ((eP->col == col) && (eP->row == row)) { return pos; } \
          pos++; eP++; \
        } \
      return posLim; \
    }

#define spmat_IMPLEMENT_add_element(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_add_element(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(row <= spmat_MAX_INDEX, "invalid row index"); \
      demand(col <= spmat_MAX_INDEX, "invalid col index"); \
      if (row >= M->rows) { M->rows = row + 1; } \
      if (col >= M->cols) { M->cols = col + 1; } \
      if (PREFIX##_elem_is_trivial(val))                                  \
        { return pos; } \
      else \
        { PREFIX##_expand(M, pos); \
          M->e[pos] = (PREFIX##_entry_t){ .row = row, .col = col, .val = val }; \
          return pos+1; \
        } \
    }

#define spmat_IMPLEMENT_add_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_add_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(row <= spmat_MAX_INDEX, "invalid row index"); \
      demand(nv <= spmat_MAX_COLS, "too many cols"); \
      if (row >= M->rows) { M->rows = row + 1; } \
      if (nv > M->cols) { M->cols = nv; } \
      int32_t k; \
      for (k = 0; k < nv; k++) \
        { ELEM_TYPE *vk = &(val[k]); \
          if (! PREFIX##_elem_is_trivial(*vk)) \
            { PREFIX##_expand(M, pos); \
              M->e[pos] = (PREFIX##_entry_t){ .row = row, .col = k, .val = (*vk) }; \
              pos++; \
            } \
        } \
      return pos; \
    }

#define spmat_IMPLEMENT_add_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_add_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(col <= spmat_MAX_INDEX, "invalid col index"); \
      demand(nv <= spmat_MAX_ROWS, "too many rows"); \
      if (nv > M->rows) { M->rows = nv; } \
      if (col >= M->cols) { M->cols = col + 1; } \
      int32_t k; \
      for (k = 0; k < nv; k++) \
        { ELEM_TYPE *vk = &(val[k]); \
          if (! PREFIX##_elem_is_trivial(*vk)) \
            { PREFIX##_expand(M, pos); \
              M->e[pos] = (PREFIX##_entry_t){ .row = k, .col = col, .val = (*vk) }; \
              pos++; \
            } \
        } \
      return pos; \
    }

#define spmat_IMPLEMENT_add_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_add_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { row = row % M->rows; \
      col = col % M->cols; \
      int32_t k; \
      for (k = 0; k < nv; k++) \
        { ELEM_TYPE *vk = &(val[k]); \
          if (! PREFIX##_elem_is_trivial(*vk)) \
            { PREFIX##_expand(M, pos); \
              M->e[pos] = (PREFIX##_entry_t){ .row = row, .col = col, .val = (*vk) }; \
              pos++; \
            } \
          row++; if (row >= M->rows) { row = 0; } \
          col++; if (col >= M->cols) { col = 0; } \
        } \
      return pos; \
     }

#define spmat_IMPLEMENT_fill_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_fill_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { if (PREFIX##_elem_is_trivial(val)) { return pos; } \
      row = row % M->rows; \
      col = col % M->cols; \
      int32_t k; \
      for (k = 0; k < nv; k++) \
        { PREFIX##_expand(M, pos); \
          M->e[pos] = (PREFIX##_entry_t){ .row = row, .col = col, .val = val }; \
          pos++; \
          row++; if (row >= M->rows) { row = 0; } \
          col++; if (col >= M->cols) { col = 0; } \
        } \
      return pos; \
    }

#define spmat_IMPLEMENT_extract_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_extract_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(nv == M->cols, "incompatible vector size"); \
      int32_t k; \
      for (k = 0; k < nv; k++) { val[k] = PREFIX##_trivial_elem; } \
      while (pos < M->ents) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          if (eP->row != row) { break; } \
          assert(eP->col < M->cols); \
          val[eP->col] = eP->val; \
          pos++; \
        } \
      return pos; \
    }

#define spmat_IMPLEMENT_extract_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_extract_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(nv == M->rows, "incompatible vector size"); \
      int32_t k; \
      for (k = 0; k < nv; k++) { val[k] = PREFIX##_trivial_elem; } \
      while (pos < M->ents) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          if (eP->col != col) { break; } \
          assert(eP->row < M->rows); \
          val[eP->row] = eP->val; \
          pos++; \
        } \
      return pos; \
    }

#define spmat_IMPLEMENT_scan_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_scan_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { while (pos < M->ents) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          if (eP->row != row) { break; } \
          assert(eP->col < M->cols); \
          if (! PREFIX##_elem_is_trivial(eP->val)) \
            { proc(row, eP->col, eP->val); } \
          pos++; \
        } \
      return pos; \
    }

#define spmat_IMPLEMENT_scan_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_scan_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { while (pos < M->ents) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          if (eP->col != col) { break; } \
          assert(eP->row < M->rows); \
          if (! PREFIX##_elem_is_trivial(eP->val)) \
            { proc(eP->row, col, eP->val); } \
          pos++; \
        } \
      return pos; \
    }

#define spmat_IMPLEMENT_sort_entries(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_sort_entries(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { if ((orow == 0) && (ocol == 0)) { return; } \
      auto int32_t cmp(const void *avP, const void *bvP); \
      int32_t cmp(const void *avP, const void *bvP) \
        { const PREFIX##_entry_t *aP = avP; \
          const PREFIX##_entry_t *bP = bvP; \
          return spmat_compare_indices \
            (aP->row, aP->col, bP->row, bP->col, orow, ocol); \
        } \
      qsort(M->e, M->ents, sizeof(PREFIX##_entry_t), &cmp); \
    }

#define spmat_IMPLEMENT_sort_entries_ins(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_sort_entries_ins(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { spmat_pos_t j; \
      for (j = posIni + 1; j < posLim; j++) \
        { spmat_pos_t k = j, p = k - 1; \
          PREFIX##_entry_t *ek = &(M->e[k]); \
          PREFIX##_entry_t *ep = &(M->e[p]); \
          int32_t cmp = spmat_compare_indices \
            (ep->row, ep->col, ek->row, ek->col, orow, ocol); \
          if (cmp > 0) \
            { /* Entry {*ek} is out of order, insert it in its proper place: */ \
              PREFIX##_entry_t etmp = (*ek); \
              do \
                { (*ek) = (*ep); \
                  k = p; ek = ep; \
                  if (k <= posIni) { break; } \
                  p--; ep--; \
                  cmp = spmat_compare_indices \
                    (ep->row, ep->col, etmp.row, etmp.col, orow, ocol); \
                } \
              while (cmp > 0); \
              (*ek) = etmp; \
            } \
        } \
    }

#define spmat_IMPLEMENT_transpose(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_transpose(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { spmat_pos_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          spmat_index_t t = eP->row; eP->row = eP->col; eP->col = t; \
        } \
    }

#define spmat_IMPLEMENT_merge(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_merge(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { demand(A->rows == B->rows, "incompatible row counts"); \
      demand(A->cols == B->cols, "incompatible column counts"); \
      spmat_pos_t posA = 0, posB = 0; \
      spmat_index_t rowPrev = spmat_NO_INDEX; \
      spmat_index_t colPrev = spmat_NO_INDEX; \
      while ((posA < A->ents) || (posB < B->ents)) \
        { PREFIX##_entry_t *aP = (posA < A->ents ? &(A->e[posA]) : NULL); \
          PREFIX##_entry_t *bP = (posB < B->ents ? &(B->e[posB]) : NULL); \
          int32_t cmp; /* Which entry comes first? */ \
          if (aP == NULL) \
            { cmp = +1; } \
          else if (bP == NULL) \
            { cmp = -1; } \
          else \
            { cmp = spmat_compare_indices \
                (aP->row, aP->col, bP->row, bP->col, +2, +1); } \
          spmat_index_t row, col; \
          ELEM_TYPE Aij, Bij; \
          if (cmp < 0) \
            { row = aP->row; col = aP->col; \
              Aij = aP->val; posA++; \
              Bij = PREFIX##_trivial_elem; \
            } \
          else if (cmp > 0) \
            { row = bP->row; col = bP->col; \
              Aij = PREFIX##_trivial_elem; \
              Bij = bP->val; posB++; \
            } \
          else \
            { assert(aP->row == bP->row); \
              assert(aP->col == bP->col); \
              row = aP->row; col = aP->col; \
              Aij = aP->val; posA++; \
              Bij = bP->val; posB++; \
            } \
          if (rowPrev != spmat_NO_INDEX) \
            { demand \
                ( (row > rowPrev) || \
                  ((row == rowPrev) && (col >= colPrev)), \
                  "entries out of order" \
                ); \
            } \
          proc(row, col, Aij, Bij); \
          rowPrev = row; colPrev = col; \
        } \
    }

#define spmat_IMPLEMENT_condense(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_condense(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { if (A->ents == 0) { return; } \
      spmat_pos_t posC = 0; PREFIX##_entry_t *cP = A->e; \
      spmat_pos_t posA = 1; PREFIX##_entry_t *aP = cP + 1; \
      while (posA < A->ents) \
        { if ((cP->col == aP->col) && (cP->row == aP->row)) \
            { cP->val = proc(cP->row, cP->col, cP->val, aP->val); } \
          else \
            { if (! PREFIX##_elem_is_trivial(cP->val)) { posC++; cP++; } \
              cP->row = aP->row; cP->col = aP->col; \
              cP->val = aP->val; \
            } \
          posA++; aP++; \
        } \
      if (! PREFIX##_elem_is_trivial(cP->val)) { posC++; cP++; } \
      if (posC < A->ents) { PREFIX##_trim(A, posC); } \
    }

#define spmat_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_new(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_expand(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_trim(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_make_desc(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_copy(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_find_element(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_add_element(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_add_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_add_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_add_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_fill_diagonal(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_extract_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_extract_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_scan_row(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_scan_col(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_sort_entries(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_sort_entries_ins(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_transpose(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_merge(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_condense(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  extern void PREFIX##_bOgUs /* To eat the semicolon. */

#endif
