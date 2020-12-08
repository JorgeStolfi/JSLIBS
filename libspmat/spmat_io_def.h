#ifndef spmat_io_def_H
#define spmat_io_def_H
/* The ugly entrails of {spmat_io.h}. */

#define spmat_io_def_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2009-08-31 21:50:09 by stolfi */

/* These inclusions are necessary if this file is included or compiled on its own: */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>

#include <spmat.h>
#include <spmat_io.h>

void spmat_read_header
  ( FILE *rd,
    char *type,
    char **cmtP,
    spmat_size_t *rowsP,
    spmat_size_t *colsP,
    spmat_count_t *entsP
  );

void spmat_read_footer(FILE *rd, char *type);

void spmat_write_header
  ( FILE *wr,
    char *type,
    char *cmt,
    spmat_size_t rows,
    spmat_size_t cols,
    spmat_count_t ents
  );

void spmat_write_footer(FILE *wr, char *type);

/* DEFINITION MACROS */

#define spmat_DECLARE_elem_write(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_elem_write(FILE *wr, ELEM_TYPE *valP)

#define spmat_DECLARE_elem_read(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_elem_read(FILE *rd, ELEM_TYPE *valP)

#define spmat_DECLARE_write(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_write(FILE *wr, MATRIX_TYPE *M)

#define spmat_DECLARE_read(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_read(FILE *rd, MATRIX_TYPE *M)

#define spmat_io_DEF(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_elem_write(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_elem_read(MATRIX_TYPE,PREFIX,ELEM_TYPE);  \
  spmat_DECLARE_write(MATRIX_TYPE,PREFIX,ELEM_TYPE); \
  spmat_DECLARE_read(MATRIX_TYPE,PREFIX,ELEM_TYPE)

/* IMPLEMENTATION MACROS */

#define spmat_IMPLEMENT_write(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_write(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { spmat_write_header(wr, #MATRIX_TYPE, NULL, M->rows, M->cols, M->ents); \
      PREFIX##_pos_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          demand(eP->row < M->rows, "invalid row index"); \
          demand(eP->col < M->cols, "invalid col index"); \
          fprintf(wr, "%7d %7d ", eP->row, eP->col); \
          PREFIX##_elem_write(wr, &(eP->val)); \
          fprintf(wr, "\n"); \
        } \
      spmat_write_footer(wr, #MATRIX_TYPE); \
    }

#define spmat_IMPLEMENT_read(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_DECLARE_read(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
    { PREFIX##_count_t ents; \
      spmat_read_header(rd, #MATRIX_TYPE, NULL, &(M->rows), &(M->cols), &ents); \
      PREFIX##_trim(M, ents); \
      PREFIX##_pos_t pos; \
      for (pos = 0; pos < M->ents; pos++) \
        { PREFIX##_entry_t *eP = &(M->e[pos]); \
          eP->row = fget_uint32(rd, 10); \
          demand(eP->row < M->rows, "invalid row index"); \
          eP->col = fget_uint32(rd, 10); \
          demand(eP->col < M->cols, "invalid col index"); \
          fget_skip_spaces(rd); \
          PREFIX##_elem_read(rd, &(eP->val)); \
          fget_eol(rd); \
        } \
      spmat_read_footer(rd, #MATRIX_TYPE); \
    }

#define spmat_io_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_write(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_IMPLEMENT_read(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  extern void PREFIX##_io_bOgUs /* To eat the semicolon. */

#endif
