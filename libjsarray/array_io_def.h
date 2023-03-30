#ifndef array_io_def_H
#define array_io_def_H
/* The ugly entrails of {array_io.h}. */

#define array_io_def_H_COPYRIGHT "Copyright © 2009 by J. Stolfi, UNICAMP"
/* Created on 2009-08-31 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 11:03:16 by stolfi */

/* These inclusions are necessary if this file is included or compiled on its own: */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>

#include <array.h>
#include <array_io.h>

void array_write_header
  ( FILE *wr,
    char *type,
    char *cmt,
    ix_descr_t *D,
    int32_t dix[]
  );

void array_write_footer(FILE *wr, char *type);

void array_read_header
  ( FILE *rd,
    char *type,
    char **cmtP,
    ix_descr_t *D
  );

void array_read_footer(FILE *rd, char *type);

/* DEFINITION MACROS */

#define array_DECLARE_elem_write(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_elem_write(FILE *wr, ELEM_TYPE *valP)

#define array_DECLARE_elem_read(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_elem_read(FILE *rd, ELEM_TYPE *valP)

#define array_DECLARE_write(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_write(FILE *wr, ARRAY_TYPE *A)

#define array_DECLARE_read(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  ARRAY_TYPE *PREFIX##_read(FILE *rd)

#define array_io_DEF(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_elem_write(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_elem_read(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_write(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_read(ARRAY_TYPE,PREFIX,ELEM_TYPE)

/* IMPLEMENTATION MACROS */

#define array_IMPLEMENT_write(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_write(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ix_descr_t *DA = &(A->ds); \
      ix_dim_t na = DA->na; \
      int32_t ia; \
      int32_t dix[na]; /* Number of columns to use for each index. */ \
      array_write_header(wr, #ARRAY_TYPE, NULL, DA, dix); \
      ix_index_t ix[na]; \
      bool_t first = TRUE; \
      if (ix_descr_indices_first(DA, ix)) \
        { ix_pos_t pA = ix_descr_position(DA, ix); \
          do { \
            if (! first) \
              { /* Print one blank line for each trailing 0 in {ix[0..na-1]}: */ \
                ia = ix_sync_level(na,DA->sz,ix,ix_order_L,FALSE);            \
                while (ia > 0) { fprintf(wr, "\n"); ia--; } \
              } \
            /* Print the index tuple: */ \
            for (ia = 0; ia < na; ia++) { fprintf(wr, " %*ld", dix[ia], ix[ia]); } \
            /* Print the value: */ \
            PREFIX##_elem_write(wr, &(A->e[pA])); \
            fprintf(wr, "\n"); \
            first = FALSE; \
          } while (! ix_descr_next(DA,ix,&pA)); \
        } \
      array_write_footer(wr, #ARRAY_TYPE); \
    }

#define array_IMPLEMENT_read(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_read(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { /* Allocate and read the array descriptor: */ \
      ARRAY_TYPE *A = PREFIX##_new_descr(); \
      ix_descr_t *DA = &(A->ds); \
      array_read_header(rd, #ARRAY_TYPE, NULL, DA); \
      ix_dim_t na = DA->na; \
      int32_t ia; \
      /* Allocate the element store: */ \
      ix_count_t ne = ix_descr_num_positions(DA); \
      if (ne == 0) \
        { A->e = (ELEM_TYPE *)NULL; } \
      else \
        { A->e = (ELEM_TYPE *)notnull(malloc(ne*sizeof(ELEM_TYPE)), "no mem"); } \
      /* Read the elements, in the specified order: */ \
      ix_index_t ix[na]; \
      uint64_t ix64[na]; \
      if (ix_descr_indices_first(DA, ix)) \
        { ix_pos_t p = ix_descr_position(DA, ix); \
          do { \
            /* Skip blank lines, if any: */ \
            fget_skip_formatting_chars(rd); \
            /* Read and check the index tuple: */ \
            for (ia = 0; ia < na; ia++) \
              { ix64[ia] = fget_uint64(rd, 10); \
                demand(ix64[ia] == ix[ia], "indices invalid or out of order"); \
              } \
            /* Read the value: */ \
            PREFIX##_elem_read(rd, &(A->e[p])); \
            fget_eol(rd); \
          } while (! ix_descr_next(DA,ix,&p)); \
        } \
      array_read_footer(rd, #ARRAY_TYPE); \
      return A; \
    }

#define array_io_IMPL(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_write(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_read(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  extern void PREFIX##_io_bOgUs /* To eat the semicolon. */

#endif
