/* See {spmat_io.h} */

#define spmat_io_C_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2019-04-09 13:01:49 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <affirm.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>

#include <spmat.h>
#include <spmat_io.h>

#define spmat_FILE_VERSION "2008-07-24"

void spmat_read_header
  ( FILE *rd, 
    char *type, 
    char **cmtP, 
    spmat_size_t *rowsP,
    spmat_size_t *colsP, 
    spmat_count_t *entsP
  )
  {
    /* Read header line: */
    filefmt_read_header(rd, type, spmat_FILE_VERSION);

    /* Read the comment lines, if any: */
    char *cmt = filefmt_read_comment(rd, '#');
    if (cmtP == NULL)
      { free(cmt); }
    else
      { (*cmtP) = cmt; }

    /* Read the matrix dimensions: */
    (*rowsP) = nget_uint32(rd, "rows", 10); fget_eol(rd);
    demand((*rowsP) <= spmat_MAX_ROWS, "too many rows");
    (*colsP) = nget_uint32(rd, "columns", 10); fget_eol(rd);
    demand((*colsP) <= spmat_MAX_COLS, "too many cols");
    (*entsP) = nget_uint32(rd, "entries", 10); fget_eol(rd);
    demand((*entsP) <= spmat_MAX_ENTS, "too many entries");
  }

void spmat_read_footer(FILE *rd, char *type)
  {
    filefmt_read_footer(rd, type);
  }

void spmat_write_header
  ( FILE *wr, 
    char *type, 
    char *cmt, 
    spmat_size_t rows, 
    spmat_size_t cols, 
    spmat_count_t ents
  )
  {
    /* Write header line: */
    filefmt_write_header(wr, type, spmat_FILE_VERSION);

    /* Write comment lines, if any: */
    int ind = 0; /* Comment indentation. */
    if (cmt != NULL) { filefmt_write_comment(wr, cmt, ind, '#'); }

    /* Write the matrix dimensions: */
    fprintf(wr, "rows = %u\n", rows);
    fprintf(wr, "columns = %u\n", cols);
    fprintf(wr, "entries = %u\n", ents);
  }

void spmat_write_footer(FILE *wr, char *type)
  {
    filefmt_write_footer(wr, type);
    fflush(wr);
  }

