/* See {array_io.h} */

#define array_io_C_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2017-01-02 21:40:42 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <affirm.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <jsmath.h>
#include <jswsize.h>

#include <array.h>
#include <array_io.h>

#define array_FILE_VERSION "2009-08-31"

void array_write_header
  ( FILE *wr, 
    char *type, 
    char *cmt, 
    ix_descr_t *D,
    ix_order_t ixor,
    int dix[]
  )
  {
    ix_dim_t na = D->na;

    /* Write header line: */
    filefmt_write_header(wr, type, array_FILE_VERSION);

    /* Write comment lines, if any: */
    if (cmt != NULL) { filefmt_write_comment(wr, cmt, '#'); }

    /* Write the effective dimension {na} and the respective sizes: */
    int ia;
    fprintf(wr, "axes = %d\n", na);
    fprintf(wr, "size =");
    for (ia = 0; ia < na; ia++)
      { ix_size_t szi = D->sz[ia];
        fprintf(wr, (" %" uint64_u_fmt), szi);
        dix[ia] = (szi == 0 ? 0 : digits(szi-1)); /* Number of digits of each index. */
      }
    fprintf(wr, "\n");
    fprintf(wr, "order = %d\n", ixor);
  }

void array_write_footer(FILE *wr, char *type)
  {
    filefmt_write_footer(wr, type);
    fflush(wr);
  }

void array_read_header
  ( FILE *rd, 
    char *type, 
    char **cmtP, 
    ix_descr_t *D,
    ix_order_t *ixorP
  )
  {
    /* Read header line: */
    filefmt_read_header(rd, type, array_FILE_VERSION);

    /* Read the comment lines, if any: */
    char *cmt = filefmt_read_comment(rd, '#');
    if (cmtP == NULL)
      { free(cmt); }
    else
      { (*cmtP) = cmt; }

    /* Read the effective number of axes {na}: */
    uint64_t d64 = nget_uint64(rd, "axes", 10); fget_eol(rd);
    demand(d64 <= array_MAX_AXES, "array has too many axes");
    ix_dim_t na = (ix_dim_t)d64;
    
    /* Read the size vector: */
    ix_size_t sz[na];
    nget_name_eq(rd, "size");
    int i;
    for (i = 0; i < na; i++)
      { uint64_t sz64 = (i < na ? fget_uint64(rd, 10) : 1);
        demand(sz64 <= ix_MAX_SIZE, "size too large");
        sz[i] = sz64;
      }
    fget_eol(rd);

    /* Read the index scan order: */
    uint64_t ixor64 = nget_uint64(rd, "order", 10); fget_eol(rd);
    demand(ixor64 <= 1, "invalid element order");
    ix_order_t ixor = ixor64;
    
    /* Compute indexing increments and return to caller: */
    (*D) = ix_descr_from_sizes(na, sz, ixor);
    (*ixorP) = ixor;
  }

void array_read_footer(FILE *rd, char *type)
  {
    filefmt_read_footer(rd, type);
  }

