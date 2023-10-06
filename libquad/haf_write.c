/* See {haf_write.h}. */
/* Last edited on 2023-10-05 20:26:01 by stolfi */
 
#define haf_write_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <affirm.h>
#include <jswsize.h>
#include <jsmath.h>
#include <filefmt.h>

#include <haf.h>
#include <haf_write.h>

void haf_write_arc(FILE *wr, haf_arc_t a, int32_t width)
  { char *fmt = ("%*" uint64_u_fmt ":%u");
    fprintf(wr, fmt, width, haf_edge_id(a), (uint32_t)haf_dir_bit(a));
  }

void haf_write_map
  ( FILE *wr,
    haf_edge_count_t ne,
    haf_arc_t a[],
    haf_edge_id_t eid0,
    haf_arc_count_t nr,
    haf_arc_t root[]
  )
  { 
    /* Write the header line: */
    filefmt_write_header(wr, haf_write_FILE_TYPE, haf_write_FILE_VERSION);
    /* Grab the number {nr} of roots and the root index width {wa}: */
    /* Compute the width in digits {wr} of the root index: */
    int32_t dr = (nr < 10 ? 1 : digits(nr-1));
    
    /* We should have zero edges if and only if we have zero roots: */
    demand((nr == 0) == (ne == 0), "need at least one root"); 
    
    /* Determine the width {eE} of edge ids in digits of the edge id: */
    int32_t dE = (ne < 10 ? 1 : digits(ne-1));

    /* Write the edge table, one edge per line: */
    fprintf(wr, "edges = %lu\n", ne);
    fprintf(wr, "edge_id_min = %lu\n", eid0);
    for (haf_edge_count_t ke = 0; ke < ne; ke++)
      { haf_arc_t ak = a[ke];
        assert(ak != NULL);
        haf_edge_id_t eid = haf_edge_id(ak);
        /* Check whether the numbering is as specified: */
        demand(eid == eid0 + ke, "edge id inconsistent with index");
        demand(haf_dir_bit(ak) == 0, "given arc is not the base arc");
        /* Write the edge's number and its {lnext} links: */
        fprintf(wr, ("%*" uint64_u_fmt), dE, eid);
        fputc(' ', wr); 
        haf_write_arc(wr, haf_lnext(ak), dE);
        fputc(' ', wr); 
        haf_write_arc(wr, haf_lnext(haf_sym(ak)), dE);
        fputc('\n', wr);
      }

    /* Write the root arcs: */
    fprintf(wr, "roots = %lu\n", nr);
    /* Write the roots, one per line: */
    for (uint32_t kr = 0; kr < nr; kr++)
      { /* Write {kr} and the root arc number {kr}: */
        fprintf(wr, "%*u ", dr, kr);
        haf_write_arc(wr, root[kr], dE);
        fputc('\n', wr);
      }
    /* Write the footer line: */
    filefmt_write_footer(wr, haf_write_FILE_TYPE);
    fflush(wr);
    return;
    
  }
