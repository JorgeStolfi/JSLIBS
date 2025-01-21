/* See {haf_write.h}. */
/* Last edited on 2025-01-09 23:57:31 by stolfi */
 
#define haf_write_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

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

void haf_write_arc(FILE *wr, haf_arc_t a, uint32_t width)
  { char *fmt = ("%*" uint64_u_fmt ":%u");
    fprintf(wr, fmt, width, haf_edge_id(haf_edge(a)), (uint32_t)haf_dir_bit(a));
  }

void haf_write_map
  ( FILE *wr,
    haf_edge_count_t NE,
    haf_edge_t ed[],
    haf_edge_id_t eid0,
    haf_arc_count_t NR,
    haf_arc_t root[]
  )
  { 
    /* Write the header line: */
    filefmt_write_header(wr, haf_write_FILE_TYPE, haf_write_FILE_VERSION);
    /* Grab the number {NR} of roots and the root index width {wa}: */
    /* Compute the width in digits {wr} of the root index: */
    uint32_t dr = (NR < 10 ? 1 : digits(NR-1));
    
    /* Determine the width {eE} of edge ids in digits of the edge id: */
    uint32_t dE = (NE < 10 ? 1 : digits(NE-1));

    /* Write the edge table, one edge per line: */
    fprintf(wr, "edges = %lu\n", NE);
    fprintf(wr, "edge_id_min = %lu\n", eid0);
    for (haf_edge_count_t ke = 0; ke < NE; ke++)
      { haf_edge_t edk = ed[ke];
        assert(edk != NULL);
        haf_edge_id_t eid = haf_edge_id(edk);
        /* Check whether the numbering is as specified: */
        demand(eid == eid0 + ke, "edge id inconsistent with index");
        /* Write the edge's number and its {lnext} links: */
        fprintf(wr, ("%*" uint64_u_fmt), dE, eid);
        fputc(' ', wr); 
        haf_arc_t ak = haf_base_arc(edk);
        haf_write_arc(wr, haf_lnext(ak), dE);
        fputc(' ', wr); 
        haf_write_arc(wr, haf_lnext(haf_sym(ak)), dE);
        fputc('\n', wr);
      }

    /* Write the root arcs: */
    fprintf(wr, "roots = %lu\n", NR);
    /* Write the roots, one per line: */
    for (uint32_t kr = 0;  kr < NR; kr++)
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
