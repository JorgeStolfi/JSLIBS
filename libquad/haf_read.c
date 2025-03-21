/* See {haf_read.h}. */
/* Last edited on 2025-01-09 23:42:34 by stolfi */
 
#define haf_read_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright
 
/* Written by J. Stolfi in October 2023. */ 

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <affirm.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>

#include <haf.h>
#include <haf_read.h>

haf_arc_t haf_read_arc(FILE *rd, haf_edge_id_t eid0, haf_edge_vec_t *A)
  { /* Parse the edge id {eid} and get the edge {ed} from the edge table {A}: */
    uint64_t eid = fget_uint64(rd, 10);
    demand(eid >= eid0, "invalid edge id - too small");
    uint64_t ke = eid - eid0;
    demand(ke < A->ne, "invalid edge id - too big");
    haf_edge_t ed = A->e[ke];
    /* Parse ":" and the direction bit {t} in binary: */
    fget_match(rd, ":");
    uint64_t t = fget_uint64(rd, 2);
    demand(t <= 1, "invalid direction bit");
    /* Put it all together: */
    return haf_orient(ed, (haf_dir_bit_t)t);
  }

#define FILE_TYPE "half-edge"
#define FILE_VERSION "2023-10-05"

void haf_read_map(FILE *rd, haf_edge_vec_t *A, haf_edge_id_t *eid0_P, haf_arc_vec_t *R)
  { /* Check and consume the header line: */
    filefmt_read_header(rd, FILE_TYPE, FILE_VERSION);
    
    demand(A != NULL, "edge table {A} is {NULL}");
    demand(R != NULL, "root table {R} is {NULL}");
    
    /* Read the edge count and the min edge id: */
    haf_edge_count_t NE = nget_uint64(rd, "edges", 10); fget_eol(rd);
    demand(NE <= haf_edge_count_MAX, "file has too many edges");
    
    haf_edge_id_t eid0 = nget_uint64(rd, "edge_id_min", 10); fget_eol(rd);
    demand(eid0 <= haf_edge_count_MAX - NE + 1, "file min edge id too large");
    
    /* Create the edge records, save their base arcs in {A->e[0..NE-1]}: */
    haf_edge_vec_expand(A, (vec_index_t)(NE-1));
    for (haf_edge_id_t ke = 0; ke < NE; ke++) 
      { haf_edge_id_t eid = eid0 + ke;
        haf_arc_t a = haf_make_stick(eid); 
        A->e[ke] = haf_edge(a);
      }

    /* Read the contents of the edge records {0..NE-1}: */
    for (haf_edge_count_t ke = 0; ke < NE; ke++) 
      { /* Parse the edge id {eid}: */
        haf_edge_id_t eid_in = fget_uint64(rd, 10);
        demand(eid_in == eid0 + ke, "edge id mismatch");
        /* Get the base arc {ak} from the edge table {A}: */
        haf_edge_t edk = A->e[ke];
        haf_arc_t ak = haf_base_arc(edk);
        /* Read the {.lnext} links of {ak} and {ak.sym}: */
        haf_arc_t b0 = haf_read_arc(rd, eid0, A);
        haf_set_lnext(ak, b0);
        haf_arc_t b1 = haf_read_arc(rd, eid0, A);
        haf_set_lnext(haf_sym(ak), b1);
        /* Skip to the next line: */
        fget_eol(rd);
      }

    /* Parse the root count {NR} and read the root arcs: */
    haf_arc_count_t NR = nget_uint64(rd, "roots", 10); fget_eol(rd);
    demand(NR <= haf_arc_count_MAX, "file has too many roots");
    haf_arc_vec_expand(R, (vec_index_t)(NR-1));
    for (haf_arc_count_t kr = 0; kr < NR; kr++)
      { /* Parse and check the root index {kr}: */
        uint64_t kr_in = fget_uint64(rd, 10);
        demand(kr_in == kr, "root index mismatch");
        /* Parse the root arc number {kr}, save it in {root}: */
        R->e[kr] = haf_read_arc(rd, eid0, A); 
        /* Skip to the next line: */
        fget_eol(rd);
      }
    
    /* Check and consume the footer line: */
    filefmt_read_footer(rd, FILE_TYPE);
    
    /* Trim tables: */
    haf_edge_vec_trim(A, (vec_size_t)NE);
    haf_arc_vec_trim(R, (vec_size_t)NR);
    
    (*eid0_P) = eid0;
  }
