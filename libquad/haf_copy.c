/* See {haf_copy.h}. */
/* Last edited on 2025-01-09 23:33:40 by stolfi */
 
#define haf_copy_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <affirm.h>

#include <haf.h>
#include <haf_enum.h>

#include <haf_copy.h>

/* IMPLEMENTATIONS */

void haf_copy
  ( haf_arc_count_t NR,
    haf_arc_t root[],
    haf_edge_id_t eid0,
    haf_edge_vec_t *new_P,
    haf_edge_vec_t *old_P
  )
  { 
    /* Provide internal edge tables if {new_P}, {old_P} are is {NULL}: */
    haf_edge_vec_t new_local = haf_edge_vec_new(0);
    haf_edge_vec_t *new = (new_P != NULL ? new_P : &new_local);

    haf_edge_vec_t old_local = haf_edge_vec_new(0);
    haf_edge_vec_t *old = (old_P != NULL ? old_P : &old_local);

    /* Internal edge renumnering table: */
    haf_edge_id_vec_t oid = haf_edge_id_vec_new(0);
    
    /* Enumerate reachable edges in the old copy, renumbering them: */
    haf_enum_edges(NR, root, eid0, old, &oid, NULL);
    uint32_t NE = old->ne;
    
    if (NE > 0)
      { /* Create copies of those edges: */
        haf_edge_vec_expand(new, (vec_index_t)NE - 1);
        for (uint64_t ke = 0; ke < NE; ke++)
          { new->e[ke] = haf_edge(haf_make_stick(eid0 + ke)); }
        haf_edge_vec_trim(new, NE);

        /* Stitch the new arcs like the old ones: */
        for (int64_t ke = 0; ke < NE; ke++)
          { haf_edge_t nek = new->e[ke];
            haf_edge_t oek = old->e[ke];
            for (int32_t dbit = 0; dbit <= 1; dbit++)
              { haf_arc_t nak = haf_orient(nek, (haf_dir_bit_t)dbit);
                haf_arc_t oak = haf_orient(oek, (haf_dir_bit_t)dbit);
                haf_arc_t oak_lprev = haf_lprev(oak);
                /* Determine the arc {nak_lprev} that should become {lprev(nak)}: */
                haf_edge_id_t oak_lprev_eid = haf_edge_id(haf_edge(oak_lprev));
                haf_dir_bit_t oak_lprev_bit = haf_dir_bit(oak_lprev);
                haf_edge_t nek_lprev = new->e[oak_lprev_eid];
                haf_arc_t nak_lprev = haf_orient(nek_lprev, oak_lprev_bit);
                assert(haf_lnext(nak) == haf_sym(nak)); /* End is still loose. */
                haf_splice(nak, nak_lprev);
              }
          }
        /* Restore the original edge IDs of the original edges: */
        for (int64_t ke = 0; ke < NE; ke++)
          { haf_edge_t oek = old->e[ke];
            haf_set_edge_id(oek, oid.e[ke]);
          }
      }
      
      
    /* Trim or recycle the {new} table: */
    if (new == &new_local) 
      { assert(new_P == NULL); free(new->e); } 
    else 
      { assert(new == new_P);  assert(new_local.e == NULL);
        haf_edge_vec_trim(new, (vec_size_t)NE);
      }
      
    /* Trim or recycle the {old} table: */
    if (old == &old_local) 
      { assert(old_P == NULL); free(old->e); } 
    else 
      { assert(old == old_P);  assert(old_local.e == NULL);
        haf_edge_vec_trim(old, (vec_size_t)NE);
      }
       
    /* Recycle the old ids table: */
    free(oid.e);
    
    return;
  }
