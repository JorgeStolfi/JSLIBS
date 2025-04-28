/* See {psgr_copy.h} */
/* Last edited on 2025-04-27 17:41:16 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <i2.h>

#include <psgr_types.h>
#include <psgr.h>

#include <psgr_copy.h>

/* SHORTER NAMES */  

#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE

#define CLEAR   psgr_mark_CLEAR
#define KILL    psgr_mark_KILL   
#define KEEP    psgr_mark_KEEP   
#define DELETED psgr_mark_DELETED

psgr_t *psgr_copy(psgr_t *gri, int32_t vj_from_vi[], int32_t ej_from_ei[])
  { uint32_t NVi = psgr_num_vertices(gri);
    uint32_t NEi = psgr_num_edges(gri);
    
    /* Count non-deleted vertices and fill {vj_from_vi}: */
    uint32_t NVj = 0;
    for (uint32_t kvi = 0; kvi < NVi; kvi++)
      { psgr_vertex_t vi = (psgr_vertex_t){ kvi };
        psgr_mark_t mkvi = psgr_vertex_mark(gri, vi);
        if (mkvi != DELETED)
          { if (vj_from_vi != NULL) { vj_from_vi[kvi] = (int32_t)NVj; }
            NVj++;
          }
        else
          { demand(vj_from_vi != NULL, "table {vi_from_vj} is required"); 
            vj_from_vi[kvi] = -1;
          }
      }
    assert(NVj <= NVi);
    assert((NVj == NVi) || (vj_from_vi != NULL));
      
    /* Count non-deleted edges and fill {ej_from_ei}: */
    uint32_t NEj = 0;
    for (uint32_t kei = 0; kei < NEi; kei++)
      { psgr_edge_t ei = (psgr_edge_t){ kei };
        psgr_mark_t mkei = psgr_edge_mark(gri, ei);
        if (mkei != DELETED)
          { if (ej_from_ei != NULL) { ej_from_ei[kei] = (int32_t)NEj; }
            NEj++;
          }
        else
          { demand(ej_from_ei != NULL, "table {vi_from_vj} is required"); 
            ej_from_ei[kei] = -1;
          }
      }
    assert(NEj <= NEi);
     
    /* Create the new graph: */
    i2_t sz =  psgr_image_size(gri);
    uint32_t NX = (uint32_t)sz.c[0], NY = (uint32_t)sz.c[1];
    psgr_t  *grj = psgr_new(NVj, NEj, NX, NY);
    
    /* Copy the non-deleted vertices: */
    for (uint32_t kvi = 0; kvi < NVi; kvi++)
      { psgr_vertex_t vi = (psgr_vertex_t){ kvi };
        psgr_mark_t mkvi = psgr_vertex_mark(gri, vi);
        if (mkvi != DELETED)
          { i2_t ipi = psgr_vertex_pixel(gri, vi);
            psgr_vertex_t vj = psgr_add_vertex(grj, &ipi);
            assert(vj.ix == psgr_num_vertices(grj) - 1);
            if (vj_from_vi != NULL) { assert(vj.ix == vj_from_vi[kvi]); }
            psgr_vertex_set_mark(grj, vj, mkvi);
          }
      }
    assert(psgr_num_vertices(grj) == NVj);
      
    /* Copy the non-deleted edges: */
    for (uint32_t kei = 0; kei < NEi; kei++)
      { psgr_edge_t ei = (psgr_edge_t){ kei };
        psgr_arc_t ai = psgr_orient_edge(ei, 0);
        psgr_mark_t mkei = psgr_edge_mark(gri, ei);
        if (mkei != DELETED)
          { psgr_vertex_t orgi = psgr_arc_org(gri, ai);
            demand(orgi.ix < NVi, "edge with invalid {org}");
            psgr_vertex_t dsti = psgr_arc_dst(gri, ai);
            demand(dsti.ix < NVi, "edge with invalid {dst}");
            demand(orgi.ix != dsti.ix, "loop edge");
            
            psgr_path_t Pj = psgr_path_copy(psgr_arc_path(gri, ai));
            
            double delta = psgr_arc_delta(gri, ai);
            double weight = psgr_arc_weight(gri, ai);
            demand(isfinite(delta), "edge delta is not finite");
            demand(isfinite(weight) && (weight >= 0), "edge weight is invalid");
            demand(weight > 0, "edge with zero weight");

            /*  Check whether org and dst are mapped correctly: */
            psgr_vertex_t orgj, dstj;
            if (vj_from_vi == NULL)
              { /* There must have been no deletions: */
                assert(NVj == NVi);
                orgj = orgi; dstj = dsti;
              }
            else
              { demand(vj_from_vi[orgi.ix] >= 0, "non-deleted edge with deleted {org}"); 
                orgj = (psgr_vertex_t){ (uint32_t)vj_from_vi[orgi.ix] }; 
                demand(vj_from_vi[dsti.ix] >= 0, "non-deleted edge with deleted {dst}");
                dstj = (psgr_vertex_t){ (uint32_t)vj_from_vi[dsti.ix] };
              }
            psgr_arc_t aj = psgr_add_edge(grj, orgj,dstj, delta,weight, Pj);
            psgr_edge_t ej = psgr_arc_edge(aj);
            assert(ej.ix == psgr_num_edges(grj) - 1);
            if (ej_from_ei != NULL) { assert(ej.ix == ej_from_ei[kei]); }
          }
      }
    assert(psgr_num_edges(grj) == NEj);
   
    return grj;
  }
