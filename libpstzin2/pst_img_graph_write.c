/* See {pst_img_graph_write.h} */
/* Last edited on 2025-01-10 07:49:32 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <haf.h>
#include <haf_write.h>
#include <filefmt.h>
#include <affirm.h>

#include <pst_img_graph.h>

#include <pst_img_graph_write.h>

void pst_img_graph_write_vertex(FILE *wr, pst_img_graph_t *g, uint32_t vi);
  /* Writes the vertex {g->vertex[vi]} to {wr} in a format compatible
    with {pst_img_graph_read_vertex}. */

void pst_img_graph_write_edge(FILE* wr,pst_img_graph_t* g, haf_arc_t e);

/* IMPLEMENTATIONS */

#define pst_img_graph_FILE_TYPE "pst_img_graph_t"
#define pst_img_graph_FILE_VERSION "2025-01-09"

void pst_img_graph_write(FILE *wr,pst_img_graph_t *g)
  {
    filefmt_write_header(wr, pst_img_graph_FILE_TYPE, pst_img_graph_FILE_VERSION);
    fprintf(wr, "NV = %d\n", g->NV);
    fprintf(wr, "NE = %d\n", g->NE);
    fprintf(wr,"\n");
    fprintf(wr, "topology\n");
    uint32_t eid0 = 0; /* Lowest edge ID. */
    uint32_t nroots = 0;
    haf_write_map(wr, g->NE, g->arc, eid0, nroots, NULL);
    fprintf(wr,"\n");
    fprintf(wr, "vertices\n");
    for (int32_t i = 0; i < g->NV; i++)
      { pst_img_graph_write_vertex(wr, g, i); }
    fprintf(wr,"\n");
    fprintf(wr,"edges\n");
    for (int32_t i = 0; i < g->NE; i++)
      { pst_img_graph_write_edge(wr, g, i, &(g->{hedge|dedge}@@[i], FALSE)); }
    fprintf(wr,"\n");
    filefmt_write_footer(wr, pst_img_graph_FILE_TYPE);
  }

void pst_img_graph_write_vertex(FILE *wr, pst_img_graph_t *g, uint32_t vi)
  { demand(vi <= g->NV, "invalid vertex index {vi}");
    fprintf(wr, "%9d", vi); 
    pst_img_graph_vertex_data_t *vd = &(g->vertex[vi]);
    if (vd == NULL)
      { fprintf(wr, " NULL\n"); }
    else
      { fprintf(wr, " [%d,%d] at (%f,%f) out", vd->x, vd->y, vd->coords.c[0], vd->coords.c[1]);
        if (vd->aout == NULL)
          { fprintf(wr, " NULL"); }
        else
          { uint32_t eout_id = (uint32_t)haf_edge_id(haf_edge(vd->aout));
            uint32_t aout_bit = (uint32_t)haf_dir_bit(vd->aout);
            fprintf(wr, " %d:%d", eout_id, aout_bit);
          }
        fprintf(wr, "\n");
      }
  }

void pst_img_graph_write_edge(FILE *wr, pst_img_graph_t *g, uint32_t eid)
  { demand(eid < g->NE, "invalid edge index {eid}");
    fprintf(wr, "%9d", eid); 
    pst_img_graph_edge_data_t *ed = g->{hedge|dedge}@@[eid].data;
    haf_arc_t *a = g->{hedge|dedge}@@[eid].arc;
    assert((ed == NULL) == (a == NULL));
    if (a == NULL) 
      { fprintf(wr, " NULL"); }
    else
      { assert(eid == (uint32_t)haf_edge_id(haf_edge(a)));
        uint32_t abit = (uint32_t)haf_dir_bit_id(a);
        int32_t org = pst_img_graph_get_arc_origin(g, a);
        int32_t dst = pst_img_graph_get_arc_origin(g, haf_sym(a));
        fprintf(wr," %d:%d (v%d --> v%d)", eid, abit, org, dst);
        
        char *label = (g->{hedge|dedge}@@[eid].data->label == NULL ? "": g->{hedge|dedge}@@[eid].data->label);
        fprintf(wr," (%d)'%s'", strlen(label), label);
        
        double   w  = pst_img_graph_get_edge_weight(g,a);
        double   d  = pst_img_graph_get_arc_delta(g,a);
        assert(pst_img_graph_get_arc_delta(g,haf_sym(a) == -d);
        fprintf(wr," d = %+.8f w = %.6f", d, w);
        
        pst_path_t p = pst_img_graph_get_edge_path(g,a);
        if (p.n > 0)
          { fprintf(wr," path = ");
            ps_img_graph_write_path(wr, p);
          }
      }
    fprintf(wr,"\n");
  }
