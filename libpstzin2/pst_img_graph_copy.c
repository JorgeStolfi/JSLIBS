/* See {pst_img_graph_copy.h} */
/* Last edited on 2025-01-10 07:50:09 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <haf.h>
#include <float_image.h>
#include <jsfile.h>
#include <filefmt.h>
#include <affirm.h>
#include <rn.h>

#include <pst_img_graph.h>

#define DELETED pst_img_graph_mark_DELETED
#define NONE pst_img_graph_mark_NONE

pst_img_graph_t *pst_img_graph_copy(pst_img_graph_t *g, haf_arc_vec_t *C)
  {
    /* Collect the valid arcs: */
    

    auto haf_arc_t linha( haf_arc_t ed); 

    int32_t *eq_vector_vt = talloc(g->NV, int32_t);
    int32_t valid_NV = 0;
    for (int32_t i = 0; i < g->NV; i++)
      { pst_img_graph_vertex_data_t *vi = &(g->vdata[i]);
        assert((vi->x == -1) == (vi->y == -1));
        if (vi->x != -1)
          { eq_vector_vt[i] = valid_NV;
            valid_NV++;
          }
        else
          { eq_vector_vt[i] = -1; }
      }
    int32_t *eq_vector_ed = talloc(g->NE, int32_t);
    haf_arc_t *ref_tab = talloc(2*g->NE, haf_arc_t);
    int32_t valid_NE  = 0;
    for (int32_t i = 0; i < g->NE;i++)
      { if (g->hedge[i] != NULL)
          { eq_vector_ed[i] = valid_NE;
            valid_NE++;
          }
        else
          { eq_vector_ed[i] = -1; }
      }
    pst_img_graph_t  *ng = pst_img_graph_new(valid_NV, valid_NE);
    /*Easy part, copy the raw data*/
    for (int32_t i = 0; i < g->NV;i++)
      { pst_img_graph_vertex_data_t *vd = &(g->vdata[i]);
        if (eq_vector_vt[i] != -1)
          {  pst_img_graph_add_vertex(ng, vd->x, vd->y, NULL, vd->coords ); }
      }
    for (int32_t i = 0; i < g->NE;i++)
      { haf_arc_t e = g->hedge[i];
        if (eq_vector_ed[i] != -1)
          { int32_t ee = haf_edge_id(haf_edge(e));
            int32_t org =  pst_img_graph_get_arc_origin(g,e);
            int32_t dst = pst_img_graph_get_arc_origin(g,haf_sym(e));
            pst_path_t p = pst_img_graph_get_edge_path(g,e);
            double delta = pst_img_graph_get_arc_delta(g,e);
            double weight = pst_img_graph_get_edge_weight(g,e);

            char *label = g->edata[ee]->label;
            char *nl = (label == NULL ? NULL : talloc(strlen(label)+1, char));
            if (label != NULL ) strcpy(nl,label);

            pst_path_t np = p;
            if (p.v != NULL )
              { np.v = talloc(p.n, r2_t);
                memcpy(np.v,p.v,sizeof(r2_t)*(p.n));
              }
            /*We dont do much with it now...*/
            int32_t new_org = eq_vector_vt[org];
            int32_t new_dst = eq_vector_vt[dst];
            haf_arc_t new_e =  pst_img_graph_edge_add(ng,new_org,new_dst,delta,weight,nl,np);
            int32_t ee_dir = haf_arc_id(e);
            int32_t ee_sym = haf_arc_id(haf_sym(e));
            ref_tab[ee_dir] = new_e;
            ref_tab[ee_sym] = haf_sym(new_e);
          }
      }
      
    haf_arc_t linha( haf_arc_t ed)
      { int32_t eid =  haf_arc_id(ed);
        return ref_tab[eid];
      }

    /*Hard part, assemble the graph data*/
    for (int32_t i = 0; i < g->NE;i++)
      { haf_arc_t e = g->hedge[i];
        int32_t eel = eq_vector_ed[i];
        if (eel != -1)
          { haf_splice(haf_oprev(linha(e)),linha(haf_oprev(e)));
            haf_splice(haf_oprev(haf_sym(linha(e))),linha(haf_oprev(haf_sym(e))));
          }
      }

    free(eq_vector_vt);
    free(eq_vector_ed);
    free(ref_tab);
    return ng;
  }
