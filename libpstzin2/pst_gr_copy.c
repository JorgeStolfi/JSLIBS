/* See {pst_gr_copy.h} */
/* Last edited on 2025-03-13 04:48:42 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>

#include <pst_gr.h>

#define DELETED pst_gr_mark_DELETED
#define NONE pst_gr_mark_NONE

pst_gr_t *pst_gr_copy(pst_gr_t *gr, haf_arc_vec_t *C)
  {
    /* Collect the valid arcs: */
    

    auto pst_gr_arc_t linha( pst_gr_arc_t ed); 

    int32_t *eq_vector_vt = talloc(gr->NV, int32_t);
    int32_t valid_NV = 0;
    for (int32_t i = 0; i < gr->NV; i++)
      { pst_gr_vertex_data_t *vi = &(gr->vdata[i]);
        assert((vi->x == -1) == (vi->y == -1));
        if (vi->x != -1)
          { eq_vector_vt[i] = valid_NV;
            valid_NV++;
          }
        else
          { eq_vector_vt[i] = -1; }
      }
    int32_t *eq_vector_ed = talloc(gr->NE, int32_t);
    pst_gr_arc_t *ref_tab = talloc(2*gr->NE, pst_gr_arc_t);
    int32_t valid_NE  = 0;
    for (int32_t i = 0; i < gr->NE;i++)
      { if (gr->hedge[i] != NULL)
          { eq_vector_ed[i] = valid_NE;
            valid_NE++;
          }
        else
          { eq_vector_ed[i] = -1; }
      }
    pst_gr_t  *ng = pst_gr_new(valid_NV, valid_NE);
    /*Easy part, copy the raw data*/
    for (int32_t i = 0; i < gr->NV;i++)
      { pst_gr_vertex_data_t *vd = &(gr->vdata[i]);
        if (eq_vector_vt[i] != -1)
          {  pst_gr_add_vertex(ng, vd->x, vd->y, NULL, vd->coords ); }
      }
    for (int32_t i = 0; i < gr->NE;i++)
      { pst_gr_arc_t e = gr->hedge[i];
        if (eq_vector_ed[i] != -1)
          { int32_t ee = haf_edge_id(haf_edge(e));
            int32_t org =  pst_gr_arc_org(gr,e);
            int32_t dst = pst_gr_arc_org(gr,pst_gr_arc_sym(e));
            pst_gr_path_t P = pst_gr_arc_path(gr,e);
            double delta = pst_gr_arc_delta(gr,e);
            double weight = pst_gr_arc_weight(gr,e);

            char *label = gr->edata[ee]->label;
            char *nl = (label == NULL ? NULL : talloc(strlen(label)+1, char));
            if (label != NULL ) strcpy(nl,label);

            pst_gr_path_t np = P;
            if (P.v != NULL )
              { np.v = talloc(P.n, r2_t);
                memcpy(np.v,P.v,sizeof(r2_t)*(P.n));
              }
            /*We dont do much with it now...*/
            int32_t new_org = eq_vector_vt[org];
            int32_t new_dst = eq_vector_vt[dst];
            pst_gr_arc_t new_e =  pst_gr_add_edge(ng,new_org,new_dst,delta,weight,nl,np);
            int32_t ee_dir = haf_arc_id(e);
            int32_t ee_sym = haf_arc_id(pst_gr_arc_sym(e));
            ref_tab[ee_dir] = new_e;
            ref_tab[ee_sym] = pst_gr_arc_sym(new_e);
          }
      }
      
    pst_gr_arc_t linha( pst_gr_arc_t ed)
      { int32_t eid =  haf_arc_id(ed);
        return ref_tab[eid];
      }

    /*Hard part, assemble the graph data*/
    for (int32_t i = 0; i < gr->NE;i++)
      { pst_gr_arc_t e = gr->hedge[i];
        int32_t eel = eq_vector_ed[i];
        if (eel != -1)
          { haf_splice(haf_oprev(linha(e)),linha(haf_oprev(e)));
            haf_splice(haf_oprev(pst_gr_arc_sym(linha(e))),linha(haf_oprev(pst_gr_arc_sym(e))));
          }
      }

    free(eq_vector_vt);
    free(eq_vector_ed);
    free(ref_tab);
    return ng;
  }
