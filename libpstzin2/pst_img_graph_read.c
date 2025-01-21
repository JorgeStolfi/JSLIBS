/* See {pst_img_graph_read.h} */
/* Last edited on 2025-01-13 07:51:07 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <haf.h>
#include <haf_read.h>
#include <filefmt.h>
#include <affirm.h>

#include <pst_img_graph.h>

#include <pst_img_graph_read.h>

#define DELETED pst_img_graph_mark_DELETED

void pst_img_graph_read_vertex(FILE* rd, pst_img_graph_t * g);

void pst_img_graph_read_edge(FILE* rd, pst_img_graph_t *g, int32_t list_onext[], int32_t list_long_bit[]);

/* IMPLEMENTATIONS */

void pst_img_graph_read_vertex(FILE *wr, pst_img_graph_t  *g)
  {
    int32_t id;
    int32_t vmark;
    r2_t coords;
    int32_t eedge;
    int32_t nsc = fscanf(wr,"%d %d %lf %lf %d",&id,&vmark,&(coords.c[0]),&(coords.c[1]),&eedge);
    demand(nsc == 5, "Cannot read vertex");
    /*Who is responsible to fill the  real edge's value is the pst_img_read_edges*/
    int32_t ix = pst_img_graph_add_vertex(g,id,NULL,coords);
    g->vertex[ix].vmark = vmark;
  }

void pst_img_graph_read_edge
  ( FILE *wr,
    pst_img_graph_t *g,
    int32_t list_onext[],
    int32_t list_long_bit[]
  )
  { int32_t eid;
    int32_t aid;
    int32_t org ;
    int32_t dst;
    double delta;
    double weight;

    demand(fscanf(wr,"%d",&eid) == 1, "Cannot read edge index");
    if (eid != -1)
      {
        int32_t ns1 = fscanf(wr,"%d %d %d %lf %lf",&aid,&org,&dst,&delta,&weight);
        demand(ns1 == 5, "Cannot read edge data");
        int32_t ns2 = fscanf(wr,"%d %d %d %d",&(list_onext[0]),&(list_long_bit[0]),&(list_onext[1]),&(list_long_bit[1]));
        demand(ns2 == 4, "Cannot read edge connectivity");

        int32_t label_len;
        int32_t ns3 = fscanf(wr,"%d",&(label_len));
        demand(ns3 == 1, "Cannot read edge label lenght");
        demand(label_len >= 0, "invalid label length");
        char *label = (label_len == 0 ? NULL : talloc(label_len +1, char));
        if (label_len > 0)
          { fgetc(wr); /*skip the first blank space*/
            for (int32_t i = 0; i < label_len; i++)
              { fscanf(wr, "%c",&(label[i])); }
            label[label_len] = '\0';
          }
        pst_path_t p = pst_img_graph_read_path(wr);
        pst_img_graph_edge_add(g,org,dst,delta,weight,label,p);
      }
    else
      { g->{hedge|dedge}@@[g->NE].arc = NULL;
        g->NE++;
      }
  }

pst_img_graph_t* pst_img_graph_read(FILE *wr)
  {
    int32_t NV, NE;
    demand( fscanf(wr,"%d %d", &NV, &NE) == 4, "parse of vertex and edge counts failed");
    demand((NV >= 0) && (NE >=0), "invalid vertex or edge counts");

    pst_img_graph_t *g = pst_img_graph_new(NV, NE);

    for (int32_t i = 0; i < NV ; i++) { pst_img_graph_read_vertex(wr,g); }
    int32_t *list_next = talloc(2*NE, int32_t);
    int32_t *list_long_bit = talloc(2*NE, int32_t);
    for (int32_t i = 0; i < NE ; i++)
      { pst_img_graph_read_edge(wr,g,&(list_next[2*i]),&(list_long_bit[2*i])); }

    auto haf_arc_t recover_edge(int32_t index, int32_t long_bit);

    haf_arc_t  recover_edge(int32_t index, int32_t long_bit)
      { assert((index >= 0 ) && (index < g->NE));
        haf_arc_t e = g->{hedge|dedge}@@[index].arc;
        assert(e != NULL);
        if (long_bit == 1) { e = haf_sym(e); }
        return e;
      }

    /*Hard part*/
    for (int32_t i = 0; i < g->NE;i++)
      { haf_arc_t e = g->{hedge|dedge}@@[i].arc;
        if (e != NULL)
          { haf_arc_t e_next = recover_edge(list_next[2*i],list_long_bit[2*i]);
            haf_arc_t e_next_sym = recover_edge(list_next[(2*i)+1],list_long_bit[(2*i)+1]);

            haf_splice( haf_oprev(e_next),e);
            haf_splice( haf_oprev(e_next_sym),haf_sym(e));
          }
      }
    free(list_next);
    free(list_long_bit);
    return g;
  }
