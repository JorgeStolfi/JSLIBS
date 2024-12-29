/* See {pst_graph.h}. */
/* Last edited on 2024-12-23 09:47:59 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <assert.h>
#include <affirm.h>
#include <jsfile.h>
#include <rn.h>

#include <pst_graph.h>

#include <pst_graph_shrink.h>

pst_graph_t *pst_graph_shrink(pst_graph_t *g)
  {
    /* We first count how much nodes we will take. This array 
      allows know the correspondence G->JG, so given a 
      vertex i in G the equivalent in JG is correspondence_vector[i]. */

    int32_t *correspondence_vector = talloc(g->NV, int32_t);
    int32_t total_vertices = 0;
     
    for(uint32_t i = 0; i < g->NV; i++)
      { int32_t xi = g->vertices[i].x;
        int32_t yi = g->vertices[i].y;
        if((xi + yi)%2 == 0)
          { correspondence_vector[i] = total_vertices; total_vertices++; }
        else
          { correspondence_vector[i] = -1; }
      }
   
    pst_graph_t *jg = pst_graph_new((uint32_t)total_vertices, (uint32_t)(2*total_vertices));
    
    for (uint32_t i = 0; i < g->NV; i++)
      {
        int32_t xi = g->vertices[i].x;
        int32_t yi = g->vertices[i].y;

        if ((xi + yi)%2 == 0)
          { /*Create vertex copy*/
            int32_t ioo = correspondence_vector[i];
            assert((ioo > 0) && (ioo < total_vertices));
            pst_graph_add_vertex(jg, g->vertices[i].id, (xi + yi)/2, (xi - yi)/2);

            double dxm, wxm;
            int32_t vtxm;
            pst_graph_vertex_get_neighbour(g, i, -1, -1, &vtxm, &dxm, &wxm);
            double dym, wym;
            int32_t vtym;
            pst_graph_vertex_get_neighbour(g, i, +1, -1, &vtym, &dym, &wym); 
            /* The last ones are computed just for referencing sake, since we dont have x and y anymore*/

            if (wxm > 0)
              { assert((vtxm >= 0) && (vtxm < total_vertices));
                int32_t imo = correspondence_vector[vtxm];
                assert((imo >= 0) && (imo < total_vertices));
                pst_graph_add_edge(jg, (uint32_t)ioo, (uint32_t)imo, dxm, wxm, 0);
              }

            if (wym > 0)
              { assert((vtym >= 0) && (vtym < total_vertices));
                int32_t iom = correspondence_vector[vtym];
                assert((iom >= 0) && (iom < total_vertices));
                pst_graph_add_edge(jg, (uint32_t)ioo, (uint32_t)iom, dym, wym, 1);
              }
          }
      }

    free(correspondence_vector);
    pst_graph_update_neighbours(jg);
    
    return jg;
  }
