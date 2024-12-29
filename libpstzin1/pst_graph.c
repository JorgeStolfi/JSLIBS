/* See {pst_graph.h}. */
/* Last edited on 2024-12-23 14:53:50 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <assert.h>
#include <affirm.h>
#include <jsfile.h>
#include <rn.h>
#include <bool.h>

#include <pst_graph.h>

void pst_graph_write_vertex(FILE* wr, pst_vertex_t* v);
  /* Writes to {wr} the vertex {v} in human-readable format. */

void pst_graph_write_edge(FILE* wr, pst_edge_t* e);
  /* Writes to {wr} the edge {e} in human-readable format. */

void pst_graph_mst_walk_recursive
  ( pst_graph_t* g,
    uint32_t vertex,
    bool_t* marked_vertices,
    uint32_t* marked_edges,
    uint32_t num_tree
  );
  /* Enumerates all edges and vertices of the graph {g} in the connected component of {g} that is reachable from 
    the vertex {vertex}.  Expects that, on entry. {marked_vertices[u]=FALSE} and {marked_edges[e]=0}
    for every vertex {u} and every edge {e} in that component.  On exit, {marked_vertices[u] = TRUE}.
    and {marked_edges[e] = num_tree} for every {u} and {e} in that component. */

void pst_graph_get_edge_difference_and_weight(pst_graph_t *g, uint32_t v, int32_t e, double *d, double *w)
  {
    demand(v < g->NV, "invalid vertex index");
    if (e == -1) 
      { *w = 0; *d = 0; }
    else
      { demand(e < g->NE, "invaid edge index");
        if (g->edges[e].u == v)
          { *d = g->edges[e].d; }
        else if (g->edges[e].v == v) 
          { *d = -g->edges[e].d; }
        else
          { demand(FALSE, "vertex does not belong to the edge"); }
        *w = g->edges[e].w;
      }
  }

void pst_graph_vertex_get_neighbour
  ( pst_graph_t *g,
    uint32_t vt0,
    int32_t dx, int32_t dy,
    int32_t *vt1,
    double *d,
    double *w
  )
  {
    int32_t vt = -1;
    int32_t ed = -1;
    double we = 0;
    double de = 0;

    if((dx*dx + dy*dy) == 1)
      { /* Single step: */
        if( (dx == 1) && (dy == 0) )
          { ed = g->vertices[vt0].vtxp;
            if (ed != -1) { vt = (int32_t)g->edges[ed].u; }
          }
        if( (dx == -1) && (dy == 0) )
          { ed = g->vertices[vt0].vtxm;
            if (ed != -1) { vt = (int32_t)g->edges[ed].v; }
          }
        if( (dx == 0) && (dy == 1) )
          { ed = g->vertices[vt0].vtyp;
            if (ed != -1) { vt = (int32_t)g->edges[ed].u; }
          }
        if( (dx == 0) && (dy == -1) )
          { ed =  g->vertices[vt0].vtym;
            if (ed != -1) { vt = (int32_t)g->edges[ed].v; }
          }
        pst_graph_get_edge_difference_and_weight(g, vt0, ed, &de, &we);
      }
    else if ((dx*dx + dy*dy) == 2)
      { /* First, takes the horizontal path */
        int32_t h0, h1;
        h0 = h1 = -1;
        double hd0, hd1, hw0, hw1;
        hd0 = hd1 = hw0 = hw1 = 0;
        /* This will never be recursive because one of the axis is always equal to 0 */
        pst_graph_vertex_get_neighbour(g, vt0, dx, 0, &h0, &hd0, &hw0); 
        if (h0 != -1) { pst_graph_vertex_get_neighbour(g, (uint32_t)h0, 0, dy, &h1, &hd1, &hw1); }

        /* Now the vertical path*/
        int32_t v0, v1;
        v0 = v1 = -1;
        double vd0, vd1, vw0, vw1;
        vd0 = vd1 = vw0 = vw1 = 0;
        /* This will never be recursive because one of the axis is always equal to 0 */
        pst_graph_vertex_get_neighbour(g, vt0, 0, dy, &v0, &vd0, &vw0); 
        if (v0 != -1) { pst_graph_vertex_get_neighbour(g, (uint32_t)v0, dx, 0, &v1, &vd1, &vw1); }

        if ((v1 != -1) && (h1 != -1)) { demand(v1 == h1, "diagonal vertices are inconsistent!"); }

        double hd = hd0+hd1;
        double hw = 0.5/((1/hw0) + (1/hw1));
        double vd = vd0+vd1;
        double vw = 0.5/((1/vw0) + (1/vw1));

        if (h1 != -1) { vt = h1; }
        if (v1 != -1) { vt = v1; }

        we = hw+vw;
        de = (hw*hd + vw*vd)/we;
        if (we == 0) { de = (hd+vd)/2.0; }
      }
    else
      { demand(FALSE, "cannot search neighbour that is not close by"); }

    *vt1 = vt;
    *d = de;
    *w = we;
  }

pst_graph_t *pst_graph_new(uint32_t NV_max, uint32_t NE_max)
  {
    pst_graph_t *g = (pst_graph_t*)malloc(sizeof(pst_graph_t));
    g->NV = 0;
    g->NE = 0;
    g->NV_max = NV_max;
    g->NE_max = NE_max;
    g->vertices = (pst_vertex_t*)malloc(sizeof(pst_vertex_t)*NV_max);
    g->edges = (pst_edge_t*)malloc(sizeof(pst_edge_t)*NE_max);
    return g;
  }


void pst_graph_add_vertex(pst_graph_t *g, uint32_t id, int32_t x, int32_t y)
  {
    assert(g->NV < g->NV_max);
    if ((x < 0) || (y < 0)) { demand((x == -1) && (y == -1), "invalid pixel indices"); }
    pst_vertex_t v;
    v.id = id;
    v.x = x;
    v.y = y;
    v.vtxp = v.vtyp = v.vtxm = v.vtym = -1;
    g->vertices[g->NV] = v;
    g->NV++;
  }

void pst_graph_add_edge(pst_graph_t *g, uint32_t u, uint32_t v, double d, double w, int32_t axis)
  {
    assert( g->NE < g->NE_max);
    assert((u >= 0) && (u < g->NV_max));
    assert((v >= 0) && (v < g->NV_max));
    pst_edge_t e;
    e.u = u;
    e.v = v;
    e.d = d;
    e.w = w;
    e.axis = axis; 
    g->edges[g->NE]  = e;
    g->NE++;
  }

void pst_graph_update_neighbours(pst_graph_t *g)
  {
    for(int32_t i = 0; i < g->NE; i++)
      { uint32_t u = g->edges[i].u;
        uint32_t v = g->edges[i].v;
        assert((u >= 0) && (u < g->NV));
        assert((v >= 0) && (v < g->NV));
        if (g->edges[i].axis == 0)
          { g->vertices[u].vtxm = i;
            g->vertices[v].vtxp = i;
          }
        else if (g->edges[i].axis == 1)
          { g->vertices[u].vtym = i;
            g->vertices[v].vtyp = i;
          }
        else
          { demand(FALSE, "invalid edge axis"); }
      }
  }

void pst_graph_write_vertex(FILE *wr, pst_vertex_t *v)
  {
    fprintf
      ( wr, "ID %d I %d J %d XM %d XP %d YM %d YP %d\n",
        v->id, v->x, v->y, v->vtxm, v->vtxp, v->vtym, v->vtyp 
      );
  }
  
void pst_graph_write_edge(FILE *wr, pst_edge_t *e)
  {
    fprintf(wr, "U %d V %d D %9.6f W %9.6f\n", e->u, e->v, e->d, e->w);
  }

void pst_graph_write(FILE *wr, pst_graph_t *g)
  {
    fprintf(wr, "%d %d\n", g->NV, g->NE);
    for (uint32_t i = 0; i < g->NV; i++)
      { fprintf(wr, "[%d] ", i);
        pst_graph_write_vertex(wr, &(g->vertices[i]));
      }
    for (uint32_t i = 0; i < g->NE; i++)
      { fprintf(wr, "[%d] ", i);
        pst_graph_write_edge(wr, &(g->edges[i]));
      }
  }

bool_t pst_graph_check_consistency(FILE *wr, pst_graph_t *g)
  {
    if ( wr != NULL) { fprintf(wr, "checking edge consistency...\n"); }
    bool_t teste = TRUE;
    for (int32_t i = 0; i < g->NE; i++)
    {
      uint32_t u = g->edges[i].u;
      uint32_t v = g->edges[i].v;
      
      /* Check u: */
      if (g->vertices[u].vtxm == i)
        { /* Must be the same in v: */
          if (g->vertices[v].vtxp != i)
            { teste = FALSE;
              if ( wr != NULL) 
                { fprintf(wr, "edge %d is VTXM at U %d but not registered properly at V %d\n", i, u, v); }
              /* Fix:  */
              g->vertices[v].vtxp = i;
            }
        }
      else if (g->vertices[u].vtym == i)
        { if (g->vertices[v].vtyp != i)
            { teste = FALSE;
              if (wr != NULL)
                { fprintf(wr, "edge %d is VTYM at U %d but not registered properly at V %d\n", i, u, v); }
              g->vertices[v].vtyp = i;
            }
        }
      else
        { teste = FALSE;
          if (wr != NULL) { fprintf(wr, "edge %d :is not registered at U =  %d\n", i, u); }
        }

      if (g->vertices[v].vtxp == i)
        { /* Must be the same in v */
          if (g->vertices[u].vtxm != i)
            { teste = FALSE;
              if (wr != NULL) 
                { fprintf(wr, "edge %d :is VTXP at V %d but not registered properly at U %d\n", i, v, u); }
              g->vertices[u].vtxm = i;
            }
        }
      else if (g->vertices[v].vtyp == i)
        { if(g->vertices[u].vtym != i)
            { teste = FALSE;
              if (wr != NULL) { fprintf(wr, "edge %d :is VTYP V %d but not registered properly at U %d\n", i, v, u); }
              g->vertices[v].vtym = i;
            }
        }
      else
        { teste = FALSE;
          if (wr != NULL) { fprintf(wr, "edge %d :is not registered at V =  %d\n", i, u); }
        }
      if( ! teste){
        if (wr != NULL) { fprintf(wr, "failed Edge %d ( %d -> %d)  \n", i, u, v);}
      }
    }

    if (teste)
      { if (wr != NULL) { fprintf(wr, "OK\n"); } }
    else
      { if (wr != NULL) { fprintf(wr, "FAILED !\n"); } }

    return teste;
  }

uint32_t *pst_graph_spanning_forest(pst_graph_t *g)
  { bool_t *marked_vertices = talloc(g->NV, bool_t);
    uint32_t *marked_edges = talloc(g->NE, uint32_t);
    for (uint32_t i = 0; i < g->NV; i++) { marked_vertices[i] = FALSE; }
    for (uint32_t i = 0; i < g->NE; i++) { marked_edges[i] = 0; }
    uint32_t num_tree = 1;
    for (uint32_t i = 0; i < g->NV; i++)
      { if (! marked_vertices[i])
          { pst_graph_mst_walk_recursive(g, i, marked_vertices, marked_edges, num_tree);
            num_tree++;
          }
      }
    free(marked_vertices);
    return marked_edges;
  }

void pst_graph_mst_walk_recursive
  ( pst_graph_t *g,
    uint32_t vertex,
    bool_t *marked_vertices,
    uint32_t *marked_edges,
    uint32_t num_tree
  )
  {
    assert(! marked_vertices[vertex]);
    marked_vertices[vertex] = TRUE; 

    double d, w;
    for (uint32_t k = 0; k < 4; k++)
      { int32_t dx, dy;
        switch (k)
          { case 0: dx = 00; dy = +1; break;
            case 1: dx = 00; dy = -1; break;
            case 2: dx = +1; dy = 00; break;
            case 3: dx = -1; dy = 00; break;
            default: assert(FALSE);
          }
        int32_t neighbour;
        pst_graph_vertex_get_neighbour(g, vertex, dx, dy, &neighbour, &d, &w);
        if (neighbour != -1)
          { if (! marked_vertices[neighbour])
              { int32_t ed;
                switch (k)
                  { case 0: ed = g->vertices[vertex].vtyp; break;
                    case 1: ed = g->vertices[vertex].vtym; break;
                    case 2: ed = g->vertices[vertex].vtxp; break;
                    case 3: ed = g->vertices[vertex].vtxm; break;
                    default: assert(FALSE);
                  }
                assert((ed >= 0) && (ed < g->NE));
                marked_edges[ed] = num_tree;
                pst_graph_mst_walk_recursive(g, (uint32_t)neighbour, marked_vertices, marked_edges, num_tree);
              }
          }
     }
  }

void pst_graph_free(pst_graph_t *g)
  { free(g->vertices);
    free(g->edges);
    free(g);
  }
