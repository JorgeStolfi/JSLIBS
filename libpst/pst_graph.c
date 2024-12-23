/* See {pst_graph.h}. */
/* Last edited on 2024-12-22 21:37:59 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <assert.h>
#include <affirm.h>
#include <jsfile.h>
#include <rn.h>

#include <pst_slope_map.h>
#include <pst_interpolate.h>

#include <pst_graph.h>

void pst_graph_update_neighbours(pst_graph_t* g);
  /* Updates the neighbours {vtxm,vtxp,vtym,vtyp} of all vertices,
    according to existing edges.  Only works for the 
    graph is derived from an image, so that each edge connects
    pixels that are adjacent horizontally or vertically. 
    Will fail if this is not the case. */

bool_t pst_graph_check_consistency(FILE* arq,pst_graph_t* g);
  /* Checks the graph consistency using the edge's list. If {arq} is not null,
    a humam-readable report is written to it. */

void pst_graph_write_vertex(FILE* arq, pst_vertex_t* v);
  /* Writes to {arq} the vertex {v} in human-readable format. */

void pst_graph_write_edge(FILE* arq, pst_edge_t* e);
  /* Writes to {arq} the edge {e} in human-readable format. */

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

/* IMPLEMENTATIONS */

int32_t pst_graph_compute_vertex_index(int32_t ix, int32_t iy, int32_t NX, int32_t NY)
  {
    if ((iy < 0 ) || (iy > NY)) { return -1; }
    if ((ix < 0 ) || (ix > NX)) { return -1; }
    return iy*NX + ix;
  }

void pst_graph_restore_vertex_index(int32_t id, uint32_t NX, uint32_t NY, int32_t *ix, int32_t *iy)
  {
    demand((NX >= 1) && (NY >= 1), "invalid image size");
    if ((id < 0 )  || (id >= NX*NY)) 
      { *ix = *iy = -1; }
    else
      { *ix = id%(int32_t)NX; *iy = id/(int32_t)NX; }
  }

void pst_graph_get_edge_difference_and_weight(pst_graph_t *g, uint32_t v, int32_t e, double *d, double *w)
  {
    demand(v < g->num_vertex, "invalid vertex index");
    if (e == -1) 
      { *w = 0; *d = 0; }
    else
      { demand(e < g->num_edge, "invaid edge index");
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

void pst_graph_interpolate_two_samples
  (  float_image_t *I, float_image_t *W,
     uint32_t c,
     int32_t x0, int32_t y0,
     int32_t x1, int32_t y1,
     double *v, double *w
   )
   {
     int32_t NX = (int32_t)I->sz[1]; 
     if (W != NULL) { assert(W->sz[1] == NX); }
     int32_t NY = (int32_t)I->sz[2]; 
     if (W != NULL) { assert(W->sz[2] == NY); }
     
     double v0 = ( (x0 < 0) || (y0 < 0) || (x0 >= NX) || (y0 >= NY) ? 0 : float_image_get_sample(I, (int32_t)c, x0, y0));
     double v1 = ( (x1 < 0) || (y1 < 0) || (x1 >= NX) || (y1 >= NY) ? 0 : float_image_get_sample(I, (int32_t)c, x1, y1));
     double w0 = ( (x0 < 0) || (y0 < 0) || (x0 >= NX) || (y0 >= NY) ? 0 : (W == NULL ? 1 : float_image_get_sample(W, 0, x0, y0)));
     double w1 = ( (x1 < 0) || (y1 < 0) || (x1 >= NX) || (y1 >= NY) ? 0 : (W == NULL ? 1 : float_image_get_sample(W, 0, x1, y1)));
   
     /* First we get the interpolation of one diagonal*/
     *v = (w0*v0+w1*v1)/2;
     *w = (w0+w1)/2;
   }
   
pst_graph_t *pst_graph_new(uint32_t max_vertex, uint32_t max_edge)
  {
    pst_graph_t *g = (pst_graph_t*)malloc(sizeof(pst_graph_t));
    g->num_vertex = 0;
    g->num_edge = 0;
    g->max_vertex = max_vertex;
    g->max_edge = max_edge;
    g->vertices = (pst_vertex_t*)malloc(sizeof(pst_vertex_t)*max_vertex);
    g->edges = (pst_edge_t*)malloc(sizeof(pst_edge_t)*max_edge);
    return g;
  }


void pst_graph_add_vertex(pst_graph_t *g, uint32_t id, int32_t x, int32_t y)
  {
    assert(g->num_vertex < g->max_vertex);
    if ((x < 0) || (y < 0)) { demand((x == -1) && (y == -1), "invalid pixel indices"); }
    pst_vertex_t v;
    v.id = id;
    v.x = x;
    v.y = y;
    v.vtxp = v.vtyp = v.vtxm = v.vtym = -1;
    g->vertices[g->num_vertex] = v;
    g->num_vertex++;
  }

void pst_graph_add_edge(pst_graph_t *g, uint32_t u, uint32_t v, double d, double w, int32_t axis)
  {
    assert( g->num_edge < g->max_edge);
    assert((u >= 0) && (u < g->max_vertex));
    assert((v >= 0) && (v < g->max_vertex));
    pst_edge_t e;
    e.u = u;
    e.v = v;
    e.d = d;
    e.w = w;
    e.axis = axis; 
    g->edges[g->num_edge]  = e;
    g->num_edge++;
  }

void pst_graph_update_neighbours(pst_graph_t *g)
  {
    for(int32_t i = 0; i < g->num_edge; i++)
      { uint32_t u = g->edges[i].u;
        uint32_t v = g->edges[i].v;
        assert((u >= 0) && (u < g->num_vertex));
        assert((v >= 0) && (v < g->num_vertex));
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

pst_graph_t *pst_graph_create_from_gradient_and_weight_maps(float_image_t *IG, float_image_t *IW)
  {
    int32_t NX = (int32_t)IG->sz[1];
    int32_t NY = (int32_t)IG->sz[2];

    int32_t NX_Z = NX+1;
    int32_t NY_Z = NY+1;

    assert(IG->sz[0] == 2);
    if (IW != NULL)
      { assert(IW->sz[0] == 1);
        assert((NX ==  IW->sz[1]) && (NY == IW->sz[2]));
      }

    uint32_t num_vertex_max = (uint32_t)(NX_Z*NY_Z);
    uint32_t num_edge_max = (uint32_t)(2*(NX_Z*NY_Z) - NX_Z - NY_Z);

    pst_graph_t *g = pst_graph_new(num_vertex_max, num_edge_max);

    for (int32_t y = 0; y < NY_Z; y++)
      { for (int32_t x = 0; x < NX_Z; x++)
          { int32_t ioo = pst_graph_compute_vertex_index(x, y, NX_Z, NY_Z);
            assert((ioo >= 0) && (ioo < num_vertex_max));
            pst_graph_add_vertex(g, (uint32_t)ioo, x, y);

            /* Edge from {(x,y)} to {(x-1,y)}: */
            double dxm, wxm;
            pst_interpolate_four_samples(IG, IW, 0, x-1, y-1, x-1, y+0, &dxm, &wxm);
            if (wxm > 0) 
              { int32_t imo = pst_graph_compute_vertex_index(x-1, y, NX_Z, NY_Z);
                assert((imo >= 0) && (imo < num_vertex_max));
                pst_graph_add_edge(g, (uint32_t)ioo, (uint32_t)imo, dxm, wxm, 0);
              }

            /* Edge from {(x,y)} to {(x,y-1)}: */
            double dym, wym;
            pst_interpolate_four_samples(IG, IW, 1, x-1, y-1, x+0, y-1, &dym, &wym);
            if (wym > 0) 
              { int32_t iom = pst_graph_compute_vertex_index(x, y-1, NX_Z, NY_Z);
                assert((iom >= 0) && (iom < num_vertex_max));
                pst_graph_add_edge(g, (uint32_t)ioo, (uint32_t)iom, dym, wym, 1);
              }
          }
      }

    pst_graph_update_neighbours(g);
    return g;
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
    fprintf(wr, "%d %d\n", g->num_vertex, g->num_edge);
    for (uint32_t i = 0; i < g->num_vertex; i++)
      { fprintf(wr, "[%d] ", i);
        pst_graph_write_vertex(wr, &(g->vertices[i]));
      }
    for (uint32_t i = 0; i < g->num_edge; i++)
      { fprintf(wr, "[%d] ", i);
        pst_graph_write_edge(wr, &(g->edges[i]));
      }
  }

bool_t pst_graph_check_consistency(FILE *wr, pst_graph_t *g)
  {
    if ( wr != NULL) { fprintf(wr, "checking edge consistency...\n"); }
    bool_t teste = TRUE;
    for (int32_t i = 0; i < g->num_edge; i++)
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

pst_graph_t *pst_graph_shrink(pst_graph_t *g)
  {
    /* We first count how much nodes we will take. This array 
      allows know the correspondence G->JG, so given a 
      vertex i in G the equivalent in JG is correspondence_vector[i]. */

    int32_t *correspondence_vector = talloc(g->num_vertex, int32_t);
    int32_t total_vertices = 0;
     
    for(uint32_t i = 0; i < g->num_vertex; i++)
      { int32_t xi = g->vertices[i].x;
        int32_t yi = g->vertices[i].y;
        if((xi + yi)%2 == 0)
          { correspondence_vector[i] = total_vertices; total_vertices++; }
        else
          { correspondence_vector[i] = -1; }
      }
   
    pst_graph_t *jg = pst_graph_new((uint32_t)total_vertices, (uint32_t)(2*total_vertices));
    
    for (uint32_t i = 0; i < g->num_vertex; i++)
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

uint32_t *pst_graph_spanning_forest(pst_graph_t *g)
  { bool_t *marked_vertices = talloc(g->num_vertex, bool_t);
    uint32_t *marked_edges = talloc(g->num_edge, uint32_t);
    for (uint32_t i = 0; i < g->num_vertex; i++) { marked_vertices[i] = FALSE; }
    for (uint32_t i = 0; i < g->num_edge; i++) { marked_edges[i] = 0; }
    uint32_t num_tree = 1;
    for (uint32_t i = 0; i < g->num_vertex; i++)
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
                assert((ed >= 0) && (ed < g->num_edge));
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
