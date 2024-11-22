/* See {pst_graph.h}. */
/* Last edited on 2023-02-25 16:09:14 by stolfi */
/* Created by Rafael F. V. Saracchini */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>
#include <affirm.h>
#include <jsfile.h>
#include <rn.h>

#include <pst_slope_map.h>
#include <pst_interpolate.h>

#include <pst_graph.h>

long int pst_graph_compute_vertex_index(long int ix, long int iy, long int NX, long int NY)
  {
    if ((iy < 0 ) || (iy > NY)) { return -1; }
    if ((ix < 0 ) || (ix > NX)) { return -1; }
    return iy*NX + ix;
  }

void pst_graph_restore_vertex_index(long int id, long int NX, long int NY, long int *ix, long int *iy)
  {
    if ((id < 0 )  || (id >= NX*NY)) { *ix = *iy = -1; }
    *ix = id%NX;
    *iy = id/NX;
  }

void pst_graph_get_edge_derivate_weight(pst_graph_t *g, long int v, long int e, double *d, double *w)
  {
    *w = 0;
    *d = 0;
    if (e == -1) return;

    if (g->edges[e].u == v) { *d = g->edges[e].d; }
    else if (g->edges[e].v == v) { *d = -g->edges[e].d; }
    else { demand(FALSE, "requested vertex does not belong to a edge"); }
    *w = g->edges[e].w;
  }

void pst_graph_vertex_get_neighbour
  ( pst_graph_t *g,
    long int vt0,
    int dx, int dy,
    long int *vt1,
    double *d,
    double *w
  )
  {
    long int vt = -1;
    long int ed = -1;
    double we = 0;
    double de = 0;

    if((dx*dx + dy*dy) == 1)
      { if( (dx == 1) && (dy == 0) )
          { ed = g->vertices[vt0].vtxp;
            if (ed != -1) { vt = g->edges[ed].u; }
          }
        if( (dx == -1) && (dy == 0) )
          { ed = g->vertices[vt0].vtxm;
            if (ed != -1) { vt = g->edges[ed].v; }
          }
        if( (dx == 0) && (dy == 1) )
          { ed = g->vertices[vt0].vtyp;
            if (ed != -1) { vt = g->edges[ed].u; }
          }
        if( (dx == 0) && (dy == -1) )
          { ed =  g->vertices[vt0].vtym;
            if (ed != -1) { vt = g->edges[ed].v; }
          }
        pst_graph_get_edge_derivate_weight(g, vt0, ed, &de, &we);
      }
    else if ((dx*dx + dy*dy) == 2)
      { /* First, takes the horizontal path */
        long int h0, h1;
        h0 = h1 = -1;
        double hd0, hd1, hw0, hw1;
        hd0 = hd1 = hw0 = hw1 = 0;
        /* This will never be recursive because one of the axis is always equal to 0 */
        pst_graph_vertex_get_neighbour(g, vt0, dx, 0, &h0, &hd0, &hw0); 
        if (h0 != -1)
          { pst_graph_vertex_get_neighbour(g, h0, 0, dy, &h1, &hd1, &hw1); }

        /* Now the vertical path*/
        long int v0, v1;
        v0 = v1 = -1;
        double vd0, vd1, vw0, vw1;
        vd0 = vd1 = vw0 = vw1 = 0;
        /* This will never be recursive because one of the axis is always equal to 0 */
        pst_graph_vertex_get_neighbour(g, vt0, 0, dy, &v0, &vd0, &vw0); 
        if (v0 != -1)
          { pst_graph_vertex_get_neighbour(g, v0, dx, 0, &v1, &vd1, &vw1); }

        if ((v1 != -1) && (h1 != -1)) 
          { demand(v1 == h1, "diagonal vertices are inconsistent!"); }

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
     int c,
     int x0, int y0,
     int x1, int y1,
     double *v, double *w
   )
   {
     int NX = (int)I->sz[1]; 
     if (W != NULL) { assert(W->sz[1] == NX); }
     int NY = (int)I->sz[2]; 
     if (W != NULL) { assert(W->sz[2] == NY); }
     
     double v0 = ( (x0 < 0) || (y0 < 0) || (x0 >= NX) || (y0 >= NY) ? 0 : float_image_get_sample(I, c, x0, y0));
     double v1 = ( (x1 < 0) || (y1 < 0) || (x1 >= NX) || (y1 >= NY) ? 0 : float_image_get_sample(I, c, x1, y1));
     double w0 = ( (x0 < 0) || (y0 < 0) || (x0 >= NX) || (y0 >= NY) ? 0 : (W == NULL ? 1 : float_image_get_sample(W, 0, x0, y0)));
     double w1 = ( (x1 < 0) || (y1 < 0) || (x1 >= NX) || (y1 >= NY) ? 0 : (W == NULL ? 1 : float_image_get_sample(W, 0, x1, y1)));
   
     /* First we get the interpolation of one diagonal*/
     *v = (w0*v0+w1*v1)/2;
     *w = (w0+w1)/2;
   }
   
pst_graph_t *pst_graph_new(long int max_vertex, long int max_edge)
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


void pst_graph_add_vertex(pst_graph_t *g, long int id, long int i, long int j)
  {
    assert(g->num_vertex < g->max_vertex);
    pst_vertex_t v;
    v.id = id;
    v.i = i;
    v.j = j;
    v.vtxp = v.vtyp = v.vtxm = v.vtym = -1;
    g->vertices[g->num_vertex] = v;
    g->num_vertex++;
  }

void pst_graph_add_edge(pst_graph_t *g, long int u, long int v, double d, double w, int axis)
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
    long int i;
    for(i = 0; i < g->num_edge; i++)
      { long int u = g->edges[i].u;
        long int v = g->edges[i].v;
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
          { demand(FALSE, "invalid axis for the graph"); }
      }
  }

pst_graph_t *pst_graph_create_from_gradient_weight_map(float_image_t *IG, float_image_t *IW)
  {
    long int NX = (int)IG->sz[1];
    long int NY = (int)IG->sz[2];

    long int NX_Z = NX+1;
    long int NY_Z = NY+1;

    assert(IG->sz[0] == 2);
    if (IW != NULL)
      { assert(IW->sz[0] == 1);
        assert((NX ==  IW->sz[1]) && (NY == IW->sz[2]));
      }

    long int num_vertex_max = NX_Z*NY_Z;
    long int num_edge_max = 2*(NX_Z*NY_Z) - NX_Z - NY_Z;

    pst_graph_t *g = pst_graph_new(num_vertex_max, num_edge_max);

    long int x, y;
    for (y = 0; y < NY_Z; y++)
      { for (x = 0; x < NX_Z; x++)
          { long int ioo = pst_graph_compute_vertex_index(x, y, NX_Z, NY_Z);
            pst_graph_add_vertex(g, ioo, x, y);

            /* South and West edges, just in case, this order ensures that the edges goes only to already created vertices*/
            double dxm, wxm;
            double dym, wym;

            pst_interpolate_four_samples(IG, IW, 1, (int)(x-1), (int)(y-1), (int)(x+0), (int)(y-1), &dym, &wym);
            pst_interpolate_four_samples(IG, IW, 0, (int)(x-1), (int)(y-1), (int)(x-1), (int)(y+0), &dxm, &wxm);

           if (wxm > 0) 
             { long int imo = pst_graph_compute_vertex_index(x-1, y, NX_Z, NY_Z);
               assert(imo != -1);
               pst_graph_add_edge(g, ioo, imo, dxm, wxm, 0);
             }

           if(wym > 0) 
             { long int iom = pst_graph_compute_vertex_index(x, y-1, NX_Z, NY_Z);
               assert(iom != -1);
               pst_graph_add_edge(g, ioo, iom, dym, wym, 1);
             }
          }
      }

    pst_graph_update_neighbours(g);
    return g;
  }

void pst_graph_write_vertex(FILE *wr, pst_vertex_t *v)
  {
    fprintf
      ( wr, "ID %ld I %ld J %ld XM %ld XP %ld YM %ld YP %ld\n",
        v->id, v->i, v->j, v->vtxm, v->vtxp, v->vtym, v->vtyp 
      );
  }
  
void pst_graph_write_edge(FILE *wr, pst_edge_t *e)
  {
    fprintf(wr, "U %ld V %ld D %9.6lf W %9.6lf\n", e->u, e->v, e->d, e->w);
  }

void pst_graph_write(FILE *wr, pst_graph_t *g)
  {
    fprintf(wr, "%ld %ld\n", g->num_vertex, g->num_edge);
    int i;
    for (i = 0 ; i < g->num_vertex; i++)
      { fprintf(wr, "[%d] ", i);
        pst_graph_write_vertex(wr, &(g->vertices[i]));
      }
    for (i = 0 ; i < g->num_edge; i++)
      { fprintf(wr, "[%d] ", i);
        pst_graph_write_edge(wr, &(g->edges[i]));
      }
  }

bool_t pst_graph_check_consistency(FILE *wr, pst_graph_t *g)
  {
    long int i;
    if ( wr != NULL) { fprintf(wr, "checking edge consistency...\n"); }
    int teste = TRUE;
    for (i = 0; i < g->num_edge; i++)
    {
      long int u = g->edges[i].u;
      long int v = g->edges[i].v;
      
      /* Check u: */
      if (g->vertices[u].vtxm == i)
        { /* Must be the same in v: */
          if (g->vertices[v].vtxp != i)
            { teste = FALSE;
              if ( wr != NULL) 
                { fprintf(wr, "edge %ld is VTXM at U %ld but not registered properly at V %ld\n", i, u, v); }
              /* Fix:  */
              g->vertices[v].vtxp = i;
            }
        }
      else if (g->vertices[u].vtym == i)
        { if (g->vertices[v].vtyp != i)
            { teste = FALSE;
              if (wr != NULL)
                { fprintf(wr, "edge %ld is VTYM at U %ld but not registered properly at V %ld\n", i, u, v); }
              g->vertices[v].vtyp = i;
            }
        }
      else
        { teste = FALSE;
          if (wr != NULL) { fprintf(wr, "edge %ld :is not registered at U =  %ld\n", i, u); }
        }

      if (g->vertices[v].vtxp == i)
        { /* Must be the same in v */
          if (g->vertices[u].vtxm != i)
            { teste = FALSE;
              if (wr != NULL) 
                { fprintf(wr, "edge %ld :is VTXP at V %ld but not registered properly at U %ld\n", i, v, u); }
              g->vertices[u].vtxm = i;
            }
        }
      else if (g->vertices[v].vtyp == i)
        { if(g->vertices[u].vtym != i)
            { teste = FALSE;
              if (wr != NULL) { fprintf(wr, "edge %ld :is VTYP V %ld but not registered properly at U %ld\n", i, v, u); }
              g->vertices[v].vtym = i;
            }
        }
      else
        { teste = FALSE;
          if (wr != NULL) { fprintf(wr, "edge %ld :is not registered at V =  %ld\n", i, u); }
        }
      if( ! teste){
        if (wr != NULL) { fprintf(wr, "failed Edge %ld ( %ld -> %ld)  \n", i, u, v);}
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
    long int *correspondence_vector = (long int*) malloc(sizeof(long int)*g->num_vertex);
    long int total_vertices = 0;
    long int i;
     
    for(i = 0; i < g->num_vertex; i++)
      { long int ii = g->vertices[i].i;
        long int ij = g->vertices[i].j;
        if((ii + ij)%2 == 0)
          { correspondence_vector[i] = total_vertices; total_vertices++; }
        else
          { correspondence_vector[i] = -1; }
      }
   
    pst_graph_t *jg = pst_graph_new(total_vertices, 2*total_vertices);
    
    for (i = 0; i < g->num_vertex; i++)
      {
        long int ii = g->vertices[i].i;
        long int ij = g->vertices[i].j;

        if ((ii + ij)%2 == 0)
          { /*Create vertex copy*/
            long int ioo = correspondence_vector[i];
            pst_graph_add_vertex(jg, g->vertices[i].id, (ii + ij)/2, (ii - ij)/2);

            double dxm, wxm;
            long int vtxm;
            pst_graph_vertex_get_neighbour(g, i, -1, -1, &vtxm, &dxm, &wxm);
            double dym, wym;
            long int vtym;
            pst_graph_vertex_get_neighbour(g, i, +1, -1, &vtym, &dym, &wym); 
            /* The last ones i compute just for referencing sake, since i dont have x and y anymore*/

            if (wxm > 0)
              { assert(vtxm != -1);
                long int imo = correspondence_vector[vtxm];
                pst_graph_add_edge(jg, ioo, imo, dxm, wxm, 0);
              }

            if (wym > 0)
              { assert(vtym != -1);
                long int iom = correspondence_vector[vtym];
                pst_graph_add_edge(jg, ioo, iom, dym, wym, 1);
              }
          }
      }

    free(correspondence_vector);
    pst_graph_update_neighbours(jg);
    
    return jg;
  }

void pst_graph_mst_walk_recursive(pst_graph_t *g, long int vertex, int *marked_vertices, int *marked_edges, int num_tree)
  {
    marked_vertices[vertex] = 1; 

    long int neighbour;
    double d, w;
    /* North: */
    pst_graph_vertex_get_neighbour(g, vertex, 0, 1, &neighbour, &d, &w);
    if (neighbour != -1)
      { if (marked_vertices[neighbour] == 0)
          { long int ed = g->vertices[vertex].vtyp;
            marked_edges[ed] = num_tree;
            pst_graph_mst_walk_recursive(g, neighbour, marked_vertices, marked_edges, num_tree);
          }
      }
    /* East; */
    pst_graph_vertex_get_neighbour(g, vertex, 1, 0, &neighbour, &d, &w);
    if (neighbour != -1)
      { if (marked_vertices[neighbour] == 0)
          { long int ed = g->vertices[vertex].vtxp;
            marked_edges[ed] = num_tree;
            pst_graph_mst_walk_recursive(g, neighbour, marked_vertices, marked_edges, num_tree);
          }
      }

    /* South: */
    pst_graph_vertex_get_neighbour(g, vertex, 0, -1, &neighbour, &d, &w);
    if (neighbour != -1)
      { if (marked_vertices[neighbour] == 0)
        { long int ed = g->vertices[vertex].vtym;
          marked_edges[ed] = num_tree;
          pst_graph_mst_walk_recursive(g, neighbour, marked_vertices, marked_edges, num_tree);
        }
      }
    /* West: */
    pst_graph_vertex_get_neighbour(g, vertex, -1, 0, &neighbour, &d, &w);
    if (neighbour != -1)
      { if (marked_vertices[neighbour] == 0)
          { long int ed = g->vertices[vertex].vtxm;
            marked_edges[ed] = num_tree;
            pst_graph_mst_walk_recursive(g, neighbour, marked_vertices, marked_edges, num_tree);
          }
      }
  }

int *pst_graph_compute_mst(pst_graph_t *g)
  { int *marked_vertices = (int*)malloc(sizeof(int)*g->num_vertex);
    int *marked_edges = (int*)malloc(sizeof(int)*g->num_edge);
    int i;
    for (i = 0; i < g->num_vertex; i++) { marked_vertices[i] = 0; }
    for (i = 0; i < g->num_edge; i++) { marked_edges[i] = 0; }
    int num_tree = 1;
    for (i = 0; i < g->num_vertex; i++)
      { if(marked_vertices[i] == 0)
          { pst_graph_mst_walk_recursive(g, i, marked_vertices, marked_edges, num_tree);
            num_tree++;
          }
      }
    free(marked_vertices);
    return marked_edges;
  }

pst_imgsys_t *pst_graph_build_integration_system(pst_graph_t *g, long int NX_Z, long int NY_Z)
  {
    long int NXY_Z = NX_Z*NY_Z;
    long int *ind_ix = (long int*)malloc(sizeof(long int)*g->num_vertex);
    pst_imgsys_equation_t *eq = (pst_imgsys_equation_t*)notnull(malloc((g->num_vertex)*sizeof(pst_imgsys_equation_t)), "no mem");
    int i;
    int N = 0;
    for (i = 0; i < g->num_vertex; i++)
      { pst_imgsys_equation_t *eqk = &(eq[N]);

        int nt = 0; /* Number of terms in equation. */
        eqk->rhs = 0.0;
        eqk->ix[nt] = i; eqk->cf[nt] = 1.00; nt++;
        ind_ix[i] = N;
        pst_vertex_t *v = &(g->vertices[i]);
        int has_xp = (v->vtxp != -1); 
        int has_xm = (v->vtxm != -1);
        int has_yp = (v->vtyp != -1);
        int has_ym = (v->vtym != -1);

        int num_neighbours = has_xp+ has_xm + has_yp+has_ym;
        if (num_neighbours > 0)
          { if (has_xm)
              { long int vxm;
                double dxm, wxm;
                pst_graph_vertex_get_neighbour(g, i, -1, 0, &vxm, &dxm, &wxm);
                assert(vxm != -1);
                eqk->ix[nt] = (int)(+vxm);
                eqk->cf[nt] = (int)(-wxm);
                eqk->rhs += +wxm*dxm;
                nt++;
              }
            if (has_xp)
              { long int vxp;
                double dxp, wxp;
                pst_graph_vertex_get_neighbour(g, i, 1, 0, &vxp, &dxp, &wxp);
                assert(vxp != -1);
                eqk->ix[nt] = (int)(+vxp);
                eqk->cf[nt] = (int)(-wxp);
                eqk->rhs += +wxp*dxp;
                nt++;
              }
            if (has_ym)
              { long int vym;
                double dym, wym;
                pst_graph_vertex_get_neighbour(g, i, 0, -1, &vym, &dym, &wym);
                assert(vym != -1);
                eqk->ix[nt] = (int)(+vym);
                eqk->cf[nt] = (int)(-wym);
                eqk->rhs += +wym*dym;
                nt++;
              }
            if (has_yp)
              { long int vyp;
                double dyp, wyp;
                pst_graph_vertex_get_neighbour(g, i, 0, 1, &vyp, &dyp, &wyp);
                assert(vyp != -1);
                eqk->ix[nt] = (int)(+vyp);
                eqk->cf[nt] = (int)(-wyp);
                eqk->rhs += +wyp*dyp;
                nt++;
              }
            double wtot = 0;
            int j;
            for (j = 1; j < nt; j++) { assert(eqk->cf[j] < 0); wtot+= -eqk->cf[j]; } 
            assert(wtot > 0);
            for (j = 1; j < nt; j++) { eqk->cf[j]/=wtot;  }
            eqk->rhs /= wtot;
            eqk->nt = nt;
            N++;
          }
        else
          { ind_ix[i] = -1; }
      }

    long int k;
    for (k = 0; k < N; k++)
      { pst_imgsys_equation_t *eqk = &(eq[k]);
        int nt = eqk->nt;
        int mt = 0;
        int i;
        for (i = 0; i < nt; i++)
          { /* Get the temporay index {xyi}: */
            int xyi = eqk->ix[i];
            /* Get the definitive index {ki}: */
            int ki = (int)(ind_ix[xyi]);
            if (ki >= 0)
              { /* Append the term to the equation: */
                int j = mt;
                eqk->ix[j] = ki;
                eqk->cf[j] = eqk->cf[i];
                assert(!isnan(eqk->cf[j]));
                mt++;
              }
          }
        eqk->nt = mt;
      }

    /* Build the inverse tables {col[0..N-1], row[0..N-1]}: */
    int *col = (int *)notnull(malloc(N*sizeof(int)), "no mem");
    int *row = (int *)notnull(malloc(N*sizeof(int)), "no mem");
    int *ix = (int *)notnull(malloc((NXY_Z)*sizeof(int)), "no mem");
      { long int xy;
        for (xy = 0; xy < NXY_Z; xy++) { ix[xy] = -1; }
        long int count_idx = 0;
        for (xy = 0; xy <g->num_vertex; xy++) 
          { if (ind_ix[xy] >= 0)
              { long int k = ind_ix[xy];
                long int x, y;
                pst_graph_restore_vertex_index(g->vertices[xy].id, NX_Z, NY_Z, &x, &y);
                col[count_idx] = (int)x;
                row[count_idx] = (int)y;
                assert((k >> 30) == 0); /* Should never happen. */
                ix[g->vertices[xy].id] = (int)k;
                count_idx++;
              }
          }
        assert(count_idx == N);
      }

    free(ind_ix);
    /* Now package the equations as a system: */
    pst_imgsys_t *S = pst_imgsys_from_eqs((int)NX_Z, (int)NY_Z, N, eq, ix, col, row);  
    return S;
  }

void pst_graph_solve_system
  ( pst_graph_t *g,
    pst_imgsys_t *S,
    float_image_t *OZ, 
    long int maxIter, 
    double convTol, 
    int para, 
    int szero, 
    bool_t verbose 
  )
  { double *Z = rn_alloc(S->N);
    pst_slope_map_copy_height_map_to_sol_vec(S, OZ, Z);
    pst_imgsys_solve(S, Z, NULL, (int)maxIter, convTol, para, szero, verbose, 0, NULL);
    pst_slope_map_copy_sol_vec_to_height_map(S, Z, OZ);
    free(Z);
  }

void pst_graph_estimate_from_shrunk(pst_graph_t *g, pst_graph_t *jg, float_image_t *OZ)
  {
    int i;
    int ind_vt = 0;
    long int NX_Z = (int)OZ->sz[1];
    long int NY_Z = (int)OZ->sz[2];
    for (i = 0; i < g->num_vertex; i++)
    {
      pst_vertex_t *v = &(g->vertices[i]);
      if (((v->i + v->j)%2) != 0) { continue; }
      while (jg->vertices[ind_vt].id != v->id)
        { ind_vt++;
          assert(ind_vt < jg->num_vertex);
        }
      long int x, y;
      pst_graph_restore_vertex_index(v->id, NX_Z, NY_Z, &x, &y);
      assert((x != -1) && (y != -1));
      double z = float_image_get_sample(OZ, 0, (int)x, (int)y);
      int has_xp = (v->vtxp != -1); 
      int has_xm = (v->vtxm != -1);
      int has_yp = (v->vtyp != -1);
      int has_ym = (v->vtym != -1);
      if (has_xp)
          { long int vxp;
            double dxp, wxp;
            pst_graph_vertex_get_neighbour(g, i, 1, 0, &vxp, &dxp, &wxp);
            assert(vxp != -1);
            long int ix, iy;
            pst_graph_restore_vertex_index(g->vertices[vxp].id, NX_Z, NY_Z, &ix, &iy);
            float_image_set_sample(OZ, 0, (int)ix, (int)iy, (float)(z - dxp));
          }
        if (has_yp)
          { long int vyp;
            double dyp, wyp;
            pst_graph_vertex_get_neighbour(g, i, 0, 1, &vyp, &dyp, &wyp);
            assert(vyp != -1);
            long int ix, iy;
            pst_graph_restore_vertex_index(g->vertices[vyp].id, NX_Z, NY_Z, &ix, &iy);
            float_image_set_sample(OZ, 0, (int)ix, (int)iy, (float)(z - dyp));
          }
        if (has_xm)
          { long int vxm;
            double dxm, wxm;
            pst_graph_vertex_get_neighbour(g, i, -1, 0, &vxm, &dxm, &wxm);
            assert(vxm != -1);
            long int ix, iy;
            pst_graph_restore_vertex_index(g->vertices[vxm].id, NX_Z, NY_Z, &ix, &iy);
            float_image_set_sample(OZ, 0, (int)ix, (int)iy, (float)(z - dxm));
          }
        if (has_ym)
          { long int vym;
            double dym, wym;
            pst_graph_vertex_get_neighbour(g, i, 0, -1, &vym, &dym, &wym);
            assert(vym != -1);
            long int ix, iy;
            pst_graph_restore_vertex_index(g->vertices[vym].id, NX_Z, NY_Z, &ix, &iy);
            float_image_set_sample(OZ, 0, (int)ix, (int)iy, (float)(z - dym));
          }
      }
  }

void pst_graph_integration_recursive
  ( pst_graph_t *g,
    float_image_t *OZ,
    long int maxIter,
    double convTol, 
    int para, 
    int szero, 
    bool_t verbose,
    int level
  )
  { char *filename_grph = NULL;
    char *filename_grph = jsprintf("graph-%02d.txt", level);
    FILE *wr_grph = open_write(filename_grph, FALSE);
    pst_graph_write(wr_grph, g);
    fclose(wr_grph);

    if (g->num_vertex > 2)
      { fprintf(stderr, "reducing Graph level[%d] with [%ld] vertices\n", level, g->num_vertex);
        pst_graph_t *jg = pst_graph_shrink(g);
        pst_graph_integration_recursive(jg, OZ, (int)ceil(((double)maxIter)*M_SQRT2), convTol/M_SQRT2, para, szero, verbose, level+1);
        pst_graph_estimate_from_shrunk(g, jg, OZ);

        char *filename = jsprintf("guess-%02d.fni", level);
        FILE *wr_dump = open_write(filename, FALSE);
        float_image_write(wr_dump, OZ);
        fclose(wr_dump);
        pst_graph_free(jg);
      }
    else
      { fprintf(stderr, "end of Recursion - returning\n");
        float_image_fill_channel(OZ, 0, 0);
      }
    fprintf(stderr, "solving Graph level[%d] with [%ld] vertices\n", level, g->num_vertex);
    pst_imgsys_t *S = pst_graph_build_integration_system(g, (int)OZ->sz[1], (int)OZ->sz[2]);

    char *filename = jsprintf("system-%02d.txt", level);
    FILE *wr_dump = open_write(filename,  FALSE);
    pst_imgsys_write(wr_dump, S);
    fclose(wr_dump);

    filename = NULL;
    char *filename = jsprintf("height-%02d.fni", level);
    wr_dump = open_write(filename, FALSE);
    float_image_write(wr_dump, OZ);
    fclose(wr_dump);

    pst_graph_solve_system(g, S, OZ, maxIter, convTol, para, szero, verbose);
    pst_imgsys_free(S);
  }

void pst_graph_free(pst_graph_t *g)
  { free(g->vertices);
    free(g->edges);
    free(g);
  }
