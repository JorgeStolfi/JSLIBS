/* See {psgr_shrink.h} */
/* Last edited on 2025-04-27 10:34:27 by stolfi */

/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <r2.h>
#include <i2.h>

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_copy.h>
#include <psgr_shrink_stacy.h>

#include <psgr_shrink.h>

/* SHORTER NAMES */  

#define VNONE psgr_vertex_NONE
#define ENONE psgr_edge_NONE
#define ANONE psgr_arc_NONE

#define CLEAR   psgr_mark_CLEAR
#define KILL    psgr_mark_KILL   
#define KEEP    psgr_mark_KEEP   
#define DELETED psgr_mark_DELETED

psgr_arc_t psgr_find_rightmost_arc(psgr_t *gr, psgr_arc_t a);
  /* Finds the arc out of {a}'s origin that goes to the vertex with
    the highest {X} pixel index. */

void psgr_choose_vertices_for_removal(psgr_t* gr, uint32_t degree[], int32_t indent);
  /* Choses the vertices that are to be deleted by
    {psgr_shrink_decimate}. The vertices of {gr} must all be marked
    {CLEAR} or {DELETED}. Vertices that are marked as {DELETED} vertices
    are ignored. The table {degree[0..gr.NV-1]} must give the
    (out)degrees of the non-deleted vertices.
    
    The chosen vertices will be a maximal subset with degrees 0 to 6
    that are not neighbors of each other. The vertices to be deleted are
    marked as {KILL}, and their neighbors with {KEEP}
    
    If {indent} is non-negative, prints diagnostics indented by {indent} columns. */
  
void psgr_remove_paralel_edges(psgr_t* gr, int32_t indent);
  /* Locates parallel edges in the graph and merges them together.
    Namely, for every maximal set of edges with same origin and
    destination, combines all the {.data} and {.weight} fields into one
    of those edges, by the parallel formula, and marks all the other
    edges of the set as {DELETED}.
    
    If {indent} is non-negative, prints diagnostics indented by {indent} columns. */

void psgr_shrink_remove_vertex_deg_2(psgr_t *gr, psgr_vertex_t v, int32_t indent);
  /* The vertex {v} must have outdegree 2.  The procedure 
    removes the two edges out of {v} and the vertex
    {v} by a single edge connecting the other endpoints {v0} and {v1} of both edges.

    Specifically, the vertex {v} and its two edges of are disconnected
    from the graph and marked {DELETED}. Then delta of the new edge is
    the sum of their deltas, and the weight is the harmonic sum of their
    weights. The path of the new edge will be
    {psgr_path_concat(P0,ip,P1}) where {P0} and {P1} are the paths of
    the two edges, both oriented from {v0} to {v1}, and {ip} is the pix
    index pair of {v}.
    
    If {indent} is non-negative, prints diagnostics indented by {indent} columns. */

void psgr_shrink_remove_vertex_general(psgr_t* gr, psgr_vertex_t v, uint32_t deg, int32_t indent);
  /* Removes a vertex {v} with degree 3 to 6, performing the star-cycle 
    transformation. */

void psgr_shrink_debug_edge(psgr_t *gr, char *class, int32_t k, psgr_arc_t a, int32_t indent); 
  /* Prints data about arc {a} which should be the member of index {k} in  a set of arcs identified 
    by {class}, such as the onext ring of a vertex that is scheduled for deletion, or the 
    bridge or cycle arcs that have been added.
    
    In the first case, it should preferrably be called before any of the arcs
    in that ring are deleted. In the second case, it should preferrably be called 
    after all new edges have been added.  
    
    The index {k} printed only if it is non-negative. */

void psgr_vertex_star_remove(psgr_t *gr, psgr_vertex_t v, int32_t indent);
  
/* IMPLEMENTATIONS */

psgr_t *psgr_shrink(psgr_t *gri, int32_t vj_from_vi[], int32_t ej_from_ei[], int32_t indent)
  { psgr_t *grt = psgr_copy(gri, NULL, NULL);
    psgr_shrink_decimate(grt, indent);
    psgr_t *grj = psgr_copy(grt, vj_from_vi, ej_from_ei);
    psgr_free(grt);
    return grj;
  }

void psgr_shrink_decimate(psgr_t *gr, int32_t indent)
  {
    uint32_t NV = psgr_num_vertices(gr);
    
    uint32_t *degree = talloc(NV, uint32_t);
    for (uint32_t v_ix = 0; v_ix < NV; v_ix++)
      { psgr_vertex_t v = (psgr_vertex_t){ v_ix };
        demand(psgr_vertex_mark(gr, v) == psgr_mark_CLEAR, "vertex marks are not clean");
        degree[v_ix] = psgr_vertex_outdegree(gr, v);
      }
      
    psgr_choose_vertices_for_removal(gr, degree, indent);

    for (uint32_t v_ix = 0; v_ix < NV; v_ix++)
      { psgr_vertex_t v = (psgr_vertex_t){ v_ix };
        switch (psgr_vertex_mark(gr, v))
          { case CLEAR:
            case KEEP:
              { psgr_vertex_set_mark(gr, v, CLEAR); break; }
            case KILL:
              { psgr_vertex_star_remove(gr, v, indent);
                psgr_vertex_set_mark(gr, v, DELETED); break;
              }
            case DELETED:
              { affirm(FALSE, "unexpected deleted vertex"); }
            default:
              { affirm(FALSE, "invalid vertex mark"); }
          }
      }
    psgr_remove_paralel_edges(gr, indent);
    free(degree);
  }

void psgr_choose_vertices_for_removal(psgr_t *gr, uint32_t degree[], int32_t indent)
  {
    uint32_t NV = psgr_num_vertices(gr);

    uint32_t DEG_MAX = 6;  /* Mex degree of removed vertices. */
    uint32_t stats[DEG_MAX+1];  /* {stats[deg]} counts removed vertices of degree [deg}. */ 
    for (uint32_t deg = 0; deg <= DEG_MAX; deg++)
      { /* Mark vertices of degree {deg} for removal, keep neighbots: */
        stats[deg] = 0;
        for (uint32_t v_ix = 0; v_ix < NV; v_ix++)
          { psgr_vertex_t v = (psgr_vertex_t){ v_ix };
            if (psgr_vertex_mark(gr, v) != CLEAR) { continue; }
            uint32_t deg_v = psgr_vertex_outdegree(gr, v);
            if (deg_v == deg)
              { stats[deg]++;
                psgr_vertex_set_mark(gr, v, KILL);
                if (deg > 0)
                  { psgr_arc_t a = psgr_vertex_out_arc(gr, v);
                    psgr_arc_t b = a;
                    do 
                      { psgr_vertex_t u = psgr_arc_dst(gr, b);
                        demand(psgr_vertex_mark(gr, u) != DELETED, "neighbor is deleted");
                        demand(psgr_vertex_mark(gr, u) != KILL, "neighbor of killed is not kept");
                        psgr_vertex_set_mark(gr, u, KEEP);
                        b = psgr_arc_onext(gr, b);
                      } while(b.ix != a.ix);
                  }
              }
          }
       if (indent >= 0) 
         { fprintf(stderr, "%*smarked %d vertices with degree %d for removal\n", indent, "", stats[deg], deg); }
    }
  }

void psgr_remove_paralel_edges(psgr_t *gr, int32_t indent)
  {
    uint32_t NE = psgr_num_edges(gr);
    
    for (uint32_t e_ix = 0; e_ix < NE ; e_ix++)
      { psgr_edge_t e = (psgr_edge_t ){ e_ix };
        if (psgr_edge_mark(gr, e) == DELETED) { continue; }
         
        psgr_arc_t a = psgr_orient_edge(e, 0);
        assert(a.ix != ANONE.ix);
        psgr_arc_t b = psgr_arc_onext(gr, a);
        while (b.ix != a.ix)
          { psgr_arc_t bnext = psgr_arc_onext(gr, b);
            bool_t same_org = (psgr_arc_org(gr,a).ix == psgr_arc_org(gr,b).ix);
            bool_t same_dst = (psgr_arc_dst(gr,a).ix == psgr_arc_dst(gr,b).ix);
            if (same_org && same_dst)
              { /* Edges {a} and {b} are parallel: */

                /* fprintf(stderr, "%*sEdge %ld (%ld - %ld) ", indent,"", ind_a, org_a, dst_a); */
                /* fprintf(stderr, "%*sEdge %ld (%ld - %ld) ", indent,"", ind_b, org_b, dst_b); */
                /* Eliminate {b} accumulating the weight on {a}: */
                double di = psgr_arc_delta(gr, a);
                double wi = psgr_arc_weight(gr, a);
                assert(wi > 0);

                double dk = psgr_arc_delta(gr, b);
                double wk = psgr_arc_weight(gr, b);
                assert(wk > 0);

                double d = (di*wi + dk*wk)/(wk+wi);
                double w = wi + wk;

                psgr_arc_set_delta(gr, a, d);
                psgr_arc_set_weight(gr, a, w);
                psgr_edge_remove(gr, psgr_arc_edge(b), indent);
              }
            b = bnext;
          }
      }
  }
    
void psgr_vertex_star_remove(psgr_t *gr, psgr_vertex_t v, int32_t indent)
  { uint32_t NV = psgr_num_vertices(gr);

    demand((v.ix != VNONE.ix) && (v.ix < NV), "invalid vertex {v}");
    demand(psgr_vertex_mark(gr, v) == KILL, "trying to remove a vertex not marked {KILL}");
    uint32_t deg_v = psgr_vertex_outdegree(gr, v);
    psgr_arc_t a = psgr_vertex_out_arc(gr, v);
    if ( a.ix != ANONE.ix) { a = psgr_find_rightmost_arc(gr, a); }
    switch (deg_v)
      { case 0: 
          { psgr_vertex_remove(gr, v, indent);
            break;
          }
        case 1:
          { psgr_edge_remove(gr, psgr_arc_edge(a), indent); 
            psgr_vertex_remove(gr, v, indent);
            break;
          }
        case 2:
          { psgr_shrink_remove_vertex_deg_2(gr, v, indent);
            break;
          }
        default:
          { psgr_shrink_remove_vertex_general(gr, v, deg_v, indent);
            break;
          }
      }
     assert(psgr_vertex_mark(gr, v) == DELETED);
     assert(psgr_vertex_out_arc(gr, v).ix == ANONE.ix);
  }
          
void psgr_shrink_remove_vertex_deg_2(psgr_t *gr, psgr_vertex_t v, int32_t indent)      
  { /* Get the two arcs {a0,a1} exiting {v}: */
    psgr_arc_t a0 = psgr_vertex_out_arc(gr, v);
    demand(a0.ix != ANONE.ix, "vertex has degree 0");
    psgr_arc_t a1 = psgr_arc_onext(gr, a0);
    demand(a0.ix != a1.ix, "vertex has degree 1");
    demand(psgr_arc_onext(gr, a1).ix == a0.ix, "vertex has degree 3 or more");
    
    /* Get the two vertices to be connected: */
    psgr_vertex_t v0 = psgr_arc_dst(gr, a0); 
    psgr_vertex_t v1 = psgr_arc_dst(gr, a1); 
    demand((v0.ix != v.ix) && (v1.ix != v.ix), "loop edge detected");
    demand(v0.ix != v1.ix, "parallel edges detected");
    
    /* Disconnect both edges from graph, mark them as deleted: */
    if (indent >= 0) { psgr_shrink_debug_edge(gr, "star", 0, a0, indent); }
    if (indent >= 0) { psgr_shrink_debug_edge(gr, "star", 1, a1, indent); }
    psgr_edge_remove(gr, psgr_arc_edge(a0), indent);
    psgr_edge_remove(gr, psgr_arc_edge(a1), indent);

    /* Compute the weight and delta: */
    double delta = psgr_arc_delta(gr, a1) - psgr_arc_delta(gr, a0);
    demand(isfinite(delta), "delta overflow or undefined");
    double w0 = psgr_arc_weight(gr, a0); assert(isfinite(w0) && (w0 > 0));
    double w1 = psgr_arc_weight(gr, a1); assert(isfinite(w1) && (w1 > 0));
    double weight = 1.0/((1/w0) + (1/w1));
    assert(isfinite(weight) && (weight >= 0));
    if (weight <= 0)
      { if (indent >= 0) { fprintf(stderr, "%*s!! underflow in weight computation", indent,""); } }
    else
      { /* Create a path by joining the existing paths: */
        i2_t ip = psgr_vertex_pixel(gr, v);
        r2_t p = (r2_t){{ ip.c[0], ip.c[1] }};
        psgr_path_t P0 = psgr_path_reverse(psgr_arc_path(gr, a0));
        psgr_path_t P1 = psgr_arc_path(gr, a1);
        psgr_path_t P = psgr_path_concatenate(P0, &p, P1);

        /* Add the new edge from {v0} to {v1} with new {.delta}, {.weight}, {.path}: */
        psgr_arc_t an = psgr_add_edge(gr, v0, v1, delta, weight, P);
        if (indent >= 0) { psgr_shrink_debug_edge(gr, "bridge", -1, an, indent); }
      }
  }

void psgr_shrink_remove_vertex_general(psgr_t *gr, psgr_vertex_t v, uint32_t deg_v, int32_t indent)
  { assert(deg_v >=3);
    psgr_arc_t a = psgr_vertex_out_arc(gr, v);
    demand(a.ix != ANONE.ix, "vertex has degree 0");

    /* Compute deltas and weights for the srat-cycle transformation: */
    double wtot;
    double d_star[deg_v], w_star[deg_v];
    psgr_shrink_stacy_get_star_data(gr, a, deg_v, d_star, w_star, &wtot);
    if (indent >= 0) { fprintf(stderr, "%*s  wtot = %12.6lf  wtot^2 = %12.6lf \n", indent,"", wtot, wtot*wtot); }
    
    double w_cycle[deg_v], d_cycle[deg_v];
    psgr_shrink_stacy_compute_cycle_data(gr, deg_v, d_star, w_star, wtot, d_cycle, w_cycle);
    
    if (indent >= 0) 
      { /* Print the edges that are to be deleted: */
        psgr_arc_t b = a;
        for (int32_t k = 0; k < deg_v; k++)
          { psgr_shrink_debug_edge(gr, "star", k, b, indent);
            b = psgr_arc_onext(gr, b);
          }
      }
      
    /* Add the cycle edges: */
    psgr_arc_t an[deg_v];
    psgr_arc_t a0 = a;
    for (int32_t k = 0; k < deg_v; k++)
      { an[k] = ANONE;
        psgr_arc_t a1 = psgr_arc_onext(gr, a0);
        demand(a0.ix != a1.ix, "vertex has degree 1");
        if (w_cycle[k] == 0)
          { if (indent >= 0) { fprintf(stderr, "%*s!! underflow in weight computation", indent,""); } }
        else
          { /* Get the two vertices to be connected: */
            psgr_vertex_t v0 = psgr_arc_dst(gr, a0);
            psgr_vertex_t v1 = psgr_arc_dst(gr, a1);
            demand((v0.ix != v.ix) && (v1.ix != v.ix), "loop edge detected");
            demand(v0.ix != v1.ix, "parallel edges detected");

            /* Create a path just inside the wedge {v0--v--v1}: */
            i2_t ip = psgr_vertex_pixel(gr, v);
            i2_t ip0 = psgr_vertex_pixel(gr, v0);
            i2_t ip1 = psgr_vertex_pixel(gr, v1);
            psgr_path_t P = psgr_path_in_star_wedge(&ip0, &ip, &ip1);

            /* Adds a new cycle edge from {v0} to {v1}: */
            an[k] = psgr_add_edge(gr, v0, v1, d_cycle[k], w_cycle[k], P);

            assert(psgr_arc_lnext(gr, a0).ix == an[k].ix);
          }
        a0 = a1;
      }
      
    /* Remove the star edges: */
    a0 = a;
    for (int32_t k = 0; k < deg_v; k++)
      { psgr_arc_t a1 = psgr_arc_onext(gr, a0);
        psgr_edge_remove(gr, psgr_arc_edge(a0), indent);
        a0 = a1;
      }
    
    if (indent >= 0) 
      { /* Print the cycle edges: */
        for (int32_t k = 0; k < deg_v; k++)
          { psgr_shrink_debug_edge(gr, "cycle", k, an[k], indent); }
      }
  }

psgr_arc_t psgr_find_rightmost_arc(psgr_t *gr, psgr_arc_t a)
  { uint32_t NV = psgr_num_vertices(gr);
    uint32_t NE = psgr_num_edges(gr);
    demand(a.ix != ANONE.ix, "invalid null arc {a}");
    psgr_edge_t e = psgr_arc_edge(a);
    demand(e.ix < NE, "invalid arc {a}");
    psgr_vertex_t org = psgr_arc_org(gr, a);
    assert(org.ix < NV);
    psgr_arc_t a_max = ANONE;
    int32_t cx_max = INT32_MIN;
    psgr_arc_t b = a;
    do 
      { psgr_vertex_t dst = psgr_arc_org(gr,psgr_arc_sym(b));
        i2_t ip = psgr_vertex_pixel(gr, dst);
        int32_t cx = ip.c[0];
        if (cx > cx_max) { a_max = b; cx_max = cx; }
        b = psgr_arc_onext(gr, b);
      } while (b.ix != a.ix);
    return a_max;
  }

void psgr_shrink_debug_edge(psgr_t *gr, char *class, int32_t k, psgr_arc_t a, int32_t indent)
  {
    fprintf(stderr,"%*s  %s", indent,"", class);
    if (k >= 0) { fprintf(stderr, "[%d]", k); }
    fprintf(stderr, " ");
    psgr_arc_print(stderr, gr, a);
    psgr_vertex_t org = psgr_arc_org(gr, a);
    psgr_vertex_t dst = psgr_arc_org(gr,psgr_arc_sym(a));
    fprintf(stderr, " v%d --> v%d", org.ix, dst.ix);
    double   d   = psgr_arc_delta(gr, a);
    double   w   = psgr_arc_weight(gr, a);
    fprintf(stderr, " d = %+12.6f w: %8.6f", d,w);
    fprintf(stderr, "\n");
  }

