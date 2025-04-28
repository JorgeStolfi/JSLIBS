/* Image Foresting Transform (IFT) - Implementation */
/* Last edited on 2025-04-24 14:14:27 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <ift.h>
#include <pqueue.h>

/* Internal procedures: */

void ift_make_arcs(double radius, uint32_t cols, uint32_t rows, ift_rel_arc_t **arcp, uint32_t *arcsp);
  /* 
    Computes a relative arc table for an Euclidean neighborhood
    of the given radius (minimum 1.0). In particular, 
      radius = 1.0 is the 4-neighbor topology, 
      radius = 1.5 is the 8-neighbor topology.
    Returns the arcs in {*arcp} and their number in {*arcsp}. */

void ift_debug_node(char *op, ift_node_t *s, ift_path_cost_t C_new);
  /* 
    Prints data about the node {s} (if not null) and the cost {C_new} (if not NAN), prefixed by {op}. */
    
/* #define DEBUG_NODE(Op,Nd,Ct) ift_debug_node((Op),(Nd),(Ct)) */
#define DEBUG_NODE(Op,Nd,Ct) {}

/* IMPLEMENTATIONS */

ift_path_cost_t ift_check_path_cost(ift_path_cost_t tC, ift_path_cost_t sC);
  /*
    Checks whether the path cost {tC} is an integer and 
    greater (but not excessively) than the cost {sC} of the 
    previous node. */

void ift_initialize_forest(ift_graph_t *G);
  /* 
    Makes every pixel {p} into a trivial tree: sets the tree root {p->R}
    to {p} itself and the predecessor {p->P} to {NULL}. */
    
/* IMPLEMENTATIONS */

ift_graph_t *ift_make_graph(uint32_t cols, uint32_t rows, double radius)
  {
    uint32_t nodes = cols*rows;
    ift_graph_t *G;
    demand(cols <= ift_MAX_COLS, "too many cols");
    demand(rows <= ift_MAX_ROWS, "too many rows");
    demand(nodes <= ift_MAX_NODES, "too many nodes");
    G = (ift_graph_t *)notnull(malloc(sizeof(ift_graph_t)), "out of memory");
    G->cols = cols;
    G->rows = rows;
    /* Allocate and initialize the ift_node_t vector: */
    G->nodes = nodes;
    G->node = talloc(nodes, ift_node_t);
    int32_t i;
    i = 0;
    for (ift_pixel_index_t row = 0; row < rows; row++)
      for (ift_pixel_index_t col = 0; col < cols; col++)
        { ift_node_t *pg = &(G->node[i]);
          pg->col = col; pg->row = row;
          /* Make {pg} into a trivial tree: */
          pg->P = NULL; pg->R = pg;
          /* Set cost to +oo: */
          pg->C = +INFINITY;
          i++;
        } 
    /* Allocate and initialize the ift_rel_arc_t vector: */
    ift_make_arcs(radius, cols, rows, &(G->arc), &(G->arcs));
    return G;
  }

void ift_make_arcs(double radius, uint32_t cols, uint32_t rows, ift_rel_arc_t **arcp, uint32_t *arcsp)
  {
    /* Compute the max col/row displacement {m} of an arc: */
    double max_radius = hypot(cols, rows); 
    int32_t m = (int32_t)floor(fabs(fmin(radius, max_radius)));
    
    double r2 = radius*radius;
    ift_rel_arc_t *arc;
    uint32_t arcs;
    /* Count the number of arcs in the neighborhood (excluding the self-loop): */
    arcs = 0;
    for (int32_t drow = -m; drow <= m; drow++)
      for (int32_t dcol = -m; dcol <= m; dcol++)
        { double d2 = (double)(dcol*dcol + drow*drow);
          if ((d2 != 0) && (d2 <= r2)) {arcs++; }
        }
        
    /* Now allocate the relative arc array and fill it: */
    arc = talloc(arcs, ift_rel_arc_t);
    arcs = 0;
    for (int32_t drow = -m; drow <= m; drow++)
      for (int32_t dcol = -m; dcol <= m; dcol++)
        { double d2 = (double)(dcol*dcol + drow*drow);
          if ((d2 != 0) && (d2 <= r2))
            { ift_rel_arc_t *pa = &(arc[arcs]);
              pa->dcol = (ift_pixel_step_t)dcol; 
              pa->drow = (ift_pixel_step_t)drow;
              pa->daddr = drow*(int32_t)cols + dcol;
              pa->len = sqrt(d2);
              arcs++;
            }
        }
    (*arcp) = arc;
    (*arcsp) = arcs;
  }

ift_node_index_t ift_node_index(ift_graph_t *G, ift_pixel_index_t col, ift_pixel_index_t row)
  {
    demand((col >= 0) && (col < G->cols), "bad col");
    demand((row >= 0) && (row < G->rows), "bad row");
    return ((int32_t)G->cols)*row + col;
  }

ift_node_t *ift_get_node(ift_graph_t *G, ift_pixel_index_t col, ift_pixel_index_t row)
  {
    demand((col >= 0) && (col < G->cols), "bad col");
    demand((row >= 0) && (row < G->rows), "bad row");
    return &(G->node[((int32_t)G->cols)*row + col]);
  }

void ift_compute_forest(
    ift_graph_t *G, 
    ift_path_cost_fn_t *pf, 
    ift_tie_breaking_t tbreak, 
    ift_scan_order_t order,
    ift_path_cost_t *maxCostp
  )
  {
    pqueue_t *Q = pqueue_new();
    pqueue_realloc(Q, (uint32_t)G->nodes, (uint32_t)G->nodes);
    
    /* Insert all seed nodes in the queue, with their trivial path costs; */
    /* Leave all other nodes with infinite cost. */
    for (uint32_t i = 0; i < G->nodes; i++)
      { assert(G->rows*G->cols == G->nodes);
        uint32_t j = (order == ift_order_UP ? i : (uint32_t)(G->nodes-1-i));
        ift_node_t *pg = &(G->node[j]);
        /* Make {pg} into a trivial forest: */
        pg->R = pg;
        pg->P = NULL;
        assert(G->rows*G->cols == G->nodes);
        /* Initialize the node cost to the cost of the trivial path {(pg)}: */
        pg->C = ift_check_path_cost(pf(NULL, NULL, pg, G), 0.0);
        DEBUG_NODE("  ini  ", pg, NAN);
        pqueue_insert(Q, (uint32_t)j, pg->C);
        assert(G->rows*G->cols == G->nodes);
     }

    /* Main loop */
    *(maxCostp) = 0;
    while (pqueue_count(Q) > 0)
      { pqueue_item_t js = pqueue_head(Q);
        ift_node_t *s = &(G->node[js]);
        assert(pqueue_value(Q, js) == s->C);
        pqueue_delete(Q, js);
        DEBUG_NODE("  pop  ", s, NAN);
        assert(s->C >= *(maxCostp));
        if (isfinite(s->C)) { *(maxCostp) = s->C; }
        for (uint32_t ia = 0; ia < G->arcs; ia++)
          { uint32_t ja = (order == ift_order_UP ? ia : (uint32_t)(G->arcs-1-ia));
            ift_rel_arc_t *a = &(G->arc[ja]);
            int32_t tcol = s->col + a->dcol;
            int32_t trow = s->row + a->drow;
            if ((tcol >= 0) && (tcol < G->cols) && (trow >= 0) && (trow < G->rows))
              { pqueue_item_t jt = js + (uint32_t)a->daddr;
                ift_node_t *t = s + a->daddr;
                if (pqueue_has(Q, jt))
                  { if ((t->C > s->C) || (tbreak == ift_tbreak_LIFO))
                      { ift_path_cost_t tC_new = ift_check_path_cost(pf(s, a, t, G), s->C);
                        if ((tC_new < t->C)  || 
                            ((tC_new == t->C) && (tbreak == ift_tbreak_LIFO)))
                          { DEBUG_NODE("    -  ",t, tC_new);
                            t->C = tC_new; t->P = s; t->R = s->R;
                            pqueue_set_value(Q, jt, t->C);
                          }
                        else
                          { DEBUG_NODE("    <  ",t, tC_new); }
                      }
                    else
                      { DEBUG_NODE("    << ",t, s->C); }
                  }
                else
                  { DEBUG_NODE("    old",t, NAN); }
              }
          }
      }
  }
  
void ift_debug_node(char *op, ift_node_t *s, ift_path_cost_t C_new)
  {
    fprintf(stderr, "  %-4s", op);
    if (s == NULL)
      { fprintf(stderr, "NULL"); }
    else
      { ift_path_cost_t sC = s->C;
        fprintf(stderr, " (%4d,%4d)", s->col, s->row);
        if (sC == +INFINITY) 
          { fprintf(stderr, " cost = +oo"); }
        else
          { fprintf(stderr, " cost = %12.7f", sC); }
      }
    if (! isnan(C_new))
      { if (C_new == +INFINITY) 
          { fprintf(stderr, " new cost = +oo"); }
        else
          { fprintf(stderr, " new cost = %12.7f", C_new); }
      }
    fprintf(stderr, "\n");
  }
  
ift_path_cost_t ift_check_path_cost(ift_path_cost_t tC, ift_path_cost_t sC)
  {
    demand(tC >= sC, "path cost decreased");
    return tC;
  }
