/* Image Foresting Transform (IFT) - Implementation */
/* Last edited on 2007-12-26 20:34:00 by stolfi */

#include "ift.h"
#include "ift_queue.h"
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <stdio.h>

/* Internal procedures: */

void ift_make_arcs(double radius, int cols, int rows, RelArc **arcp, int *arcsp);
  /* 
    Computes a relative arc table for an Euclidean neighborhood
    of the given radius (minimum 1.0). In particular, 
      radius = 1.0 is the 4-neighbor topology, 
      radius = 1.5 is the 8-neighbor topology.
    Returns the arcs in `*arcp' and their number in `*arcsp'. */

void ift_debug_node(char *op, PixelNode *s);
  /* 
    Prints data about the node `s', prefixed by `op'. */
    
/* #define DEBUG_NODE(Op,Nd) ift_debug_node((Op),(Nd)) */
#define DEBUG_NODE(Op,Nd) {}

/* IMPLEMENTATIONS */

PathCost ift_check_path_cost(PathCost tC, PathCost sC);
  /*
    Checks whether the path cost `tC' is an integer and 
    greater (but not excessively) than the cost `sC' of the 
    previous node. */

void ift_initialize_forest(ImageGraph *G);
  /* 
    Makes every pixel `p' into a trivial tree: sets the tree root `p->R'
    to `p' itself and the predecessor `p->P' to `NULL'. */
    
/* IMPLEMENTATIONS */

ImageGraph *ift_make_graph(int cols, int rows, double radius)
  {
    long nodes = ((long)cols)*((long)rows);
    int col, row, k; long i;
    ImageGraph *G;
    PixelValue zero;
    if (cols > MAX_COLS) { IFT_ERROR("too many cols"); }
    if (rows > MAX_ROWS) { IFT_ERROR("too many rows"); }
    if (nodes > MAX_NODES) { IFT_ERROR("too many nodes"); }
    G = (ImageGraph *)malloc(sizeof(ImageGraph));
    if (G == NULL) { IFT_ERROR("out of memory"); }
    G->cols = cols;
    G->rows = rows;
    G->channels = MAX_CHANNELS;
    /* Define the null pixel value: */
    for (k = 0; k < MAX_CHANNELS; k++) { zero.c[k] = 0; }
    /* Allocate and initialize the PixelNode vector: */
    G->nodes = nodes;
    G->node = (PixelNode *)malloc(nodes*sizeof(PixelNode));
    i = 0;
    for (row = 0; row < rows; row++)
      for (col = 0; col < cols; col++)
        { PixelNode *pg = &(G->node[i]);
          pg->col = col; pg->row = row;
          pg->y = zero;
          /* Make `pg' into a trivial tree: */
          pg->P = NULL; pg->R = pg;
          /* Clear out the seed label, set cost to +oo: */
          pg->L = NON_SEED_LABEL; pg->C = INFINITE_PATH_COST;
          i++;
        } 
    /* Allocate and initialize the RelArc vector: */
    ift_make_arcs(radius, cols, rows, &(G->arc), &(G->arcs));
    return G;
  }

void ift_make_arcs(double radius, int cols, int rows, RelArc **arcp, int *arcsp)
  /* 
    Computes a relative arc table for an Euclidean neighborhood
    of the given radius .
    Returns the arcs in `*arcp' and their number in `*arcsp'. */
  {
    int m = (int)floor(fabs(radius));
    double r2 = radius*radius;
    RelArc *arc;
    int dcol, drow, arcs;
    if (m > MAX_ARC_LENGTH) { IFT_ERROR("radius too big"); }
    /* Count the number of arcs in the neighborhood: */
    arcs = 0;
    for (drow = -m; drow <= m; drow++)
      for (dcol = -m; dcol <= m; dcol++)
        { double d2 = (double)(dcol*dcol + drow*drow);
          if ((d2 != 0) && (d2 <= r2)) {arcs++; }
        }
    /* Now allocate the array and fill it: */
    arc = (RelArc*)malloc(arcs*sizeof(RelArc));
    if (arc == NULL) { IFT_ERROR("out of memory"); }
    arcs = 0;
    for (drow = -m; drow <= m; drow++)
      for (dcol = -m; dcol <= m; dcol++)
        { double d2 = (double)(dcol*dcol + drow*drow);
          if ((d2 != 0) && (d2 <= r2))
            { RelArc *pa = &(arc[arcs]);
              pa->dcol = dcol; pa->drow = drow;
              pa->daddr = drow*cols + dcol;
              pa->len = sqrt(d2);
              arcs++;
            }
        }
    (*arcp) = arc;
    (*arcsp) = arcs;
  }

PixelNode *ift_get_pixel_node(ImageGraph *G, PixelIndex col, PixelIndex row)
  {
    if (col >= G->cols)  { IFT_ERROR("bad col"); }
    if (row >= G->rows)  { IFT_ERROR("bad row"); }
    return &(G->node[G->cols*row + col]);
  }
    
void ift_set_all_seed_labels(ImageGraph *G, SeedLabel L)
  { long i;
    if (L > MAX_LABEL) { IFT_ERROR("bad L"); }
    for (i = 0; i < G->nodes; i++)
      { PixelNode *pg = &(G->node[i]); 
        pg->L = L;
      }
  }

void ift_set_seed_label(ImageGraph *G, PixelIndex col, PixelIndex row, SeedLabel L)
  { PixelNode *pg = ift_get_pixel_node(G, col, row);
    if (L > MAX_LABEL) { IFT_ERROR("bad L"); }
    pg->L = L;
  }

void ift_compute_forest(
    ImageGraph *G, 
    PathCostFn pf, 
    ift_tie_breaking tbreak, 
    ift_scan_order order,
    PathCost *maxCostp
  )
  {
    long i;
    int ia;
    ift_queue *Q = new_ift_queue(tbreak);
    
    /* Initialize the path cost function: */
    (void)pf(NULL, NULL, NULL, G, 0);
    
    /* Make all paths trivial. */
    /* Insert all seed nodes in the queue, with their handicap costs; */
    /* Leave all other nodes with infinite cost. */
    for (i = 0; i < G->nodes; i++)
      { long j = (order == ift_order_UP ? i : G->nodes-1-i);
        PixelNode *p = &(G->node[j]);
        /* Make `p' into a trivial forest: */
        p->R = p;
        p->P = NULL;
        /* Clear out the queue links, just in case: */
        p->prev = NULL;
        p->next = NULL;
        /* Initialize the node cost to the handicap cost if seed, +oo otherwise: */
        if (p->L == NON_SEED_LABEL) 
          { p->C = INFINITE_PATH_COST; }
        else
          { p->C = ift_check_path_cost(pf(NULL, NULL, p, G, 1), 0.0);}
        DEBUG_NODE("      +",p);
        ift_queue_insert(Q, p);
     }

    /* Main loop */
    *(maxCostp) = 0;
    while (! ift_queue_is_empty(Q))
      { PixelNode *s = ift_queue_delete_min(Q);
        DEBUG_NODE("  <    ", s);
        if (s->C < INFINITE_PATH_COST) { *(maxCostp) = s->C; }
        for (ia = 0; ia < G->arcs; ia++)
          { int ja = (order == ift_order_UP ? ia : G->arcs-1-ia);
            RelArc *a = &(G->arc[ja]);
            int tcol = s->col + a->dcol;
            int trow = s->row + a->drow;
            if ((tcol >= 0) && (tcol < G->cols) && (trow >= 0) && (trow < G->rows))
              { PixelNode *t = s + a->daddr;
                DEBUG_NODE("    ?  ",t);
                if (ift_queue_is_enqueued(t) &&
                    ((t->C > s->C) || (tbreak == ift_tbreak_LIFO)))
                  { PathCost tC_new = ift_check_path_cost(pf(s, a, t, G, 1), s->C);
                    if ((tC_new < t->C)  || 
                        ((tC_new == t->C) && (tbreak == ift_tbreak_LIFO)))
                      { DEBUG_NODE("      -",t); ift_queue_delete(Q, t);
                        t->C = tC_new; t->P = s; t->L = s->L; t->R = s->R;
                        DEBUG_NODE("      +",t); ift_queue_insert(Q, t);
                      }
                  }
              }
          }
      }
      
    /* Finalize the path cost function: */
    (void)pf(NULL, NULL, NULL, G, 2);
  }
  
void ift_debug_node(char *op, PixelNode *s)
  {
    fprintf(stderr, "  %-4s", op);
    if (s == NULL)
      { fprintf(stderr, "NULL"); }
    else
      { PathCost sC = s->C;
        int k;
        fprintf(stderr, " (%4d,%4d)", s->col, s->row);
        fprintf(stderr, " = (");
        for (k = 0; k < MAX_CHANNELS; k++)
          { fprintf(stderr, " %7.5f", s->y.c[k]); }
        fprintf(stderr, " )");
        if (sC == INFINITE_PATH_COST) 
          { fprintf(stderr, " cost = +oo"); }
        else
          { fprintf(stderr, " cost = %.0f (%g)", sC, sC/((double)MAX_FINITE_ARC_COST)); }
      }
    fprintf(stderr, "\n");
  }
  
PathCost ift_check_path_cost(PathCost tC, PathCost sC)
  {
    if (tC != INFINITE_PATH_COST) 
      { if (tC != floor(tC)) { IFT_ERROR("path cost is not an integer"); }
        if (tC < sC) { IFT_ERROR("path cost decreased"); }
        if (tC - sC > MAX_FINITE_ARC_COST) { IFT_ERROR("path cost increment too large"); }
      }
    return tC;
  }

void ift_error(char *file, int line, char *msg)
  { 
    fprintf(stderr, "%s:%d: %s\n", file, line, msg);
    exit(1);
  }
