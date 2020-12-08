/* The Image Foresting Transform algorithm */
/* See Falcao et al., IEEE Trans. on Patt. Anal. and Mach. Intel. (TPAMI), 2004 */
/* Last edited on 2010-06-06 14:53:12 by stolfi */

#ifndef ift_H
#define ift_H

#define _GNU_SOURCE
#include <math.h>
#include <float.h>
#include <stdint.h>

/* Pixel column and row indices and their increments: */
typedef int16_t ift_pixel_index_t;
typedef int16_t ift_pixel_step_t;
#define ift_MAX_COLS (32768)
#define ift_MAX_ROWS (32768)

/* Graph node indices (pixel index pairs, linearized): */
typedef int32_t ift_node_index_t;
#define ift_MAX_NODES (((uint32_t)ift_MAX_ROWS)*((uint32_t)ift_MAX_COLS))

/* Path costs. */
typedef double ift_path_cost_t;

/* Representation of a pixel in the image graph: */
typedef struct ift_node_t
  { ift_pixel_index_t col, row; /* Column and row indices. */
    ift_path_cost_t C;          /* Cost of a path to the node. */
    struct ift_node_t *P;       /* Predecessor in path forest, or NULL if root. */
    struct ift_node_t *R;       /* Root of path forest. */
  } ift_node_t;
  
/* Representation of a relative arc out of a generic pixel: */
typedef struct ift_rel_arc_t
  { ift_pixel_step_t dcol, drow;  /* Increments in pixel indices, up to MAX_ARC_LENGTH. */
    int32_t daddr;         /* Index increment in ift_node_t vector. */
    double len;            /* Euclidean length of arc */
  } ift_rel_arc_t;
  
/* Graph representation of an image and its IFT: */
typedef struct ift_graph_t 
  { int cols;            /* Number of columns in image. */
    int rows;            /* Number of rows in image. */
    int32_t nodes;       /* Number of nodes = cols*rows */
    ift_node_t *node;    /* Nodes in row-by-row order. */
    int arcs;            /* Number of neighbors of a generic pixel. */
    ift_rel_arc_t *arc;  /* Relative arcs out of a generic pixel. */
  } ift_graph_t;
  /* A pixel in column {col} and row {row} of the image is represented 
    by the node {G.node[ip]} where {ip = col + row*G.rows}. */

ift_graph_t *ift_make_graph(int cols, int rows, double radius);
  /*  Builds a graph suitable for an image with the specified dimensions.
    The relative arcs are those of the Euclidean adjacency relation 
    with given neighborhood {radius} (minimum 1.0). In particular, 
      radius = 1.0 is the 4-neighbor topology, 
      radius = 1.5 is the 8-neighbor topology.
    Also sets the fields {col,row} of all pixels to the proper indices. */
    
ift_node_index_t ift_node_index(ift_graph_t *G, ift_pixel_index_t col, ift_pixel_index_t row);
  /* Returns the index of the node of {G} that represents pixel {[col,row]}
    of the image. */
    
ift_node_t *ift_get_node(ift_graph_t *G, ift_pixel_index_t col, ift_pixel_index_t row);
  /* Returns the address of the node of {G} that represents pixel {[col,row]}
    of the image. */
    
typedef ift_path_cost_t ift_path_cost_fn_t (ift_node_t *s, ift_rel_arc_t *a, ift_node_t *t, ift_graph_t *G);
  /* A function to be called by {ift_compute_forest} to determine
    the cost of paths.
    
    The function should return the cost of the current path to {s} (as
    defined by the {P} fields), concatenated with the arc {a} that
    goes from {s} to {t}. If {s} and/or {a} are NULL, the procedure
    should compute the cost of the trivial path {<t>}. */

typedef enum {ift_tbreak_LIFO, ift_tbreak_FIFO} ift_tie_breaking_t;
typedef enum {ift_order_UP, ift_order_DN} ift_scan_order_t;

void ift_compute_forest(
    ift_graph_t *G, 
    ift_path_cost_fn_t *pf, 
    ift_tie_breaking_t tbreak, 
    ift_scan_order_t order,
    ift_path_cost_t *maxCostp
  );
  /* Computes the IFT of the image graph {G}. Assumes that the nodes
    have been initialized with {ift_initialize} so that each node {pg}
    has been turned into a trivial path tree, with {pg.R = pg}, 
    {pr.P = NULL}, {pr.C = +oo}.
    
    The procedure begins by setting the cost {C(t)} of every
    pixel {t} to {pf(NULL, NULL, t, G, 1)}.
    
    The procedure then propagates the fields {C}, {P}, and {R}
    according to the IFT algorithm.  The path costs are computed
    by {pf(s, a, t, G, 1)} where {t} is the final node of the path,
    {s} is the next-to-last node, and {a} is the arc {<s,t>}.
    The {pf} function should return the cost of the path {P*(s)·<s,t>}.
    
    The {tbreak} parameter specifies the tie-breaking policy for selecting
    the next pixel from the queue among several minimum-cost ones
    (FIFO = most ancient, LIFO = most recent). It also specifies whether
    the predecessor map {P} should be updated when a new path is found 
    to have the same cost as the current one (FIFO = no, LIFO = yes).  
    Finally, the {tbreak} parameter also affects the handling of 
    infinite-cost nodes (FIFO = each is left as a singleton tree, LIFO = 
    they are joined into a depth-first spanning forest).
    
    The {order} parameter specifies the order in which seeds
    and arcs are to be considered.
    
    The procedure returns in {*maxCostp} the maximum finite path cost
    seen. */

#endif
