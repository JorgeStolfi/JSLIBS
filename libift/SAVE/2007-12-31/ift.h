/* The Image Foresting Transform algorithm */
/* Last edited on 2007-12-27 02:43:16 by stolfi */

#ifndef ift_H
#define ift_H

#define _GNU_SOURCE
#include <math.h>
#include <values.h>

/* Pixel indices and their increments: */
typedef unsigned short PixelIndex;
typedef signed short PixelStep;
#define MAX_COLS (32768)
#define MAX_ROWS (32768)
#define MAX_NODES (((long)MAX_ROWS)*((long)MAX_COLS))

/* Input pixel values: */
typedef double Intensity;
/* #define MAX_INTENSITY (1.0) */
#define MAX_CHANNELS 3
typedef struct { Intensity c[MAX_CHANNELS]; } PixelValue;
#define ZERO_NODE_VALUE (PixelValue){{0,0,0}}

/* Path costs */

/* 
  Path costs are integers, normalized so that the maximum increment in
  the cost of a path due to adding one arc is `MAX_FINITE_ARC_COST'.
  This constant is large enough to store exact L_1 pixel distances
  for both PPM and PGM images, but still allows a bucket-based queue.
  
  Since a single path may include all pixels in the image, even with
  this restriction a path cost may not fit in a `long' variable ---
  which is just 32 bits on most machines. Therefore, we store path
  costs in `double' variables, where integers up to 2^50 can be stored
  without any rounding. */
  
/* Increment in path cost for adding one arc. */
typedef double ArcCost;  /* Actually a (big) unsigned integer */ 
#define MAX_FINITE_ARC_COST (65535)
#define INFINITE_ARC_COST INFINITY

/* Discretized path costs: */
typedef double PathCost;  /* Actually a (big) unsigned integer */ 
#define MAX_FINITE_PATH_COST (((double)MAX_FINITE_ARC_COST)*((double)MAX_COLS)*((double)MAX_ROWS))
#define INFINITE_PATH_COST INFINITY

/* Propagated seed labels: */
typedef unsigned short SeedLabel;
#define MAX_LABEL 65535
#define NON_SEED_LABEL 0
#define DEFAULT_SEED_LABEL 1

/* Representation of a pixel in the image graph: */
typedef struct PixelNode *PixelPtr;
typedef struct PixelNode
  { PixelIndex col, row;  /* Pixel indices. */
    PixelValue y;         /* Image pixel value (if grayscale, only `y[0]' is used). */
    PathCost C;           /* Cost of a path to the node (handicap or minimum). */
    SeedLabel L;          /* Label of pixel, if seed, or `NON_SEED_LABEL'. */
    PixelPtr P;           /* Predecessor in path forest, or NULL if root. */
    PixelPtr R;           /* Root of path forest. */
    PixelPtr prev, next;  /* Prev and next node in bucket list. */
  } PixelNode;
  
/* Representation of a relative arc out of a generic pixel: */
#define MAX_ARC_LENGTH 127
typedef struct RelArc
  { PixelStep dcol, drow;  /* Increments in pixel indices, up to MAX_ARC_LENGTH. */
    long daddr;            /* Index increment in PixelNode vector. */
    double len;            /* Euclidean length of arc */
  } RelArc;
  
/* Graph representation of an image and its IFT: */
typedef struct  
  { int cols;         /* Number of columns in image. */
    int rows;         /* Number of rows in image. */
    int channels;     /* Number of channels actually used. */
    long nodes;       /* Number of nodes = cols*rows */
    PixelNode *node;  /* Nodes in row-by-row order. */
    int arcs;         /* Number of neighbors of a generic pixel. */
    RelArc *arc;      /* Arcs out of a generic pixel. */
  } ImageGraph;

ImageGraph *ift_make_graph(int cols, int rows, double radius);
  /* 
    Builds an image graph suitable for an image with the specified dimensions,
    and Euclidean adjacency relation with neighborhood of the given radius
    (minimum 1.0, maximum 127). In particular, 
      radius = 1.0 is the 4-neighbor topology, 
      radius = 1.5 is the 8-neighbor topology.
    Also sets fields `col,row' of all pixels to the proper indices,
    and field `y' to zeros. */
    
PixelNode *ift_get_pixel_node(ImageGraph *G, PixelIndex col, PixelIndex row);
  /* 
    Returns the address of the node of `G' that represents pixel `[col,row]'
    of the image. */
    
void ift_set_all_seed_labels(ImageGraph *G, SeedLabel L);
  /* 
    Sets the label of every pixel to `L'. If `L == NON_SEED_LABEL', no pixel
    is a seed; otherwise everybody is a seed. */

void ift_set_seed_label(ImageGraph *G, PixelIndex col, PixelIndex row, SeedLabel L);
  /* 
    Sets the label of pixel `[col,row]' to `L'.
    If `L' is NON_SEED_LABEL it means that the pixel is *not* a seed. */
    
typedef PathCost PathCostFn (PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage);
  /* 
    A function to be called by `ift_compute_forest' to determine
    the cost of paths.
    
    When `stage = 1' the function should return the cost of the
    current path to `s' (as defined by the `P' fields), concatenated
    with the arc `a' that goes from `s' to `t'.  If `s' and/or `a'
    are NULL, the procedure should compute the cost of the trivial 
    path `<t>'.
    
    The cost should be normalized so that the increment for adding one
    more edge (or for a trivial path) is at most MAX_FINITE_ARC_COST,
    and rounded to an integer value.
    
    Upon entry to `ift_compute_forest', the cost function will be
    called once with `stage = 0,' to give it a chance to preprocess
    the graph `G' and/or initialize its internal data structures. The
    function will be called again at the end, with `stage = 2', to
    give it a chance to release its resources. */

typedef enum {ift_tbreak_LIFO, ift_tbreak_FIFO} ift_tie_breaking;
typedef enum {ift_order_UP, ift_order_DN} ift_scan_order;

void ift_compute_forest(
    ImageGraph *G, 
    PathCostFn pf, 
    ift_tie_breaking tbreak, 
    ift_scan_order order,
    PathCost *maxCostp
  );
  /* 
    Computes the IFT of the image graph `G'. Assumes that the nodes
    have been initialized with `ift_initialize' so that each node `t'
    has been turned into a trivial path tree, with `t->R = t', `t->L
    = NON_SEED_LABEL', `t->P = NULL'.  Then some `L' fields may have
    been set with `ift_set_seed' or `ift_set_seeds_from_image'.
    
    The procedure begins by initializing the cost `C(t)' of every
    pixel `t'. A pixel `t' which has `L != NON_SEED_LABEL' is
    considered to be a seed; its cost `C' is initialized with the
    handicap cost, as computed by `pf(NULL, NULL, t, G, 1)'. If `L ==
    NON_SEED_LABEL', the pixel is *not* a seed, and its cost `C' is
    initialized with `INFINITE_PATH_COST'.
    
    The procedure then propagates the fields `C', `P', `R', and `L'
    according to the IFT algorithm.  The path costs are computed
    by `pf(s, a, t, G, 1)' where `t' is the final node of the path,
    `s' is the next-to-last node, and `a' is the arc `<s,t>'.
    The `pf' function should return the cost of the path `P*(s)·<s,t>'.
    
    The `tbreak' parameter specifies the tie-breaking policy for selecting
    the next pixel from the queue among several minimum-cost ones
    (FIFO = most ancient, LIFO = most recent). It also specifies whether
    the predecessor map `P' should be updated when a new path is found 
    to have the same cost as the current one (FIFO = no, LIFO = yes).  
    Finally, the `tbreak' parameter also affects the handling of 
    infinite-cost nodes (FIFO = each is left as a singleton tree, LIFO = 
    they are joined into a depth-first spanning forest).
    
    The `order' parameter specifies the order in which seeds
    and arcs are to be considered. Returns in `*maxCostp' the maximum
    finite path cost seen. */

#define IFT_ERROR(msg) \
  do { ift_error(__FILE__, __LINE__, msg); } while (0)

void ift_error(char *file, int line, char *msg);
  /* Prints `msg' and halts with error status. */

#endif
