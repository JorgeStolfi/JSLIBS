/* uint16_color_tree.h - choose a representative set of N colors, from color histogram */
/* Last edited on 2024-12-26 20:20:49 by stolfi */ 

#ifndef uint16_color_tree_H
#define uint16_color_tree_H

#include <stdint.h>

#include <uint16_image.h>
#include <uint16_color_table.h>

#define MAXCOLORS 32767
    
typedef struct uint16_color_tree_node_t
  { uint16_t *clr;       /* Root color. */
    int32_t csplit;      /* Channel of splitting. */
    uint16_t lim[2];     /* Extremal sample values in the subtrees. */
    struct uint16_color_tree_t *sub[2]; /* Subtree roots, or {NULL}. */
  } uint16_color_tree_node_t;
  /* A record {nd} of this type is root node of a subtree of a tree of color values. 
  
    The subtree contains at least one color, {nd.clr_root[0..chns-1]}.
    The other nodes in the subtree are partitioned into two
    sub-subtrees, {sub[0]} and {sub[1]} according to their value in
    channel {csplit}. Either subtre may be {NULL}.
    
    Specifically, any nodes {nd0} and {nd1} in {nd.sub[0]} and {nd.sub[1]} will 
    satisfy have {nd0.clr[csplit] <= nd.clr[csplit] <= nd1.clr[csplit]}.
    
    Moreover, if {nd.sub[0]} is not {NULL}, {nd.lim[0]} will be the maximm of {nd0.clr[csplit]},
    for all {nd0} in that sub-subtree.   Similarly, if {nd.sub[1]} is not {NULL},
    then {nd.lim[1]} will be the minimum of {nd1.clr[csplit]} for {nd1} in that
    subtree. */
    
typedef struct uint16_color_tree_t
  { uint32_t chns;   /* Number of channels. */
    uint16_t maxval; /* Max sample value in all nodes. */
    uint32_t NT;     /* Number of nodes in the tree. */
    uint16_color_tree_node_t *root; 
    uint16_t *min;   /* Lower bound of all colors, in each channel. */
    uint16_t *max;   /* Upper bound of all colors, in each channel. */
  } uint16_color_tree_t;
  /* Head node of a color search tree.  The {min} and {max} vectors
    are the low and high corners of the tight bounding box of all colors. */

uint16_color_tree_t *uint16_color_tree_build
  ( uint32_t NN,                       /* Number of input colors. */
    uint32_t chns,                     /* Number of channels. */
    uint16_t *smp,                     /* Lienarized color samples */
    uint16_t maxval,                   /* Max sample value. */
    uint32_t maxColors                 /* Max desired colors. */
  );
  /* Given a set of {NN} colors, chooses a subset of at most {maxColors} colors
    that is supposedly good for image dithering/quantization, and builds a
    search tree for them.
    
    Each color is a tuple of {chns} samples. The {NN} tuples are assumed to be all
    concatenated in a single vector {smp[0..NN*chsn-1]} of samples. Every
    sample must be in the range {0..maxval}.
    
    The procedure is based on Paul Heckbert's median cut algorithm, as described in
    "Color Image Quantization for Frame Buffer Display", SIGGRAPH '82
    Proceedings, page 297. */

typedef double uint16_color_tree_dist_func_t(uint32_t chns, uint16_t maxval0, int32_t clr0[], uint16_t maxval1, int32_t clr1[]);
  /* Type of a function that computes a distance between two color
    tuples {clr0[0..chns-1]} and {clr1[0..chns-1]}. 
    
    The function should assume that the `normal' range of samples is
    {0..maxval0}, and {0..maxval1}, respectively, but must accept
    elements in {clra} and {clrb} that are negative or greater than
    {maxval}, since these may be generated by error propagation
    algorithms like Floyd-Steinberg. */

uint16_color_table_node_t *uint16_color_table_closest_match
  ( uint32_t chns,                       /* Number of channels. */
    int32_t *clr,                        /* Color to search. */
    uint16_t maxval,                     /* Max sample value. */
    uint16_color_tree_t *tr,             /* Tree of colors to search. */
    uint16_color_tree_dist_func_t *dist  /* Color distance function. */
  );
  /* Finds the node of tree {tr} whose color is closest to {clr[0..chns-1]}
    by the distance {dist}.
    
    The procedure assumes that the nominal range of the samples in {clr}
    is {0..maxval}, but the samples may actually be negative or gerater
    than {maxval}.  Note that the {maxval} may be different from {tr.maxval},
    but {chns} must be equal to {tr.chns}. */

#endif