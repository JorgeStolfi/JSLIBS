/* A hash table of multi-channel color values. */
/* Last edited on 2024-12-26 20:07:26 by stolfi */ 

/* Inspires on Jef Poskanzer's ppmcmap.h */

#ifndef uint16_color_table_H
#define uint16_color_table_H

#include <stdint.h>

#include <bool.h>
#include <uint16_image.h>

/* The procedures in this interface create and manipulate a table {cht}
  of multi-channel color values, intended to construct accurate color histograms
  and speed up repeated searches for closest match to the same color. */

typedef struct uint16_color_table_t 
  { struct uint16_color_table_node_t **buckets;
    uint32_t chns;
    uint16_t maxval;
    uint32_t NN;
  } uint16_color_table_t; 
  /* The table is an array {cht.buckets[0..HASH_SIZE-1]} entries, each
    containing {NULL} or a pointer to a linked list of nodes, in
    arbitrary order.  The field {cht.NN} has the total number of nodes in
    those lists, and {NN_max} is some arbitrary limit to {NN}. 
    The field {cht.chns} has the number of samples in  each color.

    The index into {cht.buckets} is a hash function computed from the {R},
    {G}, and {B} values. The nodes in the color table are all distinct
    and their coordinates are all in {0..cht.maxval}. */

typedef struct uint16_color_table_node_t 
  { uint16_t *clr; 
    uint32_t count;
    int32_t index;
    struct uint16_color_table_node_t *next;
  } uint16_color_table_node_t; 
  /* A node {nd} in a color table, representing the color with samples
    {nd.clr[0..chns-1]}, where {chns} is specified in the table
    header. The field {nd.count} is a count of occurrences of the color,
    and an {nd.index} may be an index into some table. */

uint16_color_table_node_t *uint16_color_table_lookup(uint16_color_table_t *cht, uint32_t chns, uint16_t clr[], bool_t add);
  /* If the color {clr[0..chns-1]} is in the hash table {cht} returns a pointer
    to its node after incrementing its {count}. 
    
    If the color {clr[]} is not in the table, and {add} is true, adds to {cht} a new
    node for it, with {count=1} and {index=-1}. In either case, the
    node is moved one (*only* one) step up in its bucket list if its new count is
    greater than that of the previous node.
    
    If the color {clr[]} is not in the table, and {add} is false, returns {NULL}.
    
    In any case, each coordinate of {clr[]} must be in {0..cht.maxval}. */

uint16_color_table_t *uint16_color_table_build
  ( uint16_t **samples,
    uint32_t chns, 
    uint32_t cols, 
    uint32_t rows, 
    uint16_t maxval,
    uint32_t NN_max
  );
  /* Returns a table {cht} with all the colors in the arrays
    {samples[0..rows-1][0..cols*chns-1]}.  The number {chsn}
    must be positive.
    
    Assumes that each list {samples[0..rows-1]} contains {cols} color tuples
    each consisting of {chns} samples, all in consecutive positions. Each color tuple
    is searched in the table, with {uint16_color_table_lookup}, and added if not there.
    
    If the {samples} arrays contain at most {NN_max} distinct color
    tuples, the returned table {cht} will be not {NULL}. The field
    {cht.NN} will be the number of such triplets, and the fields
    {cht.chns} and {cht.maxval} will be the given {chns} and
    {maxval}. The field {nd.count} of each node {nd} will
    be the (positive) number of times that the color tuple
    {nd.clr[-..chns-1]} occurred in the {samples} arrays, and the field
    {nd.index} will be {-1}.
    
    The {clr} field in each node will be a pointer into the first sample
    of the color tuple in the {samples} arrays. Thus those arrays should
    NOT be reclaimed while the table {cht} or its nodes are in use.
    
    However, if the {samples} arrays contain more than {NN_max} distinct
    colors, the procedure gives up and returns {NULL}. */

uint16_color_table_node_t** uint16_color_table_nodes(uint16_color_table_t *cht);
  /* Returns a single array of {cht->NN} pointers to all the nodes of all buckets of {cht},
    in random order. */

uint16_t* uint16_color_table_samples_from_nodes(uint32_t chns, uint32_t NN, uint16_color_table_node_t** nds);
  /* Returns a single flat array of {NN*chns} samples with copies of all the color tuples from all the nodes {nds[0..NN-1]},
    in the same order.  The nodes and the original sample arrays can then be freed. */

void uint16_color_table_free(uint16_color_table_t *cht);
  /* Reclaims all the storage used by the table {cht}, INCLUDING the 
    buckets array and all the nodes. */

#endif
