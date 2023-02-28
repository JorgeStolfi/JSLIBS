/* Huffman tree and Huffman decoding. */
/* Last edited on 2023-02-19 00:37:36 by stolfi */

#ifndef codetree_huff_H
#define codetree_huff_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <codetree.h>

typedef uint64_t codetree_huff_freq_t;
  /* Type of a frequency count for this module. */
  
#define codetree_huff_MAX_FREQ UINT64_MAX
  /* Maximum value of of a frequency count. */

codetree_t *codetree_huff_build(codetree_value_t maxval, codetree_huff_freq_t freq[]);
  /* Builds the Huffman tree {T} for a set of values {V},
    based on the frequency counts {freq[0..maxval]}. Returns a pointer 
    to the root node.  
    
    The tree can be used with {codetree_decode} to decode a bit string
    into a sequence of values of {V}.  For the encoding, 
    one can use {codetree_encode} with the table returned by 
    {codetree_get_encoding_table}.
    
    The set {V} consists of every integer {val} in {0..maxval} such that
    {freq[iv]} is not zero. If {V} is empty, the procedure returns
    {NULL}.  
    
    The number {maxval} must not exceed {codetree_MAX_VALUE}.
    The sum of all entries of {freq} must not exceed {codetree_huff_MAX_FREQ}.
    
    HUFMMAN TREE INVARIANTS

    The returned tree will have the following properties, besides those
    described under {codetree_t}. With respect to weights:

     * Each node of the tree has an implicit weight. The weight of a leaf
       node {p} is the frequency associated to its value, that is,
       {freq[p.value]}. The weight of an internal node is defined as the
       sum of the weights of all leaves that descend from it.

     * In any internal node {p}, the weight of the left child {p.child[0]}
       is greater than or equal to that of the right child {p.child[1]}.

    The tree also has the following properties that invove the order of
    creation of nodes:

      * The leaves are created in order of increasing value,
        irrespective of their weight.

      * Every internal node is created after all leaf nodes.

      * Every internal node is created after its two children.

      * The two children of an internal node {p} are the two nodes 
        with smaller weight created before {p}.  
        
      * If two of the nodes {q1} and {q2} created before {p}
        have the same weight, the one that was created
        last is sorted as if it had smaller weight.

    Thus, in particular, if the two children of a node {p}
    have the same weight,
    
      * If they are both leaves, the left child has the smallest
        {value}.
        
      * If only one of them is a leaf, it is the left child.
      
      * If both are internal nodes, the left child was created
        before the right one.
        
    The order of creation of internal nodes will be recorded in their
    {value} fields. They are strictly negative and all distinct. The
    first internal node created will have its {value} set to
    {codetree_MIN_VALUE}, the next one to {codetree_MIN_VALUE+1}, and so on. It follows
    that a node {p1} (internal or leaf) was created before node {p2}
    (ditto) if and only if {((uint64_t)p1.value) <
    ((uint64_t)p2.value)}.
   */

void codetree_huff_check_huffman(codetree_t *tree);
  /* Checks the Huffman properties of the given {tree}. */ 
  
#endif
