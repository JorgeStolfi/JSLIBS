/* See {huff_tree.h}. */
/* Last edited on 2023-02-06 17:46:35 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <jsstring.h>
#include <fget.h>
#include <affirm.h>

#include <huff_tree.h>

#define MAX_VALUE codetree_MAX_VALUE
#define MIN_VALUE codetree_MIN_VALUE
#define MAX_LEAVES codetree_MAX_LEAVES
#define MAX_FREQ huff_tree_MAX_FREQ

void huff_tree_print_entry(char *tag, codetree_node_count_t iv, huff_tree_freq_t weight, codetree_node_t *node);
  /* Prints the queue entry {iv,weight[iv],node[iv]} for debugging. */

void huff_tree_sort_by_weight(codetree_node_count_t nv, codetree_node_t *node[], huff_tree_freq_t weight[]);
  /* Sorts the lists {weight[0..nv-1]} and {node[0..nv-1]} by decreasing
    {weight[i]}. In case of ties, preserves the original order of the
    nodes (which should be the order of increasing value for leaves, and the 
    order creation for internal nodes).
    
    The running time is quadratic on first call, but should be linear
    when only the last element is out of order. */

codetree_node_t *huff_tree_build(codetree_value_t maxval, huff_tree_freq_t freq[])
  {
    bool_t debug = FALSE;
    
    /* The procedure first creates one leaf node {p[val]} for each value
      {val} in {0..maxval} (that is, with {p[val].value = val})
      that has nonzero frequency {freq[val]}, which is its weight.
      Then it repeatedly finds the two nodes with
      smallest weight, and replaces them by a new node which is their
      parent, whose weight is the sum of their weights; until a single
      node remains, which is the root of the tree.
      
      The {value} field of new nodes is initialized with negative
      numbers starting with {MIN_VALUE} and incrementing. Those numbers
      are all negative. Given the way that leaves are created, comparing
      {value} fields as unsigned integers gives the creation order of
      any two nodes, whether leaf or not.

      The children are sorted so that the weight of the left child is
      greater than or equal to that of the right child.  In case of ties,
      the left child is the subtree that was created earlier.
    */
    
    /* Count valid values: */
    codetree_node_count_t nvalid = 0; /* Number of valid values in {V}: */
    for (codetree_value_t val = 0; val <= maxval; val++) { if (freq[val] != 0) { nvalid++; } }
    
    /* Working vectors: */
    huff_tree_freq_t *weight = (uint64_t *)notnull(malloc(nvalid*sizeof(huff_tree_freq_t)), "no mem");
    codetree_node_t **node = (codetree_node_t **)notnull(malloc(nvalid*sizeof(codetree_node_t *)), "no mem");

    /* Build the leaf nodes and copy their frequencies: */
    codetree_node_count_t iv = 0;
    for (codetree_value_t val = 0; val <= maxval; val++) 
      { if (freq[val] > 0)
          { weight[iv] = freq[val];
            node[iv] = codetree_new_leaf(val);
            iv++;
          }
      }
    assert(iv == nvalid);

    /* Huffman tree algorithm: */
    codetree_node_t *root = NULL;
    if (nvalid > 0)
      { /* Create the tree: */
        codetree_node_count_t nint = 0;    /* Number of internal nodes created so far. */
        codetree_node_count_t nv = nvalid; /* The nodes yet to be combined are {node[0..nv-1]} with weight {weight[0..nv-1]}. */ 
        while (nv >= 2)
          { /* (Re)sort {weight[0..nv-1],node[0..nv-1]} by decreasing {weight}: */
            huff_tree_sort_by_weight(nv, node, weight);
            if (debug)
              { if (nv > 2) 
                  { huff_tree_print_entry(" ", 0, weight[0], node[0]);
                    fprintf(stderr, "    ...\n");
                  }
                huff_tree_print_entry(" ", nv-2, weight[nv-2], node[nv-2]);
                huff_tree_print_entry("+", nv-1, weight[nv-1], node[nv-1]);
              }
            codetree_node_t *child0 = node[nv-2];
            codetree_node_t *child1 = node[nv-1];
            codetree_node_t *p = codetree_new_internal(nint, child0, child1);
            nint++;
            node[nv-2] = p;
            demand(weight[nv-2] <= UINT64_MAX - weight[nv-1], "total {freq} values too large, overflow");
            weight[nv-2] += weight[nv-1];
            if (debug)
              { huff_tree_print_entry("=", nv-2, weight[nv-2], node[nv-2]);
                fprintf(stderr, "\n");
              }
            nv--;
          }
        assert(nint == nvalid-1);
        root = node[0];
      }
      
    free(node);
    free(weight);
    return root;
  }
  
void huff_tree_sort_by_weight(codetree_node_count_t nv, codetree_node_t *node[], huff_tree_freq_t weight[])
  { bool_t debug = FALSE;

    /* Insertion sort: */
    for (int64_t iv = 1; iv < nv; iv++)
      { int64_t fri = weight[iv];
        codetree_node_t *pi = node[iv];
        int64_t jv = iv;
        while ((jv > 0) && (weight[jv-1] < fri))
          { weight[jv] = weight[jv-1];
            node[jv] = node[jv-1];
            jv--;
          }
        if (jv != iv) 
          { weight[jv] = fri;
            node[jv] = pi; 
          }
      }

    if (debug)
      { /* Check sorting: */
        for (int64_t iv = 1; iv < nv; iv++) 
          { int64_t jv = iv-1;
            assert(weight[jv] >= weight[iv]);
            if (weight[jv] == weight[iv])
              { /* Compare order of creation: */
                assert(((uint64_t)node[jv]->value) < ((uint64_t)node[iv]->value));
              }
          }
      }
  }

void huff_tree_print_entry(char *tag, codetree_node_count_t iv, huff_tree_freq_t weight, codetree_node_t *node)
  { fprintf(stderr, "  %s weight[%3u] = %lu", tag, iv, weight);
    demand(node != NULL, "null node in Huffman queue");
    if (node->value >= 0)
      { fprintf(stderr, " val = %d\n", node->value); }
    else
      { fprintf(stderr, " seq = %d\n", node->value - MIN_VALUE); }
  }
  
