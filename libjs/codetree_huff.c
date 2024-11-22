/* See {codtree_huff.h}. */
/* Last edited on 2024-11-20 03:28:10 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <bool.h>
#include <affirm.h>

#include <codetree.h>
#include <codetree_huff.h>

#define MAX_VALUE codetree_data_MAX_VALUE
#define MAX_LEAVES codetree_MAX_LEAVES
#define MAX_FREQ codetree_huff_MAX_FREQ

void codetree_huff_print_queue(char *tag, codetree_node_count_t nv, codetree_node_t *node[], codetree_huff_freq_t weight[]);
  /* Prints a few entries of the Huffman algorithm queue 
    {node[0..nv-1],weight[0..nv-1]} for debugging, on s aingle line, 
    preceded by "{tag} ". Eachentry is printed with {codetree_huff_print_queue_entry}. */

void codetree_huff_print_queue_entry(codetree_node_t *nd, codetree_huff_freq_t wt);
  /* Prints the Huffman queue entry {(nd,wt)} for debugging, as "[{nd}:{wt}]".
    The node {nd} is  printed with {codetree_huff_print_node}. */

void codetree_huff_print_node(codetree_node_t *nd);
  /* Prints the Huffman tree node {nd} for debugging.
    If {nd} is a leaf, prints its {value} in "%+d" format; otherwise 
    prints it as "@{seq}" where {seq} is the node's sequence number. */

void codetree_huff_sort_by_weight(codetree_node_count_t nv, codetree_node_t *node[], codetree_huff_freq_t weight[]);
  /* Sorts the lists {weight[0..nv-1]} and {node[0..nv-1]} by decreasing
    {weight[i]}. In case of ties, preserves the original order of the
    nodes (which should be the order of increasing value for leaves, and the 
    order creation for internal nodes).
    
    The running time is quadratic on first call, but should be linear
    when only the last element is out of order. */

codetree_t *codetree_huff_build(codetree_data_value_t maxval, codetree_huff_freq_t freq[])
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
    for (codetree_data_value_t val = 0; val <= maxval; val++) { if (freq[val] != 0) { nvalid++; } }
    
    /* Working vectors: */
    codetree_huff_freq_t *weight = (uint64_t *)notnull(malloc(nvalid*sizeof(codetree_huff_freq_t)), "no mem");
    codetree_node_t **node = (codetree_node_t **)notnull(malloc(nvalid*sizeof(codetree_node_t *)), "no mem");

    /* Build the leaf nodes and copy their frequencies: */
    codetree_node_count_t iv = 0;
    for (codetree_data_value_t val = 0; val <= maxval; val++) 
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
            codetree_huff_sort_by_weight(nv, node, weight);
            if (debug) { codetree_huff_print_queue("\nS", nv, node, weight); }
            codetree_node_t *child0 = node[nv-2];
            codetree_node_t *child1 = node[nv-1];
            codetree_node_t *p = codetree_new_internal(nint, child0, child1);
            nint++;
            node[nv-2] = p;
            demand(weight[nv-2] <= UINT64_MAX - weight[nv-1], "total {freq} values too large, overflow");
            weight[nv-2] += weight[nv-1];
            nv--;
            if (debug) { codetree_huff_print_queue("M", nv, node, weight); }
          }
        assert(nint == nvalid-1);
        root = node[0];
      }
      
    free(node);
    free(weight);
    return root;
  }
  
void codetree_huff_sort_by_weight(codetree_node_count_t nv, codetree_node_t *node[], codetree_huff_freq_t weight[])
  { bool_t debug = FALSE;

    /* Insertion sort: */
    for (int64_t iv = 1; iv < nv; iv++)
      { uint64_t fri = weight[iv];
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
            
void codetree_huff_print_queue(char *tag, codetree_node_count_t nv, codetree_node_t *node[], codetree_huff_freq_t weight[])
  { fprintf(stderr, "%s ", tag);
    for (codetree_node_count_t iv = 0; iv < nv; iv++)
      { if ((iv < 3) || (iv == nv-1))
          { if (iv > 0) { fprintf(stderr, " "); }
            codetree_huff_print_queue_entry(node[nv-1-iv], weight[nv-1-iv]);
          }
        else if ((nv > 4) && (iv == 3))
          { fprintf(stderr, " ... "); }
      }
    codetree_node_t *nu = node[nv-1]; /* Last node. */
    if (nu->value < 0)
      { /* Print children of first node: */
        fprintf(stderr, " ");
        codetree_huff_print_node(nu);
        fprintf(stderr, "=(");
        codetree_huff_print_node(nu->child[0]);
        fprintf(stderr, ",");
        codetree_huff_print_node(nu->child[1]);
        fprintf(stderr, ")");
      }
    fprintf(stderr, "\n");
  }

void codetree_huff_print_queue_entry(codetree_node_t *nd, codetree_huff_freq_t wt)
  { fprintf(stderr, "[");
    codetree_huff_print_node(nd);
    fprintf(stderr, ":%lu", wt);
    fprintf(stderr, "]");
  }
  
void codetree_huff_print_node(codetree_node_t *nd)
  { demand(nd != NULL, "null node in Huffman queue");
    if (nd->value >= 0)
      { fprintf(stderr, "%+d", nd->value); }
    else
      { codetree_node_count_t seq = (codetree_node_count_t)(nd->value - codetree_node_MIN_VALUE);
        fprintf(stderr, "@%04u", seq);
      }
  }
  
