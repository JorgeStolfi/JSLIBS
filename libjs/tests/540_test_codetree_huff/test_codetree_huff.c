#define PROG_NAME "test_codetree_huff"
#define PROG_DESC "tests the {codetree_huff.h} procedures"
#define PROG_VERS "1.1"

/* Last edited on 2023-02-15 12:24:59 by stolfi */
/* Created on 2007-01-31 by J. Stolfi, UNICAMP */

#define PROG_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <ref.h>
#include <affirm.h>

#include <codetree.h>
#include <codetree_huff.h>
     
#define MAX_VALUE codetree_MAX_VALUE
#define MIN_VALUE codetree_MIN_VALUE
#define MAX_LEAVES codetree_MAX_LEAVES
#define MAX_INTERNALS codetree_MAX_INTERNALS
#define MAX_NODES codetree_MAX_NODES
#define MAX_SAMPLES codetree_MAX_SAMPLES

#define MAX_FREQ codetree_huff_MAX_FREQ

void thuf_check_types(void);
  /* Checks if data types have sufficient size. */

void thuf_check_empty(void);
  /* Test {codetree_huff_build} with an empty set {V}. */
    
void thuf_check_single(void);
  /* Test {codetree_huff_build} with a singleton set {V}. */
   
void thuf_check_small(void);
  /* Test {codetree_huff_build} with a simple example that requires tie-breaking. */
   
void thuf_check_unif(void);
  /* Test {codetree_huff_build} with a uniform freq distribution. */
    
void thuf_check_ruler(void);
  /* Test {codetree_huff_build} with a very non-uniform freq distribution. */

void thuf_check_generic
  ( codetree_value_t maxval, 
    codetree_huff_freq_t freq[],
    codetree_t *tree, 
    char *code[]
  );
  /* Test {codetree_huff_build} for a set {V} of values given the
    frequencies {freqs[0..maxval]}. The set {V} consists of every value
    {val} in {0..maxval} with {freq[val] > 0}. 
    
    If {code} is not {NULL}, also checks whether the
    resulting tree is isomorphic to the given {tree},
    and whether it defines the bit codes {code[val]} for each {val} in {V}.
    If {code} is {NULL}, skips this part. */

codetree_node_count_t thuf_check_tree
  ( codetree_t *tree, 
    codetree_value_t maxval, 
    codetree_huff_freq_t freq[]
  );
  /* Verifies if the given {tree} is well-formed, including the child order 
    property with tie-breaking by smallest leaf value.
    
    Returns the total number of nodes in the tree. */

codetree_huff_freq_t thuf_ruler_func(codetree_node_count_t iv);
  /* Returns a "ruler function" of the non-negative integer {iv}. */
  
int32_t main (int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { thuf_check_types();
    thuf_check_empty();
    thuf_check_single();
    thuf_check_small();
    thuf_check_unif();
    thuf_check_ruler();
    return 0;
  }
  
#define check_type(TYPE,VAR,ITYPE,VAL) \
  TYPE VAR = (TYPE)(VAL); \
  demand(((ITYPE)VAR) == ((ITYPE)(VAL)), "value " #VAL " doesn't fit in type " #TYPE)

void thuf_check_types(void)
  { 
    fprintf(stderr, "checking types...\n");
    check_type(codetree_value_t,        x1, int64_t,  codetree_MAX_VALUE);     
    check_type(codetree_value_t,        x2, int64_t,  codetree_MIN_VALUE);     
    check_type(codetree_node_count_t,   x3, uint64_t, codetree_MAX_NODES);     
    check_type(codetree_sample_count_t, x4, uint64_t, codetree_MAX_SAMPLES);   
    check_type(codetree_bit_count_t,    x5, uint64_t, codetree_MAX_BYTES * 8); 
    check_type(codetree_byte_count_t,   x6, uint64_t, codetree_MAX_BYTES);     
    check_type(byte_t,                  x7, uint64_t, 255);                    
    check_type(codetree_delta_t,        x8, uint64_t, codetree_MAX_DELTA);     

    check_type(codetree_huff_freq_t,        x9, uint64_t, codetree_huff_MAX_FREQ);     
  }
   
void thuf_check_empty(void)
  { 
    fprintf(stderr, "testing with empty {V}...\n");
    codetree_value_t maxval = 40; /* Some upper bound on valid values ({V}). */
    codetree_huff_freq_t freq[maxval+1];
    char *code[maxval+1];
    for (codetree_value_t val = 0; val <= maxval; val++) { freq[val] = 0; code[val] = NULL; }
    
    /* Expected tree: */
    codetree_t *tree = NULL;

    thuf_check_generic(maxval, freq, tree, code);
  }
   
void thuf_check_single(void)
  { 
    fprintf(stderr, "testing with singleton {V}...\n");
    codetree_value_t maxval = 40; /* Some upper bound on valid values ({V}). */
    codetree_huff_freq_t freq[maxval+1];
    char *code[maxval+1];
    for (codetree_value_t val = 0; val <= maxval; val++) { freq[val] = 0; code[val] = NULL; }

    /* Define the elements of {V}, their freqs, and their expected encodings: */
    freq[27] = 100; code[27] = "";

    /* Build the expected {tree}: */
    codetree_node_t *L = codetree_new_leaf(27);
    codetree_t *tree = L;

    thuf_check_generic(maxval, freq, tree, code);
 }
    
void thuf_check_small(void)
  { 
    fprintf(stderr, "testing with small fixed tree...\n");
    codetree_value_t maxval = 40; /* Some upper bound on valid values ({V}). */
    codetree_huff_freq_t freq[maxval+1];
    char *code[maxval+1];
    for (codetree_value_t val = 0; val <= maxval; val++) { freq[val] = 0; code[val] = NULL; }
    
    /* Define the elements of {V}, their freqs, and their expected encodings: */
    codetree_node_count_t nv = 5; /* Number of leaves. */
    codetree_value_t vals[nv];
    vals[0] =  7; freq[vals[0]] = 120; code[vals[0]] = "01";
    vals[1] = 14; freq[vals[1]] =  70; code[vals[1]] = "100";
    vals[2] = 21; freq[vals[2]] = 130; code[vals[2]] = "00";
    vals[3] = 28; freq[vals[3]] = 110; code[vals[3]] = "11";
    vals[4] = 35; freq[vals[4]] =  50; code[vals[4]] = "101";

    /* A simple tree with 5 values: */
    codetree_node_count_t seq = 0;
    codetree_node_t *L0 = codetree_new_leaf(vals[0]);
    codetree_node_t *L1 = codetree_new_leaf(vals[1]);
    codetree_node_t *L2 = codetree_new_leaf(vals[2]);
    codetree_node_t *L3 = codetree_new_leaf(vals[3]);
    codetree_node_t *L4 = codetree_new_leaf(vals[4]);
    
    codetree_node_t *M20 =    codetree_new_internal(seq, L2, L0); seq++;
    codetree_node_t *M14 =    codetree_new_internal(seq, L1, L4); seq++;
    codetree_node_t *M143 =   codetree_new_internal(seq, M14, L3); seq++;
    codetree_node_t *M20143 = codetree_new_internal(seq, M20, M143); seq++;
    
    codetree_t *tree = M20143;
    
    thuf_check_generic(maxval, freq, tree, code);
  }

void thuf_check_unif(void)
  { 
    fprintf(stderr, "testing with uniform frequencies...\n");
    codetree_value_t maxval = 100; /* Some upper bound on valid values ({V}). */
    codetree_huff_freq_t freq[maxval+1];
    char **code = NULL;
    for (codetree_value_t val = 0; val <= maxval; val++) { freq[val] = 0; }
    
    /* Define the elements of {V} and their freqs: */
    codetree_node_count_t nv = 64; /* Number of values in {V}. Better be power of 2. */
    for (codetree_node_count_t iv = 0; iv < nv; iv++)
      { codetree_value_t val = maxval*iv/nv;
        freq[val] = 100;
      }

    codetree_t *tree = NULL; /* No expected tree. */

    thuf_check_generic(maxval, freq, tree, code);
  }
   
void thuf_check_ruler(void)
  { fprintf(stderr, "testing with 'ruler' frequencies...\n");
    codetree_value_t maxval = 100; /* Some upper bound on valid values ({V}). */
    codetree_huff_freq_t freq[maxval+1];
    char **code = NULL;
    for (codetree_value_t val = 0; val <= maxval; val++) { freq[val] = 0; }

    /* Define the elements of {V} and their freqs: */
    codetree_node_count_t nv = 63; /* Number of values in {V}. Better be power of 2 minus 1. */
    for (codetree_node_count_t iv = 0; iv < nv; iv++) 
      { codetree_value_t val = maxval*iv/nv;
        freq[val] = 100*thuf_ruler_func(iv);
      }
      
    codetree_t *tree = NULL; /* No expected tree. */

    thuf_check_generic(maxval, freq, tree, code);
  }
        
codetree_huff_freq_t thuf_ruler_func(codetree_node_count_t iv)
  { demand(iv >= 0, "bad {iv}");
    uint32_t kv = iv+1;
    uint64_t fr = 1; /* Frequency of leaf. */
    uint64_t sf = 1; /* Total freq of subtree. */
    while (kv % 2 == 0) 
      { fr = 2*sf + 1; 
        sf = 2*sf + fr; 
        kv = kv/2;
      }
    return fr;
  }
    
void thuf_check_generic
  ( codetree_value_t maxval, 
    codetree_huff_freq_t freq[],
    codetree_t *tree, 
    char *code[]
  )
  {
    fprintf(stderr, "building a tree {tree2} for the given {freq} vector...\n");
    codetree_t *tree2 = codetree_huff_build(maxval, freq);
    
    fprintf(stderr, "counting leaves of {tree2}...\n");
    codetree_node_count_t nv2 = codetree_num_leaves(tree2);
    fprintf(stderr, "found %u leaves\n", nv2);
    
    fprintf(stderr, "values and codes in {tree2}\n");
    codetree_print_codes(stderr, tree2);
    
    fprintf(stderr, "checking {tree2} for general codetree properties...\n");
    codetree_check_tree(tree2, maxval);
    
    fprintf(stderr, "checking {tree2}for Huffman-specific properties...\n");
    codetree_node_count_t nv3 = thuf_check_tree(tree2, maxval, freq);
    demand(nv3 == nv2, "inconsistent leaf counts {nv3,nv2}");

    if (code != NULL)
      { /* We have a reference tree {tree} and bit codes {code[0...maxval]}: */
        fprintf(stderr, "checking isomorphism of {tree2} and given tree {tree}...\n");
        codetree_check_iso(tree2, tree);
        
        fprintf(stderr, "checking if the bit encoding of {tree2} is the expected one...\n");
        codetree_node_count_t nv1 = codetree_check_codes(tree2, maxval, code);
        demand(nv1 == nv2, "inconsistent leaf counts {nv1,nv2}");
      }
    
    if (tree2 != NULL)
      { fprintf(stderr, "doing {codetree_free} on {tree2}...\n");
        codetree_node_count_t nf2 = codetree_free(tree2);
        assert(nf2 == 2*nv2 - 1);
      }

    return;
  }

codetree_node_count_t thuf_check_tree
  ( codetree_t *tree, 
    codetree_value_t maxval, 
    codetree_huff_freq_t freq[]
  )
  { 
    /* Get count of representable values: */
    codetree_node_count_t nvalid = 0;
    for (codetree_value_t val = 0; val <= maxval; val++) { if (freq[val] > 0) { nvalid++; } }
    
    /* Check Huffman tree invariants on child order: */

    auto codetree_huff_freq_t checkit(codetree_node_t *nd);
      /* Checks the child ordering of all internal nodes of the 
        tree rooted at {nd}. Returns its total total leaf weight.
        Also increments {nleaf}. */
        
    codetree_node_count_t nleaf = 0; /* Number of leaves in tree. */
    if (tree == NULL) 
      { fprintf(stderr, "tree is empty\n"); }
    else
      { codetree_huff_freq_t tree_weight = checkit(tree);
        fprintf(stderr, "total tree weight = %lu\n", tree_weight);
      }
    demand(nleaf == nvalid, "leaf count does not match {V} size");
      
    return nvalid;
  
    /* IMPS */
    
    codetree_huff_freq_t checkit(codetree_node_t *nd)
      { demand(nd != NULL, "child pointer is null");
        if (nd->value < 0)
          { /* Internal node: */
            uint64_t wt_ch[2]; /* Total leaf weight of each subtree. */
            uint64_t ct_ch[2]; /* Creation "time" for tie-breaking. */
            for (uint32_t ich = 0; ich < 2; ich++)
              { codetree_node_t *ch = nd->child[ich];
                wt_ch[ich] = checkit(ch);
                ct_ch[ich] = (uint64_t)ch->value;
              }
            if (wt_ch[0] < wt_ch[1]) 
              { fprintf(stderr, "** node has children in wrong order, weights = %lu %lu\n", wt_ch[0], wt_ch[1]); 
                demand(FALSE, "** aborted");
              }
            else if ((wt_ch[0] == wt_ch[1]) && (ct_ch[0] > ct_ch[1]))
              { fprintf(stderr, "** bad child tie-breaking = %lu %lu\n", ct_ch[0], ct_ch[1]);
                demand(FALSE, "** aborted");
              }
            if((uint64_t)(MAX_FREQ - wt_ch[0]) < (uint64_t)(wt_ch[1]))
              { fprintf(stderr, "** overflow of tree weight = %lu + %lu > %lu\n", wt_ch[0], wt_ch[1], MAX_FREQ);
                demand(FALSE, "** aborted");
              }
            return wt_ch[0] + wt_ch[1];
          }
        else
          { /* Leaf node: */
            nleaf++;
            codetree_value_t val = nd->value;
            if (val > maxval) 
              { fprintf(stderr, "** leaf value = %d out of range 0..%d\n", val, maxval);
                demand(FALSE, "** aborted");
              }
            if (freq[val] == 0)
              { fprintf(stderr, "** leaf value = %d has zero frequency\n", val);
                demand(FALSE, "** aborted");
              }
            return freq[val];
          }
      }
  }
