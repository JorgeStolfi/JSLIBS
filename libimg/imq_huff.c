/* See {imq_huff.h}. */
/* Last edited on 2024-12-05 07:54:11 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <jsstring.h>
#include <jsprintf.h>
#include <fget.h>
#include <affirm.h>

#include <codetree_huff.h>
#include <imq_huff.h>

#define MAX_VALUE codetree_MAX_VALUE
#define MIN_VALUE codetree_MIN_VALUE
#define MAX_LEAVES codetree_MAX_LEAVES
#define MAX_FREQ codetree_huff_MAX_FREQ

codetree_t *imq_huff_build_tree(codetree_huff_freq_t freq[])
  {
    uint32_t nh = 511;
    codetree_huff_freq_t hist[nh];
    
    /* Complement the diffs to get proper tie-breaking in Huffman tree: */
    for (int32_t ih = 0;  ih < nh; ih++) { hist[ih] = freq[(int32_t)nh-1-ih]; }
    
    /* Build the Huffman tree for the set {V} = {0..510}: */
    codetree_data_value_t maxval = nh-1;
    codetree_t *tree = codetree_huff_build(maxval, hist);
    return tree;
  }
  
codetree_bit_count_t imq_huff_decode
  ( uint32_t line_num,
    codetree_byte_count_t nb, 
    byte_t buf[], 
    codetree_t *tree, 
    codetree_byte_count_t np, 
    byte_t pix[]
  )
  {
    bool_t verbose = FALSE;
    bool_t debug = verbose && (line_num == 1);
    
    codetree_data_value_t maxval = 510; /* Max decoded value: {0..510} -> {-255..+255}. */
    
    codetree_sample_count_t nd = np-1; /* Number of pixel diffs. */
    codetree_byte_count_t nb1 = nb-1;  /* Num bytes of Huffman codes. */
    byte_t *buf1 = buf+1;              /* The Huffman code bytes. */
    codetree_data_value_t diff[nd];
    if (verbose) 
      { fprintf(stderr, "Huffman-decoding %lu bytes", nb1);
        fprintf(stderr, " into %lu complemented differences...\n", nd);
        fprintf(stderr, "first few Huffman code bits = ");
        codetree_print_bits(stderr, (nb1 < 3 ? nb1 : 3), buf1, " ");
        if (nb1 > 3) { fprintf(stderr, " ..."); }
        fprintf(stderr, "\n");
      }
      
    /* We must complement the Huffman bytes due way IMQ interprets them: */
    for (int32_t ib = 0;  ib < nb1; ib++) { buf1[ib] ^= (byte_t)255; }
    codetree_bit_count_t nbits_decoded = codetree_decode(nb1, buf1, tree, maxval, nd, diff);
    codetree_byte_count_t nbytes_decoded = (nbits_decoded + 7)/8;
    assert(nbytes_decoded <= nb);
    if (verbose)
      { fprintf(stderr, "parsed %lu bits (%lu bytes)", nbits_decoded, nbytes_decoded);
        fprintf(stderr, " out of %lu (%lu)\n", (nb-1)*8, nb-1);
      }
    
    if (verbose) { fprintf(stderr, "integrating the differences...\n"); }
    pix[0] = buf[0];
    for (int32_t i = 0;  i < nd; i++)
      { /* Note that we reversed signs when building the tree. */
        /* So each {diff[i]} is {255 - (pix[i+1]-pix[i])}. */
        int32_t di = (int32_t)diff[i] - 255; 
        int32_t val = pix[i] + di;
        if (debug) 
          { fprintf(stderr, " %4u %+4d = %4d\n", pix[i], di, val); }
        demand((val >= 0) && (val <= 255), "integrated pixel outside {0..255}");
        pix[i+1] = (byte_t)(val & 255);
      }
    return nbits_decoded;
  }
  
void imq_huff_print_codes(FILE *stderr, codetree_t *tree)
  { 
    codetree_data_value_t maxval = 510;
    demand(tree != NULL, "empty tree");

    int32_t nv = 0;
    
    auto void enum_leaves(char *pref, codetree_node_t *p);
      /* Prints the values and codes of the subtree rooted at {p},
        assuming that its code is the string {pref}. Also increments {nv}
        for each leaf found. Note that in the IMQ encoding
        '0' means go right and '1' means go left. */
        
    enum_leaves("", tree);
    
    demand(nv == 511, "wrong number of leaves");
    
    return;
    
    void enum_leaves(char *pref, codetree_node_t *p)
      { if (p->value < 0)
          { for (int8_t ich = 0; ich < 2; ich++)
              { char *pref_ch = jsprintf("%s%d", pref, 1-ich);
                enum_leaves(pref_ch, p->child[ich]);
                free(pref_ch);
              }
          }
        else
          { demand(p->value <= maxval, "invalid leaf value");
            /* Note that we complemented the diffs before building the tree: */
            int32_t diff = p->value - 255;
            fprintf(stderr, "%12d = %+4d (%s)\n", p->value, diff, pref);
            nv++;
          }
      }
  }

