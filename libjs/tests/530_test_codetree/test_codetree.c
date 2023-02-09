#define PROG_NAME "test_codetree"
#define PROG_DESC "tests the {codetree.h} procedures"
#define PROG_VERS "1.1"

/* Last edited on 2023-02-07 17:11:19 by stolfi */
/* Created on 2007-01-31 by J. Stolfi, UNICAMP */

#define PROG_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <ref.h>
#include <affirm.h>
#include <jsrandom.h>

#include <codetree_limits.h>
#include <codetree.h>

#define MAX_VALUE codetree_MAX_VALUE
#define MIN_VALUE codetree_MIN_VALUE
#define MAX_LEAVES codetree_MAX_LEAVES
#define MAX_INTERNALS codetree_MAX_INTERNALS
#define MAX_NODES codetree_MAX_NODES
#define MAX_SAMPLES codetree_MAX_SAMPLES

void tcot_check_types(void);
  /* Checks if data types have sufficient size. */

void tcot_check_empty(void);
  /* Checks {codetree.h} functions with an empty ({NULL}) tree. */

void tcot_check_single(void);
  /* Checks {codetree.h} functions with a trivial tree (single leaf). */

void tcot_check_small(void);
  /* Checks {codetree.h} functions with a small hand-built tree. */

void tcot_check_tree(codetree_t *tree, uint32_t nv, uint64_t freq[]);
  /* Verifies if the tree is well-formed, including the child order 
    property with tie-breaking by smallest leaf value. */

void tcot_check_generic
  ( codetree_t *tree, 
    codetree_value_t maxval, 
    codetree_node_count_t nv, 
    codetree_value_t vals[], 
    char *code[]
  );
  /* Test {codetree.h} functions for the given {tree},
    whose value set {tree.V} is expected to be {vals[0..nv-1]} (that 
    should be a subset of {0..maxval}), and is expected to 
    define the bit code {code[val]} for each {val} in {V}. */

int32_t main (int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { tcot_check_types();
    tcot_check_empty();
    tcot_check_single();
    tcot_check_small();
    return 0;
  }
  
#define check_type(TYPE,VAR,ITYPE,VAL) \
  TYPE VAR = (TYPE)(VAL); \
  demand(((ITYPE)VAR) == ((ITYPE)(VAL)), "value " #VAL " doesn't fit in type " #TYPE)

void tcot_check_types(void)
  { 
    fprintf(stderr, "MAX_DELTA = %lu UINT32_MAX = %lu\n", (uint64_t)codetree_MAX_DELTA, (uint64_t)UINT32_MAX);
    
    check_type(codetree_value_t,        x1, int64_t,  codetree_MAX_VALUE);     
    check_type(codetree_value_t,        x2, int64_t,  codetree_MIN_VALUE);     
    check_type(codetree_node_count_t,   x3, uint64_t, codetree_MAX_NODES);     
    check_type(codetree_sample_count_t, x4, uint64_t, codetree_MAX_SAMPLES);   
    check_type(codetree_bit_count_t,    x5, uint64_t, codetree_MAX_BYTES * 8); 
    check_type(codetree_byte_count_t,   x6, uint64_t, codetree_MAX_BYTES);     
    check_type(byte_t,                  x7, uint64_t, 255);                    
    check_type(codetree_delta_t,        x8, uint64_t, codetree_MAX_DELTA);     
  }
   
void tcot_check_empty(void)
  { 
    codetree_value_t maxval = 40; /* Some upper bound on valid values ({V}). */
    char *code[maxval+1];         /* Expected bit code of each value. */
    for (codetree_value_t val = 0; val <= maxval; val++) { code[val] = NULL; }
    
    /* Expected tree and encoding: */
    codetree_node_count_t nv = 0; /* Number of valid values. */
    codetree_value_t vals[nv];

    codetree_t *tree = NULL;

    tcot_check_generic(tree, maxval, nv, vals, code);
  }
   
void tcot_check_single(void)
  { 
    codetree_value_t maxval = 40; /* Some upper bound on valid values ({V}). */
    char *code[maxval+1];
    for (codetree_value_t val = 0; val <= maxval; val++) { code[val] = NULL; }

    /* Expected tree and encoding: */
    codetree_node_count_t nv = 1; /* Number of valid values. */
    codetree_value_t vals[nv];

    vals[0] = 27; code[vals[0]] = "";
    
    codetree_node_t *L = codetree_new_leaf(vals[0]);
    codetree_t *tree = L;

    tcot_check_generic(tree, maxval, nv, vals, code);
 }
    
void tcot_check_small(void)
  { 
    codetree_value_t maxval = 40; /* Some upper bound on valid values ({V}). */
    char *code[maxval+1];
    for (codetree_value_t val = 0; val <= maxval; val++) { code[val] = NULL; }

    /* Expected tree and encodings: */
    codetree_node_count_t nv = 5; /* Number of valid values. */
    codetree_value_t vals[nv];

    vals[0] =  7; code[vals[0]] = "000";
    vals[1] = 14; code[vals[1]] = "001";
    vals[2] = 21; code[vals[2]] = "10";
    vals[3] = 28; code[vals[3]] = "11";
    vals[4] = 35; code[vals[4]] = "01";

    codetree_node_count_t seq = 0;
    codetree_node_t *L1 = codetree_new_leaf(vals[0]);
    codetree_node_t *L2 = codetree_new_leaf(vals[1]);
    codetree_node_t *L3 = codetree_new_leaf(vals[2]);
    codetree_node_t *L4 = codetree_new_leaf(vals[3]);
    codetree_node_t *L5 = codetree_new_leaf(vals[4]);
    
    codetree_node_t *M12 = codetree_new_internal(seq, L1, L2); seq++;
    codetree_node_t *M125 = codetree_new_internal(seq, M12, L5); seq++;
    codetree_node_t *M34 = codetree_new_internal(seq, L3, L4); seq++;
    codetree_node_t *M12534 = codetree_new_internal(seq, M125, M34); seq++;
    
    codetree_t *tree = M12534;
  
    tcot_check_generic(tree, maxval, nv, vals, code);
  }
 
void tcot_check_generic
  ( codetree_t *tree, 
    codetree_value_t maxval, 
    codetree_node_count_t nv, 
    codetree_value_t vals[], 
    char *code[]
  )
  { 
    bool_t prtables = (nv < 100); /* True to print the tables. */
    fprintf(stderr, "######################################################################\n");
    fprintf(stderr, "### begin {tcot_check_generic} - {maxval} = %u {nv} = %u\n", maxval, nv);
    demand((maxval >= 0) && (maxval <= MAX_VALUE), "bad {maxval} given");
    
    if (prtables)
      { fprintf(stderr, "given values and codes:\n");
        for (int32_t iv = 0; iv < nv; iv++)
          { fprintf(stderr, "%12d (%s)\n", vals[iv], code[vals[iv]]); }
        fprintf(stderr, "testing {codetree_print_codes}:\n");
        codetree_print_codes(stderr, tree);
      }
    
    fprintf(stderr, "testing {codetree_num_leaves}...\n");
    demand(nv == codetree_num_leaves(tree), "num leaves does not check");
    
    fprintf(stderr, "checking the codes implied by the tree...\n");
    codetree_node_count_t nv2 = codetree_check_codes(tree, maxval, code);
    demand(nv2 == nv, "inconsistent {nv2}");

    /* Get max code length: */
    int32_t maxlen = 0;
    for (codetree_value_t val = 0; val <= maxval; val++)
      { if (code[val] != NULL)
          { int32_t len = (int32_t)strlen(code[val]);
            if (len > maxlen) { maxlen = len; }
          }
      }
    
    if (tree != NULL)
      {   
        fprintf(stderr, "generating a sample sequence to encode...\n");
        codetree_sample_count_t ns = 20; /* Number of samples to encode. */
        bool_t prsamples = (ns < 100);   /* True to print the sample sequence. */
        codetree_value_t smp[ns];        /* Sample values to encode. */
        char *smpcode[ns];               /* Their codes, as given. */
        for(int32_t ks = 0; ks < ns; ks++) 
          { int32_t iv = int32_abrandom(0, nv-1);
            smp[ks] = vals[iv];
            smpcode[ks] = code[vals[iv]];
            if (prsamples) { fprintf(stderr, " %d", smp[ks]); }
          }
        if (prsamples) { fprintf(stderr, "\n"); }
        
        fprintf(stderr, "hand-encoding the samples with given codes...\n");
        codetree_byte_count_t nb = (ns*maxlen + 7)*8; /* Num of bytes in buf. Should be enough. */
        byte_t buf[nb];
        codetree_byte_count_t ib = 0; /* Next unused byte is {buf[ib]}. */
        byte_t mask = 0; /* Next bit in {buf[ib-1]} or 0. */
        int64_t nu = 0; /* Number of bits used. */
        for(int32_t ks = 0; ks < ns; ks++) 
          { char *ck = smpcode[ks];
            while ((*ck) != 0)
              { int8_t bit = ((*ck) - '0');
                if (mask == 0)
                  { assert(ib < nb);
                    ib++;
                    buf[ib-1] = 0;
                    mask = 128;
                  }
                if (bit == 1) { buf[ib-1] |= mask; }
                mask = mask/2;
                nu++;
                ck++;
              }
          }
        int64_t mb = (nu + 7)/8; /* Actual number of bytes used. */
        fprintf(stderr, "used %lu bits in %lu bytes ({ib} = %lu)\n", nu, mb, ib);
        assert(mb == ib);
        demand(mb < nb, "impossible used bit count {nu}");
        
        fprintf(stderr, "hand-encoded bits = "); 
        codetree_print_bits(stderr, (mb < 3 ? mb : 3), buf, " ");
        if (mb > 3) { fprintf(stderr, " ..."); }
        fprintf(stderr, "\n"); 
        
        fprintf(stderr, "testing {codetree_decode} on the encoded string...\n");
        codetree_value_t smp2[ns];
        int64_t nu2 = codetree_decode(nb, buf, tree, maxval, ns, smp2);
        demand(nu2 == nu, "decoded bit counts do not match");
        for (int32_t ks = 0; ks < ns; ks++) 
          { demand(smp2[ks] == smp[ks], "decoded samples do not match"); }
        
        fprintf(stderr, "testing {codetree_get_encoding_table}...\n");
        int32_t nd = maxval + nv; /* Size of encoding table, {(maxval+1)+(|V|-1)}. */
        codetree_delta_t delta[nd];
        int64_t md = codetree_get_encoding_table(tree, maxval, nd, delta);
        demand(md == nd, "unexpected encoding table size");
        if (prtables)
          { for (int32_t id = 0; id < nd; id++)
              { fprintf(stderr, "%3d", id);
                codetree_delta_t del = delta[id];
                if (del == 0)
                  { fprintf(stderr, "   -\n"); }
                else
                  { int8_t ich = del % 2;
                    codetree_node_count_t ip = id + del/2;
                    fprintf(stderr, " %3u:%d\n", ip, ich);
                  }
               }
           }
        
        fprintf(stderr, "encoding the samples with {codetree_encode} and the encoding table...\n");
        byte_t buf3[nb];
        for (int32_t ib = 0; ib < nb; ib++) { buf3[ib] = 17; /* Arbitrary nonzero value: */ }
        int64_t nu3 = codetree_encode(ns, smp, maxval, nd, delta, nb, buf3);
        demand(nu3 == nu, "num bits encoded with {delta} does not match"); 
        
        fprintf(stderr, "comparing the two bit strings...\n");
        for (int32_t ib = 0; ib < nb; ib++)
          { if (ib < mb) 
              { demand(buf3[ib] == buf[ib], "encoded bit strings differ"); }
            else
              { demand(buf3[ib] == 17, "wrote too much into {buf3}"); }
          }
      }
    fprintf(stderr, "calling {codetree_check_tree}...\n");
    codetree_check_tree(tree, maxval);
    fprintf(stderr, "### end {tcot_check_generic}\n");
    fprintf(stderr, "######################################################################\n");
  }

