/* See {codetree.h}. */
/* Last edited on 2024-11-20 06:51:50 by stolfi */

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

#include <codetree.h>

#define MAX_VALUE codetree_data_MAX_VALUE
#define MAX_LEAVES codetree_MAX_LEAVES
#define MAX_INTERNALS codetree_MAX_INTERNALS
#define MAX_NODES codetree_MAX_NODES
#define MAX_SAMPLES codetree_MAX_SAMPLES
#define MAX_BYTES codetree_MAX_BYTES
#define MAX_BITS codetree_MAX_BITS

/* IMPS */

codetree_node_count_t codetree_free(codetree_t *tree)
  { if (tree == NULL) 
      { return 0; }
    else 
      { codetree_node_count_t nf = 1;
        if (tree->value < 0)
          { nf += codetree_free(tree->child[0]);
            nf += codetree_free(tree->child[1]);
          }
        free(tree);
        return nf;
      }
  }
 
codetree_node_count_t codetree_num_leaves(codetree_t *tree)
  { if (tree == NULL) 
      { return 0; }
    else 
      { codetree_node_count_t nl;
        if (tree->value >= 0)
          { nl = 1; }
        else
          { nl = 0;
            nl += codetree_num_leaves(tree->child[0]);
            nl += codetree_num_leaves(tree->child[1]);
          }
        return nl;
      }
  }
  
codetree_node_t *codetree_new_leaf(codetree_data_value_t val)
  {
    demand(val >= 0, "invalid value {val}");
    demand(val <= MAX_VALUE, "invalid value {val}");
    size_t sz = (size_t)(sizeof(codetree_node_t) - 2*sizeof(codetree_node_t*));
    assert(sz >= sizeof(uint32_t));
    codetree_node_t *p = (codetree_node_t *)notnull(malloc(sz), "no mem");
    p->value = (codetree_node_value_t)val;
    return p;
  }

codetree_node_t *codetree_new_internal
  ( codetree_node_count_t seq,
    codetree_node_t *child0,
    codetree_node_t *child1
  )
  { demand((seq >= 0) && (seq < MAX_INTERNALS), "invalid internal node sequence num");
    size_t sz = sizeof(codetree_node_t);
    codetree_node_t *p = (codetree_node_t *)notnull(malloc(sz), "no mem");
    p->value = (codetree_node_value_t)(codetree_node_MIN_VALUE + (int32_t)seq);
    p->child[0] = child0;
    p->child[1] = child1;
    return p;
  }
 
codetree_bit_count_t codetree_decode
  ( codetree_byte_count_t nb, 
    byte_t buf[], 
    codetree_t *tree, 
    codetree_data_value_t maxval, 
    codetree_sample_count_t ns, 
    codetree_data_value_t smp[]
  )
  {
    bool_t debug = FALSE;
    
    if (ns == 0) { return 0; }
    demand((nb >= 0) && (nb <= MAX_BYTES), "invalid byte count {nb}");
    demand((maxval >= 1) && (maxval <= MAX_VALUE), "invalid {maxval}");
    demand(tree != NULL, "decoding tree is {NULL}");
    demand((ns >= 0) && (ns <= MAX_SAMPLES), "too many samples");

    int64_t ib = 0;    /* Index of byte of {buf} being parsed. */
    byte_t mask = 128; /* Mask of next unprocessed by of {buf[ib]}. */
    uint64_t nu = 0;   /* Number of bits used. */
    for (codetree_sample_count_t ks = 0; ks < ns; ks++)
      { /* Decode next sample {smp[ks]}: */
        codetree_node_t *p = tree;
        while (p->value < 0)
          { /* Get next bit: */
            demand(ib < nb, "buffer exhausted before {ns} samples");
            int8_t ich = ((buf[ib] & mask) != 0 ? 1 : 0);
            p = p->child[ich];
            mask >>= 1;
            if (mask == 0)
              { ib++;
                mask = 128;
              }
            demand(nu < MAX_BITS, "too many bits in encoded string");
            nu++;
          }
        codetree_data_value_t val = (codetree_data_value_t)p->value;
        if (debug) 
          { fprintf(stderr, " smp[%3lu] = %d", ks, val);
            fprintf(stderr, "\n");
          }
        demand(val <= maxval, "invalid value in tree leaf");
        smp[ks] = val;
      }
    assert(ib == nu/8);
    return nu;
  }

codetree_bit_count_t codetree_encode
  ( codetree_sample_count_t ns,
    codetree_data_value_t smp[], 
    codetree_data_value_t maxval, 
    codetree_node_count_t nd,
    codetree_delta_t delta[], 
    codetree_byte_count_t nb, 
    byte_t buf[]
  )
  {
    if (ns == 0) { return 0; }
    demand((nb >= 0) && (nb <= MAX_BYTES), "invalid byte count {nb}");
    demand((maxval >= 0) && (maxval <= MAX_VALUE), "invalid {maxval}");
    demand((ns >= 0) && (ns <= MAX_SAMPLES), "sample count {ns} too large");
    demand(nd > maxval, "table too short");
    demand(nd <= MAX_NODES, "invalid table size {nd}");

    uint64_t nb_used = 0; /* Index of next unused byte in {buf}. */
    byte_t mask = 0; /* Mask of next unused by of {buf[nb_used-1]}, or 0. */
    uint64_t nu = 0; /* Number of bits stored. */
    
    auto void enc(codetree_node_count_t id);
      /* Encode starting at {delta[id]}. */
    
    for (codetree_sample_count_t ks = 0; ks < ns; ks++)
      { /* Encode {smp[ks]}. */
        codetree_data_value_t val = smp[ks];
        demand((val >= 0) && (val <= maxval), "invalid sample value");
        codetree_node_count_t id = (uint32_t)val; /* Index into {delta} */
        enc(id);
      }
    assert(nb_used == (nu + 7)/8);
    return nu;
    
    /* -------------------------------------------------- */
    
    void enc(codetree_node_count_t id)
      { if (delta[id] != 0)
          { codetree_delta_t del = delta[id];
            int8_t ich = del % 2;
            codetree_node_count_t ip = id + del/2;
            demand((ip > id) && (ip < nd), "inconsistent {delta} table");
            enc(ip);
            if (mask == 0)
              { demand(nb_used < nb, "buffer exhausted before {ns} samples");
                nb_used++;
                buf[nb_used-1] = 0;
                mask = 128;
              }
            assert(nb_used > 0);
            if (ich != 0) { buf[nb_used-1] |= mask; }
            demand(nu < MAX_BITS, "too many bits in encoded string");
            nu++;
            mask >>= 1;
          }
      }
  }
  
codetree_node_count_t codetree_get_encoding_table
  ( codetree_t *tree, 
    codetree_data_value_t maxval, 
    codetree_node_count_t nd,
    codetree_delta_t delta[]
  )
  { 
    bool_t debug = FALSE;
    if (tree == NULL) { return 0; }
    demand((maxval >= 0) && (maxval <= MAX_VALUE), "invalid maxval");
    demand((nd > maxval) && (nd <= MAX_NODES), "invalid table size {nd}");
    if (debug) { fprintf(stderr, "    clearing table (nd = %u)...\n", nd); }
    for (codetree_node_count_t id = 0; id < nd; id++) { delta[id]= 0; }
    
    codetree_node_count_t nvalid = 0; /* Number of leaves of the tree. */
    codetree_node_count_t nint = 0;   /* Number of internal nodes. */
    
    auto codetree_node_count_t trans(codetree_node_t *q, int32_t level);
      /* Sets the pointers {delta[j]} for all nodes in the 
        subtree of {q}, except {q} itself.  Returns the index
        {iq} of the {delta} entry assigned to {q}, but leaves 
        {delta[iq] = 0}.  Also increments {nvalid} and {nint}. */
    
    if (debug) { fprintf(stderr, "    recursing on tree...\n"); }
    codetree_node_count_t ir = trans(tree, 0);
    
    assert((nvalid > 0) && (nvalid <= maxval+1));
    assert(nint == nvalid - 1);
    codetree_node_count_t md = maxval + 1 + nint; /* {delta} entries actually used. */
    assert(md <= nd);
    assert((ir <= maxval) || (ir == md-1));
    return md;
    
    /* -------------------------------------------------- */
    
    codetree_node_count_t trans(codetree_node_t *q, int32_t level)
      { codetree_node_value_t nval = q->value;
        codetree_node_count_t id; /* Index in {delta} assigned to {q}. */
        if (nval >= 0)
          { codetree_data_value_t  val = (codetree_data_value_t)nval;
            if (debug) { fprintf(stderr, "    | %*sval = %d\n", 2*level, "", val); }
            demand(val <= maxval, "invalid value in tree leaf");
            demand(nvalid < MAX_LEAVES, "too many leaves in tree");
            nvalid++;
            /* The {delta} entry index for a leaf is its value: */
            id = (codetree_node_count_t)val;
          }
        else
          { codetree_node_count_t seq = (codetree_node_count_t)(nval - codetree_node_MIN_VALUE); /* Sequence number in tree.*/
            if (debug) { fprintf(stderr, "    | %*sseq = %d\n", 2*level, "", seq); }
            /* Setup table for subtrees: */
            codetree_node_count_t id_ch[2];
            for (int8_t ich = 0; ich < 2; ich++) 
              { id_ch[ich] = trans(q->child[ich], level+1); }

            /* Assign index {id} for this node: */
            id = (codetree_node_count_t)(maxval + 1 + nint);
            if (debug) { fprintf(stderr, "    | %*sseq = %d  id = %u\n", 2*level, "", seq, id); }
            nint++;

            /* Set parent pointers of children: */
            demand(id_ch[0] != id_ch[1], "repeated value in table");
            for (uint8_t ich = 0; ich < 2; ich++) 
              { assert(id_ch[ich] < id);
                assert(delta[id_ch[ich]] == 0);
                delta[id_ch[ich]] = 2*(id - id_ch[ich]) + ich;
                if (debug) { fprintf(stderr, "    | %*sdelta[%3u] = %3u:%d\n", 2*level, "", id_ch[ich], id, ich); }
              }
          }
        demand(id < nd, "not enough space in table");
        demand(delta[id] == 0, "repeated leaf value in tree");
        return id;
      }
  }

codetree_t *codetree_get_decoding_tree
  ( codetree_node_count_t nd,
    codetree_delta_t delta[],
    codetree_data_value_t maxval
  )
  {
    if (nd == 0) { return NULL; }
    demand((maxval >= 0) && (maxval <= MAX_VALUE), "invalid maxval");
    demand((nd >= maxval) && (nd <= MAX_NODES), "invalid table size {nd}");
    
    /* Count valid values: */
    codetree_node_count_t nvalid = 0; /* Number of valid values. */
    codetree_data_value_t aval; /* Some valid value. */ 
    for (codetree_node_count_t id = 0; id <= maxval; id++) 
      { if (delta[id] != 0) { aval = (codetree_data_value_t)id; nvalid++; } }
    
    /* If there are no leaf nodes, the tree is {NULL} and no samples can be encoded. */
    if (nvalid == 0) { return NULL; }
    
    /* If there is only one valid value, the tree is just that leaf node: */
    if (nvalid == 1) { return codetree_new_leaf(aval); }

    /* We have some internal nodes: */
    codetree_node_count_t nint = nvalid-1;
    codetree_node_count_t md = maxval + 1 + nint; /* Expected number of table entries. */
    demand(nd >= md, "table is too short");
    
    /* Allocate node pointer table: */
    codetree_node_t **nodes = (codetree_node_t **)notnull(malloc(md*sizeof(codetree_node_t *)), "no mem");
    
    /* Create the leaf nodes: */
    for (codetree_node_count_t id = 0; id <= maxval; id++)
      { if (delta[id] != 0)
          { nodes[id] = codetree_new_leaf((codetree_data_value_t)id); }
        else
          { nodes[id] = NULL; }
      }
    
    /* Create the internal nodes, still without children: */
    for (codetree_node_count_t ii = 0; ii < nint; ii++)
      { codetree_node_count_t id = maxval + 1 + ii;
        nodes[id] = codetree_new_internal(ii, NULL, NULL);
      }
      
    /* Set children of internal nodes based on {delta}: */
    for (codetree_node_count_t id = 0; id < md; id++)
      { if (delta[id] != 0)
          { codetree_node_t *q = nodes[id];
            assert(q != NULL);
            int8_t ich = delta[id] % 2;
            codetree_node_count_t ip = id + delta[id]/2;
            demand(ip > id, "table order bug");
            demand((ip > maxval) && (ip < md), "bad table");
            codetree_node_t *p = nodes[ip]; /* Parent node of {q}. */
            demand(p->child[ich] == NULL, "bad table");
            p->child[ich] = q;
          }
        else
          { assert((id <= maxval) || (id == md-1)); }
      }
    
    /* Get the root: */
    codetree_node_t *root = nodes[md-1];
    free(nodes);
    return root;
  }

void codetree_check_iso(codetree_t *p, codetree_t *q)
  { demand((p == NULL) == (q == NULL), "only one is null");
    if (p == NULL) { return; }
    demand((p->value >= 0) == (q->value >= 0), "only one is leaf");
    if (p->value >= 0)
      { demand(p->value == q->value, "leaves have different values"); }
    else
      { for (int8_t ich = 0; ich < 2; ich++)
          { demand(p->child[ich] != NULL, "null child in {p}");
            demand(q->child[ich] != NULL, "null child in {q}");
            codetree_check_iso(p->child[ich], q->child[ich]);
          }
      }
  }

void codetree_check_tree(codetree_t *tree, codetree_data_value_t maxval)
  { 
    fprintf(stderr, "\n");
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, "--- begin {codetree_check_tree} - {maxval} = %u\n", maxval);
        
    demand((maxval >= 0) && (maxval <= MAX_VALUE), "bad {mxval}");
    
    fprintf(stderr, "  counting the leaves...\n");
    codetree_node_count_t nvalid = codetree_num_leaves(tree);
    fprintf(stderr, "  found %u leaves\n", nvalid);
    
    if (nvalid == 0)
      { demand(tree == NULL, "non-empty tree with empty set {V}"); }
    else if (nvalid == 1)
      { demand((tree != NULL) && (tree->value >= 0), "tree for single {V} must be a leaf");
        codetree_data_value_t val = (codetree_data_value_t)tree->value;
        demand(val <= maxval, "invalid value in leaf node");
      }
    else
      { fprintf(stderr, "  building the encoding table from the tree...\n");
        codetree_node_count_t nd = maxval + nvalid;
        demand(nd <= MAX_NODES, "too many nodes in tree");
        codetree_delta_t *delta = (codetree_delta_t*)notnull(malloc(nd*sizeof(codetree_delta_t)), "no mem");
        codetree_get_encoding_table(tree, maxval, nd, delta);
        
        fprintf(stderr, "  rebuilding the tree from the encoding table...\n");
        codetree_t *tree2 = codetree_get_decoding_tree(nd, delta, maxval);
        
        fprintf(stderr, "  comparing the two trees...\n");
        codetree_check_iso(tree, tree2);
        
        free(delta);
        codetree_free(tree2);
      }
    fprintf(stderr, "--- end {codetree_check_tree}\n");
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, "\n");
  }
          
void codetree_check_table(codetree_node_count_t nd, codetree_delta_t delta[], codetree_data_value_t maxval)
  {
    fprintf(stderr, "\n");
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, "--- begin {codetree_check_table} - {nd} = %u {maxval} = %u\n", nd, maxval);

    demand(maxval <= MAX_VALUE, "bad {mxval}");
    demand((nd >= maxval) && (nd <= MAX_NODES), "invalid table size {nd}");
    
    fprintf(stderr, "  obtaining the decoding tree from the table...\n");
    codetree_t *tree = codetree_get_decoding_tree(nd, delta, maxval);
    
    fprintf(stderr, "  getting back the table from the tree...\n");
    codetree_delta_t *delta2 = (codetree_delta_t*)notnull(malloc(nd*sizeof(codetree_delta_t)), "no mem");
    codetree_node_count_t md = codetree_get_encoding_table(tree, maxval, nd, delta);
    assert(md <= nd);
    
    fprintf(stderr, "  comparing the two tables...\n");
    for(codetree_data_value_t val = 0; val <= maxval; val++)
      { demand((delta[val] == 0) == (delta2[val] == 0), "only one of {delta,delta2} is zero");
        if (delta[val] != 0)
          { codetree_node_count_t id1 = (uint32_t)val;
            codetree_node_count_t id2 = (uint32_t)val;
            while ((id1 < md-1) || (id2 < md-1))
              { demand((id1 < md-1) == (id2 < md-1), "only one reached the root");
                int8_t ich1 = delta[id1] % 2;
                int8_t ich2 = delta2[id2] %2;
                demand(ich1 == ich2, "codes differ");
                codetree_node_count_t ip1 = id1 + delta[id1] / 2;
                demand(ip1 > id1, "{delta} not monotonic");
                demand((ip1 > maxval) && (ip1 < md), "{delta} buggy");
                codetree_node_count_t ip2 = id2 + delta2[id2] / 2;
                demand(ip2 > id2, "{delta2} not monotonic");
                demand((ip2 > maxval) && (ip2 < md), "{delta2} buggy");
                id1 = ip1;
                id2 = ip2;
              }
            demand(delta[id1] == 0, "bad root entry in {delta}");
            demand(delta2[id2] == 0, "bad root entry in {delta2}");
          }
      }
    codetree_free(tree);
    free(delta2);
    fprintf(stderr, "--- end {codetree_check_table}\n");
    fprintf(stderr, "----------------------------------------------------------------------\n");
    fprintf(stderr, "\n");
  }

codetree_node_count_t codetree_print_codes(FILE *wr, codetree_t *tree)
  { 
    if (tree == NULL) { fprintf(wr, "empty tree\n"); return 0; }

    uint32_t nv = 0;
    
    auto void enum_leaves(char *pref, codetree_node_t *p);
      /* Prints the values and codes of the subtree rooted at {p},
        assuming that its code is the string {pref}. Also increments {nv}
        for each leaf found. */
        
    enum_leaves("", tree);
    
    return nv;
    
    void enum_leaves(char *pref, codetree_node_t *p)
      { if (p->value < 0)
          { for (int8_t ich = 0; ich < 2; ich++)
              { char *pref_ch = jsprintf("%s%d", pref, ich);
                enum_leaves(pref_ch, p->child[ich]);
                free(pref_ch);
              }
          }
        else
          { fprintf(stderr, "%12d (%s)\n", p->value, pref);
            nv++;
          }
      }
  }

codetree_node_count_t codetree_check_codes(codetree_t *tree, codetree_data_value_t maxval, char *code[])
  { /* Check if the codes lead to the desired values, note max code length {maxlen}: */
    uint32_t ncode = 0; /* Number of non-null codes. */
    for (codetree_data_value_t val = 0; val <= maxval; val++)
      { char *ck = code[val];
        if (ck != NULL)
          { ncode++;
            /* Check if {code[val]} leads to leaf value {val}: */
            codetree_node_t *p = tree;
            while ((p->value < 0) && ((*ck) != 0))
              { demand(p->child[0] != NULL, "null child 0 in internal node");
                demand(p->child[1] != NULL, "null child 1 in internal node");
                if ((*ck) == '0')
                  { p = p->child[0]; }
                else if ((*ck) == '1')
                  { p = p->child[1]; }
                else
                  { assert(FALSE); }
                ck++;
              }
            demand(p->value >= 0, "code leads to non-leaf");
            demand((*ck) == 0, "code tries to go below leaf");
            demand(p->value == val, "leaf value does not check");
          }
      }
    demand(ncode == codetree_num_leaves(tree), "{code} count does not match leaf count");
    if (ncode == 0) { assert(tree == NULL); }
    return ncode;
  }

void codetree_print_bits(FILE *wr, codetree_byte_count_t nb, byte_t buf[], char *sep)
  { 
    if (nb ==  0) { return; }
    int64_t ib = 0;    /* Byte currently being printed. */
    byte_t mask = 128; /* Mask of next unprocessed by of {buf[ib]}, or 0. */
    while (TRUE)
      { char ch = ((buf[ib] & mask) != 0 ? '1' : '0');
        fputc(ch, wr);
        mask >>= 1;
        if (mask == 0)
          { ib++;
            if (ib >= nb) { break; }
            if (sep != NULL) { fputs(sep, wr); }
            mask = 128;
          }
      }     
    fflush(wr);
  }
