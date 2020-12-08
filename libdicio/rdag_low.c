/* Procedures of {rdag.h} that require access to the internal rep. */
/* Last edited on 2017-01-03 02:04:15 by jstolfi */

#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include <affirm.h>
#include <bool.h>

#include <rdag.h>
#include <rdag_def.h>

/* 
  INTERNAL PROTOTYPES */
  
uint32_t rdag_hash_table_size(rdag_node_t max_alloc_node);
  /* Computes a suitable hash table size for a dag with nodes
    numbered from 0 to {max_alloc_node}. */

/*
  INTERNAL PROCEDURES

  These procedures require access to {rdag_def.h} */

uint32_t extract64(uint64_t pk, uint32_t shift, uint32_t mask);
  /* Extracts the bitfield of {pk} (up to 32 bits) with the given
    {shift} and {mask}. */

uint64_t insert64(uint64_t pk, uint32_t val, uint32_t shift, uint32_t mask);
  /* Inserts the value {val} into the bitfield of {pk} (up to 32 bits) with
    given {shift} and {mask}. */

uint32_t extract64(uint64_t pk, uint32_t shift, uint32_t mask)
  {
    return mask & (uint32_t)(pk >> shift);
  }

uint64_t insert64(uint64_t pk, uint32_t val, uint32_t shift, uint32_t mask)
  {
    return (pk & (~ (((uint64_t)mask) << shift))) | (((uint64_t)val) << shift);
  }

void rdag_node_unpack(rdag_t *D, rdag_node_packed_t *pk, rdag_node_data_t *dt)
  {
    dt->f_link = extract64((*pk), D->shift_f_link, D->mask_node);
    dt->i_mark = extract64((*pk), D->shift_i_mark, D->mask_i_symbol);
    dt->o_mark = extract64((*pk), D->shift_o_mark, D->mask_o_symbol);
    dt->p_link = extract64((*pk), D->shift_p_link, D->mask_node);
  }
  
void rdag_node_pack(rdag_t *D, rdag_node_data_t *dt, rdag_node_packed_t *pk)
  {
    rdag_node_packed_t pkd = 0;
    pkd = insert64(pkd, dt->f_link, D->shift_f_link, D->mask_node);
    pkd = insert64(pkd, dt->i_mark, D->shift_i_mark, D->mask_i_symbol);
    pkd = insert64(pkd, dt->o_mark, D->shift_o_mark, D->mask_o_symbol);
    pkd = insert64(pkd, dt->p_link, D->shift_p_link, D->mask_node);
    (*pk) = pkd;
  }

/*
  LOW-LEVEL PROCEDURES

  These procedures require access to {rdag_def.h} */

rdag_node_t rdag_node_max(rdag_t *D)
  {
    return D->max_node;
  }

rdag_symbol_t rdag_i_mark_max(rdag_t *D)
  {
    return (rdag_symbol_t)(D->mask_i_symbol);
  }
  
rdag_symbol_t rdag_o_mark_max(rdag_t *D)
  {
    return (rdag_symbol_t)(D->mask_o_symbol);
  }

rdag_symbol_t rdag_i_mark(rdag_t *D, rdag_node_t s)
  {
    assert(s <= D->max_node);
    assert(s > 0);
    rdag_node_packed_t pk = D->pk[s-1];
    return (rdag_symbol_t)extract64(pk, D->shift_i_mark, D->mask_i_symbol);
  }

rdag_symbol_t rdag_o_mark(rdag_t *D, rdag_node_t s)
  {
    assert(s <= D->max_node);
    assert(s > 0);
    rdag_node_packed_t pk = D->pk[s-1];
    return (rdag_symbol_t)extract64(pk, D->shift_o_mark, D->mask_o_symbol);
  }

rdag_node_t rdag_f_link(rdag_t *D, rdag_node_t s)
  {
    assert(s <= D->max_node);
    assert(s > 0);
    rdag_node_packed_t pk = D->pk[s-1];
    return (rdag_node_t)extract64(pk, D->shift_f_link, D->mask_node);
  }

rdag_node_t rdag_p_link(rdag_t *D, rdag_node_t s)
  {
    assert(s <= D->max_node);
    assert(s > 0);
    rdag_node_packed_t pk = D->pk[s-1];
    return (rdag_node_t)extract64(pk, D->shift_p_link, D->mask_node);
  }

void rdag_node_data_get(rdag_t *D, rdag_node_t s, rdag_node_data_t *dt)
  {
    demand(s > 0, "null node has no data");
    demand(s <= rdag_node_max(D), "no such node");
    rdag_node_unpack(D, &(D->pk[s-1]), dt);
  }

void rdag_node_data_set(rdag_t *D, rdag_node_t s, rdag_node_data_t *dt)
  {
    demand(s > 0, "null node has no data");
    demand(s <= rdag_node_max(D), "no such node");
    demand(dt->f_link < s, "invalid f-link");
    demand(dt->i_mark <= rdag_i_mark_max(D), "invalid i-mark");
    demand(dt->o_mark <= rdag_o_mark_max(D), "invalid o-mark");
    demand(dt->p_link < s, "invalid p-link");
    rdag_node_pack(D, dt, &(D->pk[s-1]));
    D->hash_valid = FALSE;
  }

void rdag_discard(rdag_t *D, rdag_node_t s)
  {
    demand(s > 0, "cannot discard the null node");
    if (s <= D->max_node)
      { D->max_node = s-1;
        D->hash_valid = FALSE;
      }
  }

void rdag_check_invariants(rdag_t *D, rdag_node_t lo, rdag_node_t hi)
  {
    /* If there are no nodes in the interval, then OK: */
    if (lo > D->max_node) { return; }
    
    /* Check the nodes in the interval: */
    rdag_node_t s;
    for (s = lo; (s < hi) && (s <= D->max_node); s++) 
      {
        if (s == rdag_node_NULL)
          { /* The null node is always OK: */ }
        else 
          { /* Get the node attributes: */
            rdag_node_t f = rdag_f_link(D, s);
            rdag_node_t p = rdag_p_link(D, s);
            /* Check topological ordering of nodes: */
            demand(f < s, "invalid f-link");
            demand(p < s, "invalid p-link");
            /* Check lexicographic ordering of subnodes: */
            if (f != rdag_node_NULL)
              { rdag_symbol_t i_s = rdag_i_mark(D, s);
                rdag_symbol_t i_f = rdag_i_mark(D, f);
                demand(i_s > i_f, "i-marks repeated or out of order");
              }
          }
      }
  }

uint32_t rdag_node_hash_packed(rdag_t *D, rdag_node_packed_t *pk)
  { 
    rdag_node_data_t dt;
    rdag_node_unpack(D, pk, &dt);
    return rdag_node_hash(D, &dt);
  }

uint32_t rdag_node_hash(rdag_t *D, rdag_node_data_t *dt)
  {
    /* This is the original "dicio" hash function (J.Stolfi, 1992). */
    /* The following are standard random numbers: 8-) */
    uint64_t M_i = 418 - 1;
    uint64_t M_o = 1003;
    uint64_t M_f = 4615;
    uint64_t M_p = 480317;
    /* Extract the fields and multiply: */
    uint64_t H_i = M_i * (uint64_t)dt->i_mark;
    uint64_t H_o = M_o * (uint64_t)dt->o_mark;
    uint64_t H_f = M_f * (uint64_t)dt->f_link;
    uint64_t H_p = M_p * (uint64_t)dt->p_link;
    uint64_t modulus = (uint64_t)(D->hash_size);
    return (uint32_t)((H_i + H_o + H_f + H_p) % modulus);
  }

rdag_node_t rdag_node_max_alloc(rdag_t *D)
  {
    return D->max_alloc_node;
  }

uint32_t rdag_hash_table_size(rdag_node_t max_alloc_node)
  {
    return (max_alloc_node / 2) | 1UL;
  }

void rdag_rehash(rdag_t *D)
  {
    /* Reallocate the hash tables of proper size: */
    uint32_t hash_size = rdag_hash_table_size(D->max_alloc_node);
    if (hash_size != D->hash_size)
      { if (D->debug) { fprintf(stderr, "%s: new hash size = %u\n", __FUNCTION__, hash_size); }
        D->hinit = (rdag_node_t *)notnull(realloc(D->hinit, hash_size * sizeof(rdag_node_t)), "no mem");
        D->hnext = (rdag_node_t *)notnull(realloc(D->hnext, D->max_alloc_node * sizeof(rdag_node_t)), "no mem");
        D->hash_size = hash_size;
      }

    /* Initialize all buckets to the empty list: */
    uint32_t k;
    for (k = 0; k < D->hash_size; k++) { D->hinit[k] = rdag_node_NULL; }
    
    /* Rehash all nodes and insert in hash table: */
    for (k = 0; k < D->max_node; k++)
      { rdag_node_t s = k+1;
        uint32_t h = rdag_node_hash_packed(D, &(D->pk[s-1]));
        D->hnext[k] = D->hinit[h];
        D->hinit[h] = s;
      }

    /* Now the hash tables are valid: */
    D->hash_valid = TRUE;
  }
  
#define rdag_ONE32 ((uint32_t)1)

rdag_t *rdag_new(uint32_t nn, uint32_t ni, uint32_t no, rdag_node_t max_alloc_node)
  {
    demand(nn <= rdag_nn_MAX, "too many bits in node number");
    demand(ni <= rdag_ni_MAX, "too many bits in input symbol");
    demand(no <= rdag_no_MAX, "too many bits in output symbol");
    demand(ni + no + 2*nn <=rdag_ntot_MAX, "too many bits per node");
    
    rdag_t *D = (rdag_t *)notnull(malloc(sizeof(rdag_t)), "no mem");
    D->debug = TRUE; /* For now. */
    D->nn = nn;
    D->ni = ni;
    D->no = no;
    
    /* Create the field masks: */
    D->mask_node = (rdag_ONE32 << nn) - rdag_ONE32;
    D->mask_i_symbol = (rdag_ONE32 << ni) - rdag_ONE32;
    D->mask_o_symbol = (rdag_ONE32 << no) - rdag_ONE32;
    
    /* Compute the field shifts: */
    D->shift_f_link = 0;
    D->shift_i_mark = D->shift_f_link + nn;
    D->shift_o_mark = D->shift_i_mark + ni;
    D->shift_p_link = D->shift_o_mark + no;
    
    D->max_node = 0; 

    /* Null data vector, to be allocated by {rdag-expand}: */
    D->max_alloc_node = 0; 
    D->pk = NULL;
    
    /* Null hash table, to be created by {rdag_rehash}: */
    D->hash_size = 0;
    D->hinit = NULL;
    D->hnext = NULL;
    D->hash_valid = FALSE;
    
    /* Now allocate to requested size: */
    rdag_expand(D, max_alloc_node);
    
    return D;
  }

rdag_node_t rdag_node_find_packed(rdag_t *D, rdag_node_packed_t *pk, uint32_t h)
  {
    /* Make sure that the hash tables are valid: */
    if (! D->hash_valid) { rdag_rehash(D); }

    /* Look it up in the hash table: */
    rdag_node_t *hinit = D->hinit;
    rdag_node_t *hnext = D->hnext;
    rdag_node_packed_t *data = D->pk;
    rdag_node_t s = hinit[h];
    rdag_node_t s1 = rdag_node_NULL; /* Node before {s} in bucket. */
    rdag_node_t s2 = rdag_node_NULL; /* Node before {s1} in bucket. */
    if (s != rdag_node_NULL)
      { assert(s <= D->max_node);
        while ((s != rdag_node_NULL) && (data[s-1] != (*pk)))
          { s2 = s1; s1 = s; s = hnext[s-1];assert(s <= D->max_node); }
      }
    if (s != rdag_node_NULL)
      { /* Node already exist ({s}): */
        if (s1 != rdag_node_NULL)
          { /* Not first in bucket, move it up: */
            /* (Is this really helpful?) */
            if (s2 == rdag_node_NULL)
              { hinit[h] = s; }
            else
              { hnext[s2-1] = s; }
            hnext[s1-1] = hnext[s-1];
            hnext[s-1] = s1;
          }
      }
    return s;
  }

rdag_node_t rdag_node_append_packed(rdag_t *D, rdag_node_packed_t *pk, uint32_t h)
  {
    /* Get the next available node number: */
    demand (D->max_node < D->mask_node, "too many nodes");
    rdag_node_t s = D->max_node + 1;
    assert(s <= D->mask_node);
    /* No luck, must create a new entry: */
    if (s > D->max_alloc_node)
      { /* More bad luck, ran out of space: */
        /* Choose new max allocated node: */
        rdag_node_t old_max_alloc_node = D->max_alloc_node;
        rdag_node_t new_max_alloc_node;
        rdag_node_t max_possible_node = D->mask_node;
        if (old_max_alloc_node < (max_possible_node - 1)/2)
          { new_max_alloc_node = 2*old_max_alloc_node + 1; }
        else
          { new_max_alloc_node = max_possible_node; }
        /* This should be true in either case: */
        assert(s <= new_max_alloc_node);
        rdag_expand(D, new_max_alloc_node);
      }
    /* Store its data: */
    assert(s <= D->max_alloc_node);
    D->pk[s-1] = (*pk);
    D->max_node = s;
    if (D->hash_valid)
      { /* Add it to the hash bucket: */
        assert(h < D->hash_size);
        D->hnext[s-1] = D->hinit[h];
        D->hinit[h] = s;
      }
    return s;
  }

void rdag_expand(rdag_t *D, rdag_node_t max_alloc_node)
  {
    /* Check number of bits in node number: */
    demand(max_alloc_node <= D->mask_node, "max_alloc_node too large");

    uint32_t old_size = D->max_alloc_node;
    uint32_t new_size = max_alloc_node;
    
    /* Check whether expansion is really necessary: */
    if (new_size <= old_size){ return; }
    
    /* Expand the old node data vector: */
    if (D->debug) { fprintf(stderr, "%s: old size = %u  new size = %u\n", __FUNCTION__, old_size, new_size); }
    D->pk = (rdag_node_packed_t *)notnull(realloc(D->pk, new_size * sizeof(rdag_node_packed_t)), "no mem");
    D->max_alloc_node = max_alloc_node;
    
    /* The hash tables are no longer valid: */
    D->hash_valid = FALSE;
  }

uint32_t rdag_reachable_node_count(rdag_t *D, uint32_t nroots, rdag_node_t root[])
  {
    /* Find the maximum root node {m}: */
    rdag_node_t max_root = 0;
    uint32_t i;
    for (i = 0; i < nroots; i++) 
      { rdag_node_t ri = root[i];
        demand(ri <= D->max_node, "invalid root node");
        if (ri > max_root) { max_root = ri; }
      }
   
    /* Allocate a boolean table {rch} such that {rch[s-1]} tells whether {s} is reachable: */
    bool_t *rch = (bool_t *)notnull(malloc(max_root*sizeof(bool_t)), "no mem");
    
    /* Mark all nodes unreachable: */
    uint32_t k;
    for (k = 0; k < max_root; k++) { rch[k] = FALSE; }

    /* Mark all proper roots as reachable: */
    for (i = 0; i < nroots; i++) 
      { rdag_node_t ri = root[i];
        assert(ri <= max_root);
        if (ri > 0) { rch[ri-1] = TRUE; }
      }
   
    /* Scan nodes propagating their reachability and counting: */
    uint32_t nreach = 0; /* Count of reachable nodes. */
    rdag_node_t s = max_root;
    while (s > 0)
      { if (rch[s-1])
          { /* One more reachable node: */
            nreach++;
            /* Get its fields: */
            rdag_node_data_t dt;
            rdag_node_data_get(D, s, &dt);
            /* Propagate through its f-link: */
            rdag_node_t f = dt.f_link;
            assert(f < s);
            if (f > 0) { rch[f-1] = TRUE; }
            /* Propagate through its p-link: */
            rdag_node_t p = dt.p_link;
            assert(p < s);
            if (p > 0) { rch[p-1] = TRUE; }
          }
        s--;
      }
    /* Free the table: */
    free(rch);
    
    return nreach;
  }

void rdag_crunch(rdag_t *D, uint32_t nroots, rdag_node_t root[])
  {
    /*
      rdag_crunch is done in in four passes:

      rdag_first pass: unmark all nodes ("mark[s-1] = FALSE").

      Second pass: mark all roots as reachable ("mark[s-1] = TRUE"),
      then scan the nodes from high to low, mark the descendants of
      every marked node found.  (This works because the dag is
      acyclic and the nodes are topologically sorted).

      Third pass: scan all nodes from 1 up, and move each marked
      node {s} to the lowest unused position, updating its {dest} and
      {rest} pointers, and saving that position in {map[s-1]}.  Also
      update the root pointers.

      Fourth pass: rebuild the hash table for the compacted nodes.

      We actually store the marks {mark[s-1]} and the pointers
      {map[s-1]} in {hnext[s-1]}, since the hash table will be
      rebuilt anyway in step four.
    */

    /* Mark the hash tables as invalid: */
    D->hash_valid = FALSE;
    
    /* Unmark all nodes */
    rdag_node_t *mark = D->hnext; /* Using {D.hnext[]} to mark nodes. */
    uint32_t k;
    for (k = 0; k < D->max_node; k++) { mark[k] = rdag_node_NULL; }


    /* Mark root nodes, and remember highest one: */
    /* (convention: a node is marked iff it has {mark[t-1] == 1}.) */
    rdag_node_t old_max_node = 0;
    uint32_t i;
    for (i = 0; i < nroots; i++)
      { rdag_node_t t = root[i];
        if (t != rdag_node_NULL)
          { if (t > old_max_node) { old_max_node =  t; }
            mark[t-1] = 1;
          }
      }

    /* Scan from {old_max_node} down to {1}, marking the descendants of marked nodes: */
    /* Note that {rdag_node_t} is unsigned so we can't count down to 0 with a {for}. */
    uint32_t s;
    for (s = old_max_node; s > 0; s--)
      {
        if (mark[s-1] == 1)
          { /* Node {s} is marked, mark its sucessors: */
            rdag_node_t f = rdag_f_link(D, s);
            if (f != rdag_node_NULL) { mark[f-1] = 1; }
            rdag_node_t p = rdag_p_link(D, s);
            if (p != rdag_node_NULL) { mark[p-1] = 1; }
          }
      }

    /* Move every reachable node {t} to the lowest free place, record it {map[t-1]}: */
    rdag_node_t *map = D->hnext; /* Using {D.hnext[]} as map. */
    
    auto rdag_node_t apply_map(rdag_node_t old_t);
      /* Gets the new number of node {old_s}, from the map. 
        If {old-s} is null, returns it unchanged. 
        Also does some checks. */
        
    rdag_node_t apply_map(rdag_node_t old_t)
      { if (old_t == rdag_node_NULL) { return rdag_node_NULL; }
        rdag_node_t new_t = map[old_t - 1];
        assert(new_t != rdag_node_NULL);
        assert(new_t <= old_t);
        return new_t;
      }
      
    rdag_node_t new_max_node = rdag_node_NULL;
    for (k = 0; k < old_max_node; k++)
      {
        rdag_node_t old_s = k+1;
        if (mark[old_s-1] == 1)
          { /* Node {old_s} is marked, pull it down to {new_s}: */
            new_max_node++;
            rdag_node_t new_s = new_max_node;
            assert(new_s <= old_s);
            if (new_s < old_s)
              { /* Successors may have changed too: */
                /* Get old node data: */
                rdag_node_packed_t *old_pk = &(D->pk[old_s-1]);
                rdag_node_data_t old_dt;
                rdag_node_unpack(D, old_pk, &old_dt);
                /* Update the f-link and p-link: */
                rdag_node_data_t new_dt;
                new_dt.f_link = apply_map(old_dt.f_link);
                new_dt.i_mark = old_dt.i_mark;
                new_dt.o_mark = old_dt.o_mark;
                new_dt.p_link = apply_map(old_dt.p_link);
                /* Store the new data under the new node number: */
                rdag_node_packed_t new_pk;
                rdag_node_pack(D, &new_dt, &new_pk);
                D->pk[new_s-1] = new_pk;
              }
            /* Record the new node number in the map: */
            map[old_s-1] = new_s;
          }
      }

    /* Update the root nodes and {D->max_node}: */
    for (i = 0; i < nroots; i++) { root[i] = apply_map(root[i]); };
    D->max_node = new_max_node;

    /* Rehash nodes */
    rdag_rehash(D);
  }

rdag_node_t rdag_node_from_fields(rdag_t *D, rdag_node_t f, rdag_symbol_t i, rdag_symbol_t o, rdag_node_t p)
  { 
    rdag_node_data_t dt = (rdag_node_data_t){ .f_link = f, .i_mark = i, .o_mark = o, .p_link = p };
    return rdag_node_from_data(D, &dt);
  }

rdag_node_t rdag_node_from_data(rdag_t *D, rdag_node_data_t *dt)
  {
    if (dt->i_mark > rdag_i_mark_max(D)) { fprintf(stderr, "%s: i_mark = %u", __FUNCTION__, dt->i_mark); }
    
    demand(dt->f_link <= rdag_node_max(D), "invalid f-link");
    demand(dt->i_mark <= rdag_i_mark_max(D), "invalid i-mark");
    demand(dt->o_mark <= rdag_o_mark_max(D), "invalid o-mark");
    demand(dt->p_link <= rdag_node_max(D), "invalid p-link");
    
    /* Make sure that we have the hash tables: */
    if (! D->hash_valid) { rdag_rehash(D); }
    
    /* See if node already exists: */
    rdag_node_packed_t pk; rdag_node_pack(D, dt, &pk);
    uint32_t h = rdag_node_hash(D, dt);
    rdag_node_t s = rdag_node_find_packed(D, &pk, h);
    if (s == rdag_node_NULL)
      { /* Append as a new node: */
        s = rdag_node_append_packed(D, &pk, h);
      }
    return s;
  }

void rdag_free(rdag_t *D)
  {
    free(D->pk);
    free(D->hinit);
    free(D->hnext);
    free(D);
  }
