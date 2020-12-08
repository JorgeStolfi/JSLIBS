#ifndef rdag_def_H
#define rdag_def_H

#define rdag_def_DESC "Private definitions of {rdag_t}"

/* Last edited on 2009-10-29 13:07:45 by stolfi */
/*©*/

#include <stdint.h>

#include <bool.h>

#include <rdag.h>

typedef uint64_t rdag_node_packed_t;
  /* Packed node data. */

struct rdag_t 
  { 
    /* Basic data: */
    uint32_t nn; /* Number of bits per node field. */
    uint32_t ni; /* Number of bits per input mark. */
    uint32_t no; /* Number of bits per output mark. */
    
    /* Node data: */
    rdag_node_t max_node;       /* Max node currently in dag. */
    rdag_node_t max_alloc_node; /* Max node that can be stored without expansion. */
    rdag_node_packed_t *pk;     /* {pk[s-1]} is the packed data of proper node {s}. */
    
    /* Node hash table to find nodes with given data (only for proper nodes): */
    uint32_t hash_size;  /* Size of hash table. */
    rdag_node_t *hinit;  /* {hinit[h]} is the first node with hash value {h}. */
    rdag_node_t *hnext;  /* {hnext[s-1]} is the next node with same hash value as {s}. */
    bool_t hash_valid;   /* TRUE iff the contents of {hinit} and {hnext} are valid. */

    /* Bit shifts for each field: */
    uint32_t shift_f_link;  /* Shift of f-link field. */
    uint32_t shift_p_link;  /* Shift of p-link field. */
    uint32_t shift_i_mark; /* Shift of i-mark field. */
    uint32_t shift_o_mark; /* Shift of o-mark field. */
    
    /* Masks for each field type: */
    uint32_t mask_node;     /* Mask for a node number. */
    uint32_t mask_i_symbol; /* Mask for an i-symbol. */
    uint32_t mask_o_symbol; /* Mask for an o-symbol. */
    
    /* Debugging info: */
    bool_t debug;  /* Set to TRUE to print diagnostics. */
    
  };
  /* The representation of an automaton. 
  
    The {.pk} and {.hnext} vectors have {max_alloc_node} entries.
    The {.hinit} vector has {.hash_size} entries.
    
    The valid nodes are {0 .. .max_node}. The attributes of proper
    node {s} are stored in {pk[s-1]}, packed into a 64-bit word. */

/* LOW-LEVEL PROCEDURES */

void rdag_node_unpack(rdag_t *D, rdag_node_packed_t *pk, rdag_node_data_t *dt);
  /* Loosens the packed node data {*pk} into the unpacked data record {*dt}. */

void rdag_node_pack(rdag_t *D, rdag_node_data_t *dt, rdag_node_packed_t *pk);
  /* Squeezes the node data record {*dt} into the packed data {*pk}. */

rdag_node_t rdag_node_find_packed(rdag_t *D, rdag_node_packed_t *pk, uint32_t h);
  /* Returns a node in {D} with the packed node data {*pk} and hash
    value {h}. If there is no such node, return {rdag_node_SKIP}.
    Rebuilds the hash tables if they are not valid.
    
    Cost: if the hash tables are valid, {O(1)} time and space. If the
    hash tables are not valid, {O(rdag_max_node(D))} time and
    space. */

rdag_node_t rdag_node_append_packed(rdag_t *D, rdag_node_packed_t *pk, uint32_t h);
  /* Appends to {D} a node with the packed node data {*pk} and hash
    value {h}, and returns its number. May expand {D} if necessary,
    usually to about twice its current size. Does not perform any
    checking. If the hash tables of {D} are valid, inserts the node
    into them.  Cost: {O(1)} time and space, apart from eventual 
    expansion and rehashing costs. */

void rdag_check_invariants(rdag_t *D, rdag_node_t lo, rdag_node_t hi);
  /* Checks representation invariants for all nodes {s} in {[lo .. hi-1]},
    except node uniqueness. */

uint32_t rdag_node_hash(rdag_t *D, rdag_node_data_t *dt);
  /* Computes a hash value in {[0..D->hash_size-1]} for a node whose 
    data is {*dt}. The hash table size {D->hash_size} must be valid. */

uint32_t rdag_node_hash_packed(rdag_t *D, rdag_node_packed_t *pk);
  /* Computes a hash value in {[0..D->hash_size-1]} for a node whose
    data record is {*dt}. The hash table size {D->hash_size} must be valid. */
    
void rdag_rehash(rdag_t *D);
  /* Rebuild the hash tables of {D}, from the node data {D->pk}. */

#endif
