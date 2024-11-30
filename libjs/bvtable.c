/* See {bvtable.h} */
/* Last edited on 2024-11-15 19:49:36 by stolfi */
  
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <bool.h>
#include <affirm.h>

#include <bvtable.h>

typedef unsigned char ubyte;

/* !!! Allow tables > 2^31 entries. !!! */

/* INTERNAL REPRESENTATION */

#define bvtable_nh_MAX ((uint32_t)1u << 31)
  /* Max number of slots in the hash table. Allows 
    indices to be any non-negative signed 32-bit integers.
    Also allows {UINT32_MAX} to be used as a null value. */

typedef struct bvtable_t
  { 
    size_t sz;       /* Size of each entry in bytes. */
    uint32_t ne;     /* Count of defined entries. */
    ubyte *pe;       /* Pointer to the sequential list of entries. */  
    uint32_t ne_max; /* Count of available entry slots in {*pe}. */
    
    uint32_t nh;     /* Number of slots in the hash table. */
    uint32_t *ih;    /* Hash table of indices into {*pe}. */
  } bvtable_t;
  /* At any moment the table holds a list of {ne} distinct items,
    each {sz} bytes long, stored in consecutive positions starting 
    at address {pe}.  The latter points to an area that has space
    for at least {ne_max} such entries (i.e. {ne_max*sz} bytes total).
    
    The lookup functions use an auxiliary table {hh[0..nh-1]} of integers.
    Each element of {ih} is either {UINT32_MAX} or an integer in {0..ne-1}.
    If an entry {X} is stored in the table, it is stored at address
    {ad(k) = pe + sz*ie(k)} where {ix(k) = ih[(a + b*k)%nh]}, 
    {a} and {b} are hash values computed from the entry {X},
    and {k} is some natural number.  The entry is not in the table
    iff {ix(k)} is {UINT32_MAX} for some {k}, and the indices
    {ix(0..k-1)} are all different from {UINT32_MAX}, and 
    the entries at {ad(0..k-1)} are all different from {X}. */
    

/* INTERNAL PROTOTYPES */

void bvtable_alloc_entries(bvtable_t *tb, uint32_t ne_max);
  /* Allocates the entry storage area {tb->pe} for maximum of {ne_max} elements, 
    and saves {ne_max} in {tb->ne_max}. */
  
void bvtable_alloc_hash(bvtable_t *tb);
  /* Chooses an appropriate size {tb->nh} for the hash vector {tb->ih},
    allocates it, and fills it with {UINT32_MAX}. 
    Assumes that {tb->ih} is NULL. */
    
uint32_t bvtable_get_hash_slot_index
  ( bvtable_t *tb, 
    ubyte *X,
    bvtable_hash_proc_t *hash,
    bvtable_cmp_proc_t *cmp
  );
  /* Fids the index {j} of the hash table slot {tb->ih[j]}
    where the index of {*X} in the list is or should be.  */
    
void bvtable_expand
  ( bvtable_t *tb, 
    bvtable_hash_proc_t *hash,
    bvtable_cmp_proc_t *cmp 
  );
  /* Expand the allocated area {tb->pe} of {tb} to hold up about twice
    as many elements as currently in it. Expands the hash table {tb->ih}
    accordingly. Re-inserts all the old elements in the table. */

void bvtable_hash_params(uint64_t h, uint32_t nh, uint32_t *aP, uint32_t *bP);
  /* Compute a starting hash index {a} and an odd hash increment {b}, both
    in the range {0..nh-1}, from a 64-bit hash value {h}.
    Expects {nh} to be a power of 2, greater than 1. */

#define bvtable_MIN_NEW 256
  /* At each expansion, double the current allocation of space ,
    for entries, but reserve at least this many entries. */

/* EXPORTED IMPLEMENTATIONS */

bvtable_t *bvtable_new(size_t sz, uint32_t ne_guess)
  { 
    bvtable_t *tb = notnull(malloc(sizeof(bvtable_t)), "no mem for {tb}");
    demand(sz > 0, "entry size must be positive");
    
    /* Set the scalar fields and allocate the {*pe} area: */
    tb->sz = sz;
    tb->ne = 0;
    tb->pe = NULL;
    bvtable_alloc_entries(tb, ne_guess);
    
    /* Allocate and initialize the hash table {tb->ih}: */
    tb->ih = NULL;
    bvtable_alloc_hash(tb);
    
    return tb;
  }
  
uint32_t bvtable_item_count(bvtable_t *tb)
  { 
    return tb->ne;
  }

uint32_t bvtable_alloc_size(bvtable_t *tb)
  { 
    return tb->ne_max;
  }

uint32_t bvtable_hash_size(bvtable_t *tb)
  { 
    return tb->nh;
  }

void *bvtable_item_address(bvtable_t *tb, uint32_t ie)
  { 
    demand(ie < tb->ne, "invalid element index");
    ubyte *pi = tb->pe + tb->sz*ie;
    return pi;
  }

uint32_t bvtable_get_index
  ( bvtable_t *tb, 
    void *X,
    bvtable_hash_proc_t *hash,
    bvtable_cmp_proc_t *cmp
  )
  { 
    /* Get the index {j} of the hash table slot where the index of {*X} should be: */
    uint32_t jh = bvtable_get_hash_slot_index(tb, (ubyte *)X, hash, cmp);
    return tb->ih[jh];
  }
   
uint32_t bvtable_add
  ( bvtable_t *tb, 
    void *X,
    bvtable_hash_proc_t *hash,
    bvtable_cmp_proc_t *cmp
  )
  { uint32_t jh = bvtable_get_hash_slot_index(tb, (ubyte *)X, hash, cmp);
    uint32_t ie = tb->ih[jh];
    if (ie != UINT32_MAX) 
      { return ie; }
    else
      { if (tb->ne >= tb->ne_max) 
          { /* Storage exausted, get some more: */
            bvtable_expand(tb, hash, cmp);
            /* Since the hash table size changed, must recompute the hash slot index: */
            jh = bvtable_get_hash_slot_index(tb, (ubyte *)X, hash, cmp);
          }
        assert(tb->ne < tb->ne_max);
        ie = tb->ne;
        ubyte *pi = tb->pe + tb->sz*ie;
        memcpy(pi, X, tb->sz);
        tb->ne++;
        tb->ih[jh] = ie;
        return ie;
      }
  }

void bvtable_close(bvtable_t *tb, uint32_t *neP, void **peP)
  { 
    tb->pe = realloc(tb->pe, tb->ne*tb->sz);
    if (tb->ne > 0) { affirm(tb->pe != NULL, "no mem for realloc"); }
    (*neP) = tb->ne;
    (*peP) = (void *)tb->pe;
    free(tb->ih);
    free(tb);
  }
   
/* INTERNAL IMPLEMENTATIONS */

#define bvtable_MIN_HASH_SIZE 1024
  /* Minimum size of an hash table.  Must be a power of 2. */
    
void bvtable_alloc_entries(bvtable_t *tb, uint32_t ne_max)
  { 
    demand(ne_max <= bvtable_ne_MAX, "attempt to alocate too many entries");
    bool_t verbose = FALSE;
    assert(tb->pe == NULL);
    if (ne_max == 0)
      { tb->ne_max = 0; }
    else
      { tb->ne_max = ne_max;
        tb->pe = notnull(malloc(ne_max*tb->sz), "no mem for {pe}");
      }
    if (verbose) { fprintf(stderr, "allocated table for %u entries\n", tb->ne_max); }
  }    
    
void bvtable_alloc_hash(bvtable_t *tb)
  { 
    bool_t verbose = FALSE;
    assert(tb->ih == NULL);
    /* Choose a hash size {2^k}, at least twice the max num of entries: */
    uint64_t nh_min = 2 * (uint64_t)tb->ne_max;
    uint64_t nh_tmp = (uint64_t)bvtable_MIN_HASH_SIZE;
    while (nh_tmp < nh_min) { nh_tmp = 2 * nh_tmp; }
    demand(nh_tmp <= (uint64_t)bvtable_nh_MAX, "attempt to allocate too many hash slots");
    tb->nh = (uint32_t)nh_tmp;
    tb->ih = notnull(malloc(tb->nh*sizeof(int32_t)), "no mem for {ih}");
    uint32_t jh;
    for (jh = 0; jh < tb->nh; jh++) { tb->ih[jh] = UINT32_MAX; }
    if (verbose) { fprintf(stderr, "allocated hash table with %u slots\n", tb->nh); }
  }    

uint32_t bvtable_get_hash_slot_index
  ( bvtable_t *tb, 
    ubyte *X,
    bvtable_hash_proc_t *hash,
    bvtable_cmp_proc_t *cmp
  )
  {
    bool_t verbose = FALSE;
    
    size_t sz = tb->sz;
    uint32_t ne = tb->ne;
    uint32_t nh = tb->nh;
    assert(nh >= 2u);
    assert((nh & (uint32_t)(nh - 1u)) == 0);
    assert(nh > ne);
    assert(tb->ih != NULL);
    
    /* Linear hash probing: */
    uint64_t h = hash(X, sz);
    uint32_t a, b; /* Hash probe start and step. */
    bvtable_hash_params(h, nh, &a, &b);
    assert(a < nh);         /* Start {a} must be in the hash table. */
    assert((b & 1u) == 1u); /* Increment {b} must be odd. */
    uint32_t jh = a;
    uint32_t npr = 0; /* Number of probes. */
    while (TRUE)
      { npr++;
        uint32_t ie = tb->ih[jh];
        if (ie == UINT32_MAX) { /* Found an empty slot: */ break; }
        /* Check whether the item at {ie} is {X}: */
        assert(ie < tb->ne);
        ubyte *pi = tb->pe + tb->sz*ie;
        if (cmp(X, pi, sz) == 0) { /* Found it: */ break; }
        /* Get next slot: */
        jh = (jh + b) % nh;
        assert(jh != a); /* There must be an empty slot. */
      }
    if (verbose && (npr > 3)) { fprintf(stderr, "a = %u b = %u probes = %u\n", a, b, npr); }
    return jh;
  }
    
void bvtable_expand
  ( bvtable_t *tb, 
    bvtable_hash_proc_t *hash,
    bvtable_cmp_proc_t *cmp 
  )
  { 
    bool_t verbose = TRUE;

    demand(tb->ne_max < bvtable_ne_MAX, "too many items, cannot expand further");
    
    /* Compute the new table size {ne_new} -- about double, within limits: */
    uint32_t ne_new; /* New max number of entries in table. */
    int64_t ne_tmp = 2*(uint64_t)tb->ne_max;
    /* Not worth allocating too few items: */
    if (ne_tmp < bvtable_MIN_NEW) { ne_tmp = bvtable_MIN_NEW; }
    /* Limit allocation to {bvtable_ne_MAX}: */
    if (ne_tmp <= (uint64_t)bvtable_ne_MAX) 
      { ne_new = (uint32_t)ne_tmp; }
    else
      { ne_new = bvtable_ne_MAX;
        fprintf(stderr, "could not double hash table, increasing from %u to %u", tb->ne_max, ne_new);
      }
    assert(ne_new > 0);

    /* Save the old entry list: */
    uint32_t no = tb->ne;
    ubyte *po = tb->pe;
    
    /* Allocate a new area, initially empty: */
    tb->ne = 0;
    tb->pe = NULL;
    if (verbose) { fprintf(stderr, "expanding table from %u to %u entries\n", tb->ne_max, ne_new); }
    bvtable_alloc_entries(tb, ne_new);
    
    /* Reallocate the hash table: */
    free(tb->ih);
    tb->ih = NULL;
    bvtable_alloc_hash(tb);
    if (verbose) { fprintf(stderr, "hash table now has %u hash slots\n", tb->nh); }
    
    /* Insert all entries again: */
    ubyte *pi = po;
    for (uint32_t ie = 0;  ie < no; ie++)
      { uint32_t iadd = bvtable_add(tb, (void *)pi, hash, cmp);
        if (iadd != ie) { fprintf(stderr, "** BUG: inconsistent rehash %u %u\n", ie, iadd); }
        assert(iadd == ie);
        pi += tb->sz;
      }
    assert(tb->ne == no);
    free(po);
  }
    
void bvtable_hash_params(uint64_t h, uint32_t nh, uint32_t *aP, uint32_t *bP)
  { 
    assert(nh >= 2);
    assert(nh <= UINT32_MAX);
    /* Ensure that {nh} is power of 2. */
    assert((nh & (uint32_t)(nh - 1)) == 0);
    /* This should be OK if {nh} is a power of 2 not exceeding {2^32}: */
    uint32_t a = (uint32_t)(h & (uint64_t)(nh - 1));
    uint32_t b = (uint32_t)((h >> 32) & (uint64_t)(nh - 1)) | 1u;
    (*aP) = a;
    (*bP) = b;
  }
    
