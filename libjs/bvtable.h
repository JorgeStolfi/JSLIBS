/* An append_only list of distinct arbitrary data with fast lookup. */
/* Last edited on 2016-04-01 00:43:13 by stolfilocal */

#ifndef bvtable_H
#define bvtable_H

#define _GNU_SOURCE
#include <stdio.h>

#include <vec.h>
#include <bool.h>

typedef struct bvtable_t bvtable_t;
  /* A hash tabe for items of arbitrary size. */
    
bvtable_t *bvtable_new(size_t sz, uint32_t ne_guess);
  /* Creates a lookup table for entries of size {sz}.
  
    The table is preallocated with space for {ne_guess} entries. This
    parameter can be any non-negative value, even zero; but the table is
    more efficient if {ne_guess} is equal to the actual number of
    entries that will be stored, or slightly larger. */
 
uint32_t bvtable_item_count(bvtable_t *tb);
  /* Returns the current count of distinct items stored in {tb}. */
   
typedef int bvtable_cmp_proc_t(void *x, void *y, size_t sz);
  /* Type of a procedure that compares two items, each {sz} bytes long,
    that start at {*x} and {*y}, respectively. */
    
/* !!! Should use boolean {eq} instead of signed {cmp}. !!! */

typedef uint64_t bvtable_hash_proc_t(void *p, size_t sz);
  /* Type of a hashing procedure that computes a 64-bit unsigned
    hash value from the {sz} bytes starting at {*p}. */

uint32_t bvtable_get_index
  ( bvtable_t *tb, 
    void *X,
    bvtable_hash_proc_t *hash,
    bvtable_cmp_proc_t *cmp
  );
  /* Checks whether the entry {*X} is in the lookup table {tb}. The
    entry {X} is assumed to have {sz} bytes, where {sz} is the size
    passed to to {bvtable_new}. If {*X} is present,
    returns its index in the list; otherwise returns {UINT32_MAX}.
    
    The procedure uses {cmp(x,y,sz)} to compare two items in the 
    table, where {sz} is the item size as defined in the table's creation.
    The {cmp} function should return {-1}, {0}, or {+1}, and
    must implement an ordering relation with equivalent elements.  Namely,
    it must satisfy {cmp(x,x,sz) == 0}, {cmp(x,y,sz) = -cmp(y,x,sz)},
    and {cmp(x,z,sz) = +1} if {cmp(x,y,sz) == cmp(y,z,sz) == +1}.
    It may return 0 even when {x} and {y} are distinct, though.
    
    The procedure uses {hash(x,sz)} to compute a 64-bit hash value 
    from each item {x} in the table.  The {hash} procedure must be
    consistent with {cmp}, in the sense that {cmp(x,y,sz)==0}
    implies {hash(x,sz)==hash(y,sz)}. */
    
uint32_t bvtable_add
  ( bvtable_t *tb, 
    void *X,
    bvtable_hash_proc_t *hash,
    bvtable_cmp_proc_t *cmp
  );
  /* Same as {bvtable_get_index}; however, if the item {*X} is not
    present, appends it to the list, reallocating the internal storage
    as needed, without changing the order of the other items. In any
    case, returns the index of {*X} in the list. */

void bvtable_close(bvtable_t *tb, uint32_t *neP, void **peP);
  /* Terminates the addition of new elements to {tb}
    and releases all internally allocated storage
    (including the table header {*tb} itself).
    
    The number {ne} of elements present at the time is stored in {*neP},
    and the address {pe} of the area containing the elements, trimmed to
    {ne} entries, is stored in {*peP}. If {ne} is zero, {*peP} is set to
    {NULL}. */

/* LIMITS */

#define bvtable_ne_MAX ((uint32_t)1u << 30)
  /* Max number of entries in a lookup table. Allows 
    indices to be any non-negative signed 32-bit integers.
    Also allows {UINT32_MAX} to be used as a null value. */

/* UNSAFE AND SPECIALIZED PROCEDURES */

uint32_t bvtable_alloc_size(bvtable_t *tb);
  /* Returns the number of items that can be stored in 
    {tb} (including those already there) without
    causing a table expansion. */
    
uint32_t bvtable_hash_size(bvtable_t *tb);
  /* Returns the number of slots in the internal
    hash table.  Will change if and when 
    the table gets automatically 
    expanded by {bvtable_add}. */

void *bvtable_item_address(bvtable_t *tb, uint32_t ie);
  /* Returns the current address of the element at index {ie}.
    Fails if the index is not in {0..ne-1} where
    {ne} is the current element count.
  
    WARNING: this address will change if the table
    gets automatically expanded by {bvtable_add}, 
    or after the table is closed. */

#endif
