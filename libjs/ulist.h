#ifndef ulist_H
#define ulist_H

/* A finite list of 64-bit strings, with fast lookup by sequential hashing. */
/* Last edited on 2009-03-06 20:37:32 by stolfi */

#include <stdint.h>

#include <vec.h>
#include <bool.h>

/* 
  !!! Should generalize it to arbitrary item types (see {vec_typedef}) !!!
  
  SETS OF 64-BIT STRINGS

  This module implements a data structure {ulist_t} where one can
  store a finite set of 64-bit strings.  The module provides efficient
  procedures that locate, add, count, and enumerate the items of
  such a set. The sister interface {ulist_del.h} provides 
  additional operations that delete and replace entries of the list.

  A {ulist_t} can also store a finite set of items of any class that
  can be cast into 64-bit strings unambiguously -- such as signed
  integers, {double}s, array indices, addresses, address offsets,
  strings of up to 7 characters, etc..
  
  In the comments below we denote by {ct(S)} the current cardinality
  of the set {S}, namely the /count/ of items currently in it. 

  ITEM ORDERING
  
  An {ulist_t S} also maintains an arbitrary ordering of the items
  that are stored in it. More precisely, each item {a} stored in {S}
  has a unique /index/, denoted by {ix(S,a)}, which is an integer in
  the range {0..ct(S)-1}. This ordering can be queried and arbitrarily
  modified by the client, through the procedures {ulist_index} and
  {ulist_swap} below.
  
  Thus, a {uint64_t} actually stores a finite /list/ of distinct items
  rather than just a finite /set/ thereof; hence its name.
  (Technically, it stores an /arrangement/, rather than a
  /combination/).

  ITEM EQUIVALENCE

  The module also allows the client to specify an arbitrary
  equivalence relation over the items. This feature allows a
  {ulist_t} to be used as an efficient dictionary or
  associative array, as explained further on.

  OPERATION COSTS

  The operation costs are approximate upper bounds; the factors
  {K,K1,K2} etc. appearing in those formulas are constants that
  depend on the platform and implementation, and are usually
  different for each operation. If not said otherwise, the
  costs are both worst-case and expected (each with different
  constants).  The formulas assume that the cost of a {malloc(n)}
  call is bounded by a constant independent of {n}.
  
  SET CAPACITY

  The number of items that can be stored in a given {ulist_t S} is
  limited to a certain maximum called the the /capacity/ of {S} and
  denoted here by {sz(S)}. The capacity of {S} is defined when {S} is
  created, but it may be increased (or decreased) when needed by
  calling {ulist_resize} below.

  A {ulist_t S} occupies about {20*sz(S)} bytes of memory, i.e. 20
  bytes per item slot (whether filled or not): 8 bytes for the item
  itself, and 12 bytes for data structure overhead). */

typedef struct ulist_t ulist_t;
  /* A {ulist_t} is a dynamic finite list of {ulist_item_t} values. */

typedef uint64_t ulist_item_t;
  /* A {ulist_item_t} is an item that can be stored in an {ulist_t},
    namey a 64-bit string (actually a 64-bit {unsigned int}). */

typedef uint32_t ulist_count_t;
  /* A count of items in an {ulist_t}. */

typedef uint32_t ulist_index_t;
  /* A {ulist_index_t} value specifes the position of some item in the
    ordering of items of a set {S}. The valid values range from
    {0} (first item) to {ct(S)-1} (last item). The value {ct(S)}
    may be used to mean  `no such item'. */
    
#define ulist_item_count_MAX (2u << 30)
  /* The maximum capacity (hence maximum item count) of any {ulist_t}. */

ulist_t *ulist_new(ulist_count_t n); 
  /* Creates a new {ulist_t} with capacity {n}, initially empty.
    Cost: {K*n}. */

void ulist_free(ulist_t *S); 
  /* Releases all the internal storage used by {S}, including 
    the record {*S} itself. Cost: {K}. */

ulist_count_t ulist_count(ulist_t *S);
  /* The cardinality {ct(S)} of {S}, namely the number of
    items currently in {S}. Cost: {K}. */

bool_t ulist_is_empty(ulist_t *S);
  /*  Same as {ulist_count(S) == 0}.  Cost: {K}. */

ulist_count_t ulist_capacity(ulist_t *S);
  /* The capacity {sz(S)} of {S}, namely the maximum number of items that
    can be stored in {S} without resizing it. Cost: {K}. */

void ulist_clear(ulist_t *S);
  /* Makes {S} empty, i.e. deletes all items from {S}. Does not change
    its capacity. Cost: {K}. */

ulist_item_t ulist_item_at(ulist_t *S, ulist_index_t i);
  /* Returns the item in {S} with index {i}. 
    
    Fails if {i} is not in the range {0..ct(S)-1}.
    Cost: {K} expected, {K*ct(S)} worst-case. */
  
ulist_index_t ulist_index_of(ulist_t *S, ulist_item_t a);
  /* If {S} has the item {a}, or any item equivalent to it, returns
    the current index of that item, in the range {0..ct(S)-1};
    otherwise returns {ct(S)}.
    
    Cost: {K} expected, {K*ct(S)} worst-case. */

bool_t ulist_has(ulist_t *S, ulist_item_t a);
  /* Returns TRUE if {S} has the item {a}, or any item equivalent to
    it; returns FALSE otherwise. 
    
    Equivalent to {ulist_index_of(S,a) < ulist_count(S)}.
    Cost: {K} expected, {K*ct(S)} worst-case. */
    
ulist_index_t ulist_insert_last(ulist_t *S, ulist_item_t a);
  /* Adds the item {a} to the list {S}, as the last item of the list..
    The other elements remain in {S}, with their current indices.
    Returns the new index of {a} in {S}, namely {ct(S)-1}. 
    
    Fails if {S} is full (i.e. if. {ct(S) = sz(S)}), 
    or if {S} already contains {a} (or any item equivalent to it). 
    Cost: {K} expected, {K*sz(S)} worst-case. */
    
void ulist_swap(ulist_t *S, ulist_index_t i, ulist_index_t j);
  /* Swaps the positions of the two items of {S} with indices {i} and
    {j}. (This is a no-op if {i==j}.) 
    
    Fails if either index lies outside the range {0..ct(S)-1}.
    Cost: {K}. */

vec_typedef(ulist_item_vec_t, ulist_item_vec, ulist_item_t);
  /* Defines the type {ulist_item_vec_t}, a vector of {ulist_item_t}
    values, and associated operations. */

ulist_item_vec_t ulist_items(ulist_t *S);
  /* Returns a vector that contains all the items currently in {S},
    in arbitrary order. Cost: {K*ct(S)}. */
    
/* ITEM EQUIVALENCE

  The client may attach to a set {S} with an arbitrary /item
  equivalence predicate/ {eq(a,b)}, a function that compares two
  {uint64_t} values and returns a boolean result.
  
  The procedures in this module will then consider two items {a,b}
  /equivalent/, for the purposes of membership in {S}, if and only if
  {eq(a,b)} is TRUE. In that case, only one of those items can be in {S} at
  the same time; and that item will be found when searching {S} for
  either of them.

  Thus, for example, if the items are actually pointers to C strings,
  one can use the equivalence predicate {eq(a,b) ==
  (strcmp((char*)a,(char*)b)==0)} to treat all strings with same
  contents as the same item. Then {S} becomes effectively a 
  set of abstract character sequences.
  
  As another example, if the items are pairs {{key,val}} (or pointers
  to such pairs), one can use {eq(a,b) == (a.key == b.key)}; in that
  case, the {ulist_t} will function as a dynamic dictionary
  indexed by {key}.
  
  The default equivalence function (used when the client does not specify 
  {eq}, or specifies {eq == NULL}) is the trivial equality of {uint64_t}
  values,  namely {eq(a,b) == (a == b)}.
  
  In any case, the predicate {eq} must be an equivalence relation:
  that is, the formulas {eq(a,a)}, {eq(a,b)==eq(b,a)} and
  {(eq(a,b)&eq(b,c)) <= eq(a,c)} must be true for any items {a,b,c}
  that may be inserted in {S}.
  
  Also, the equivalence function of a set {S} should not change its
  behavior over time in any way that may cause two distinct items
  {a,b} that are currently in {S} to become equivalent. On the other
  hand, it is OK to modify or replace {eq} when {S} is empty, or to
  make {eq} more demanding (i.e. to change {eq(a,b)} from TRUE to
  FALSE for some pairs {a,b}). */

typedef bool_t ulist_eq_func_t(ulist_item_t a, ulist_item_t b);
  /* Type of a client-provided equivalence predicate.  */

void ulist_install_eq(ulist_t *S, ulist_eq_func_t *eq);
  /* Installs {eq} as the item equivalence predicate for subsequent
    operations on {S}.  
    
    Fails if {S} currently contains more than one item. Cost: {K}. */

ulist_eq_func_t *ulist_get_eq(ulist_t *S);
  /* Returns the current equivalence predicate of {S}. Cost: {K}. */

/* HASHING FUNCTION

  To speed up the main set operations, a {ulist_t} is implemented
  as a hash table with sequential search, whose hash function {hash} 
  may be provided by the client.
  
  The call {hash(a,n)} should map any {item_t} value {a} to an
  unsigned integer in the range {0..n-1}. The function should be
  chosen so as to distribute all the items that may be inserted into
  {S} as uniformly as possible over that range.
  
  The default hashing function (which is used when the client does not
  specify {hash}, or specifies {hash == NULL}) is a general purpose
  hashing function that looks at {a} as a meaningless 64-bit string;
  namely, {ulist_uint64_hash} below.
  
  The hashing function of a {ulist_t S} may be replaced only while the
  set is empty, or when it is about to be resized. Between such
  events, {hash} must be consistent --- always return the same result
  when given the same arguments {a} and {n}.
  
  Moreover, the hashing function must always be consistent with the
  item equivalence relation {eq} in use; that is, {eq(a,b)} must imply
  {hash(a,n) == hash(b,n)} for any {a,b,n}. */

typedef uint32_t ulist_hash_size_t; 
  /* Type of the size of some hash table. */

#define ulist_hash_size_MAX (1u << 31)
  /* The maximum valid value for a {ulist_hash_size_t}. */

typedef uint32_t ulist_hash_val_t; 
  /* Type of the result of some hashing function. */

#define ulist_hash_val_MAX ((unsigned)(ulist_hash_size_MAX - 1u))
  /* The maximum valid value for a {ulist_hash_val_t}. */

typedef ulist_hash_val_t ulist_hash_func_t(ulist_item_t a, ulist_hash_size_t n);
  /* Type of a client-provided item hashing function. */

void ulist_install_hash(ulist_t *S, ulist_hash_func_t *hash);
  /* Installs {hash} as the item hash function for {S}. 
    
    Fails if {S} is not empty. Cost: {K}. */

ulist_hash_func_t *ulist_get_hash(ulist_t *S);
  /* Returns the current hash function of {S}. Cost: {K}. */

ulist_hash_size_t ulist_hash_size(ulist_t *S);
  /* Returns the hash table size of {S}, namely the parameter {n} that will
    be passed to the item hash function of {S} by various operations. 
    Note that this parameter is changed by {ulist_resize} below. */

/* RESIZING A LIST */

void ulist_resize(ulist_t *S, int n, ulist_hash_func_t *hash);
  /* Changes the capacity {sz(S)} of {S} to be exactly {n}. Also
    changes the hash function to be {hash}. All the items stored in
    {S} are preserved, in their current ordering.
    
    Fails if {S} already contains more than {n} items.
    Cost: {K1*ct(S) + K2}.
    
    This function is typically called before adding an item to a list
    which has reached its capacity. It may also be called to reclaim
    unused space when the final item count of an over-allocated list
    becomes known; or when the item count becomes much smaller than
    the capacity as a result of deletions.
    
    The capacity of {S} should be expanded or reduced by a fixed
    percentage (*not* a fixed amount) every time the item count
    {ct(S)} becomes too large or too small relative to {sz(S)}. That
    way, the cost of all {ulist_resize} calls, divided by the number
    of {ulist_add} and {ulist_delete} operations, will amount to only
    a constant overhead per operation.
    
    The optimal resize policy depends on the application.  However, 
    a rule that ensures the constant amortized overhead is 
    to double {sz(S)} when {ct(S) >= sz(S)} before a {ulist_add}, 
    and to reduce {sz(S)} by half when {ct(S) < sz(S)/4} after
    a {ulist_delete}. */

/* TRIVIAL BITSTRING EQUIVALENCE AND HASHING */

ulist_hash_val_t ulist_uint64_hash(ulist_item_t a, ulist_hash_size_t n);
  /* Returns a hash value in {0..n-1}, computed from the bits of {a} alone.
    Fails if {n} is 0. Cost: {K}. */

bool_t ulist_uint64_eq(ulist_item_t a, ulist_item_t b);
  /* Returns TRUE iff {a} and {b} are the same 64-bit strings. Cost: {K}. */

/* DEBUGGING */

bool_t ulist_verify(ulist_t *S, bool_t die);
  /* Checks whether {S} satisfies the main data structure invariants.
    If it does, prints nothing and returns TRUE. If it fails any invariant,
    prints a diagnostic  message; then either aborts (if {die} is TRUE)
    or returns FALSE (if {die} is FALSE).  
    
    Does not check whether all items in {S} are pairwise
    non-equivalent. Does not check whether {eq} is an equivalence
    relation, nor whether {hash} is consistent with {eq}.
    Cost: {K*ct(S)}. */

void ulist_stats_print(ulist_t *S);
  /* Prints usage and performance statistics of {S}, such as operation
    counts, probes, etc. */

void ulist_stats_clear(ulist_t *S);
  /* Resets the statistics of {S} to zero. */

#endif
