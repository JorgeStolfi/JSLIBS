#ifndef ulist_impl_H
#define ulist_impl_H

/* Internal definitions for an implementation of {ulist.h}. */
/* Last edited on 2024-11-15 19:16:46 by stolfi */

#include <stdint.h>

#include <vec.h>
#include <ulist.h>

#define ulist_DEBUG 1
  /* Define this to 1 to cause the printout of debugging info
    and performance statistics. */

#define ulist_PARANOIA 1
  /* Define this to 1 to call {ulist_verify} after main operations. */

vec_typedef(ulist_hash_val_vec_t, ulist_hash_val_vec, ulist_hash_val_t);
vec_typedef(ulist_index_vec_t, ulist_index_vec, ulist_index_t);

struct ulist_t {
  /* The item list: */
  ulist_item_vec_t it;        /* The item storage slots. */
  ulist_count_t ct;           /* The current count of items in the {ulist_t}. */
  /* The hashing tables: */
  ulist_hash_val_vec_t hx;    /* Maps item indices to (displaced) hash values. */
  ulist_index_vec_t ix;       /* Maps (displaced) hash values to item indices. */
  /* The methods: */
  ulist_eq_func_t *eq;        /* The item equivalence predicate. */
  ulist_hash_func_t *hash;    /* The item hash function. */
  #if (ulist_DEBUG)
  /* Debugging data: */
  uint64_t ct_add;            /* Count of {ulist_insert_last} operations. */
  uint64_t ct_del;            /* Count of {ulist_delete_last} operations. */
  uint64_t ct_ind;            /* Count of {ulist_index_of} operations. */
  uint64_t ct_rsz;            /* Count of {ulist_resize} operations. */
  uint64_t ct_prb;            /* Count of probes in hash table search. */
  #endif /* (ulist_DEBUG) */
};

/* INTERNAL NOTATION
  
  In the following comments, `{~a}' means `any {ulist_item_t} value {b} such
  that {eq(a,b)} is TRUE'. 
  
  We also write 
  
    {ct} instead of {ct(S)}. 
    
    {sz} instead of {sz(S)} (that is, {S->it.ne}).
    
    {nh} instead of {nh(S)}, the size of the hash table. 
    
    {it[i]} for the variable {S->it.e[i]}, for any {i} in {0..ct-1}.
    
    {hx[i]} for the variable {S->hx.e[i]}, for any {i} in {0..ct-1}.
    
    {ix[k]} for the variable {S->ix[imod(k,nh)]}, for any integer {k}.
  
  Note that {ix} is conceptually a circular array with period {nh},
  while {it} and {hx} are ordinary (non-circular) arrays.
  
  Also, by {i..j(%n)} we denote all integers {r} in the range {i..i+d}
  where {d = imod(j-i,n)}; namely, the {d+1} consecutive integers that
  start with {i} and continue *forward* until the first integer
  congruent to {j} modulo {n} --- even if {i < j}. In particular, if
  {j==i} (or {j == i+k*n} for any {k}) then {d == 0} and {i..j(%n)} is
  just one integer {i == j}. Otherwise {i..j(%n)} has between 2 and
  {n-1} integers.
  
  The notation {ix[i..j]} is then defined as {ix[r]} for all {r} in {i..j(%nh)}.

  STRUCTURE INVARIANTS:

    I0: {ct == S->ct} and {sz == S->it.ne} and {nh = S->ix.ne}.
    
    I1: {ct <= sz <= ulist_item_count_MAX}.
    
    I2: {S->xh.ne == sz}.
    
    I3 (list): The item of {S} with index {i} is {it[i]}, for {i} in {0..ct-1}.
    
    I4 (bilink): for any {i} in {0..ct-1}, {(hx[i] < nh) && (ix[hx[i]] == i)}. 
    
    I5 (seq-search): for any {i} in {0..ct-1}, {occ(r)} is true for any
       {r} in {h..k(%nh)}, where {h = hash(it[i],nh)} and {k = hx[i]}.

    I6 (distinctness): for any {i,j} in {0..ct-1}, {eq(it[i],it[j])} implies {i==j}.
    
    I7. Either {sz == nh == 0}, or {sz < nh <= ulist_hash_size_MAX}
  
  Invariant I4 means that an entry {ix[k]} is associated to an item
  {it[i]} iff {ix[k] == i} and {hx[i] == k}. Entries of {ix[0..nh-1]}
  that pass this test are `occupied', all other are `vacant'.
  Therefore we define the predicate {occ(k)}, for {k} in {0..nh-1}, as
  {(ix[k] < ct) && (k == hx[ix[k]])}; and {vac(k)} as {! occ(k)}.

  Invariants {I4} and {I5} imply that, if {~a} is in {S}, its (unique) index
  {i} will be found in the {ix} table by searching consecutive entries
  {ix[h]}, {ix[h+1]}, {ix[h+2]}, ..., starting at its hash value {h =
  hash(a,nh)}, before finding a vacant entry.
 
 OCCUPANCY RATIO AND PERFORMANCE

  Some {ulist_t} operations call the {ulist_locate} function, that 
  finds the entry in the {ix} table for an element {a}, given the 
  value of {hash(a,nh)}.  The expected cost of the {locate} 
  function generally depends on the /relative occupancy/ of the 
  {ix} table, that is, {ro == ct/nh}.
  
  If {ro} is well below 1, the expected cost of {locate} is about {K1
  + K2*ro}. As {ro} approaches 1, the expected cost increases to
  {K2*ct}.
  
  The cost of deleting an item is about {K1 + K2*ro^2} if {ro} is
  small, and can be very high ( perhaps {K1 + K2*ct*log(ct)}? ) as {ro}
  approaches {1}. Therefore, to preserve the efficiency of those
  operations, we keep {nh} approximately equal to {2*sz}. */

void ulist_locate(ulist_t *S, ulist_item_t a, ulist_index_t *ip, ulist_hash_val_t *kp);
  /* If {~a} is in {S}, the procedure returns in {*ip} the index {i} 
    of {~a}, and in {*kp} the position {k} of that index in {ix[0..nh-1]}
    (that is, returns {i} and {k} such that {eq(it[i],a)} and {ix[k] == i}).
    
    If {~a} is not in {S}, returns in {*ip} the current count {ct}, and
    in {*kp} the position {k} of the vacant entry {ix[k]} where the index
    of {a} would have to be stored if it were appended to the list. 
    
    Fails if {ct >= nh}. */

#endif
