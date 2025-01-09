/* See ulist.h */
/* Last edited on 2020-10-01 20:49:55 by jstolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <ref.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>

#include <ulist.h>
#include <ulist_del.h>
#include <ulist_impl.h>
  
ulist_item_t ulist_delete_last(ulist_t *S)
  { 
    #if (ulist_DEBUG)
    S->ct_del++;
    #endif /* (ulist_DEBUG) */
    /* Get the list count {ct}: */
    ulist_count_t ct = S->ct;
    /* We cannot delete from an empty list: */
    demand(ct > 0, "list is empty");
    /* Get the item {a} to be deleted: */
    ulist_item_t a = S->it.e[ct-1]; 
    /* Locate the index {i} and hash table slot {k} of {~a}: */
    ulist_index_t i;
    ulist_hash_val_t k;
    ulist_locate(S, a, &i, &k);
    assert(i == ct-1);
    /* Grab the item {b} to be deleted: */
    /* bool_t debug = FALSE; */
    /* bool_t debug = ((i == 749) && (k == 0)); */
    bool_t debug = ((S->ix.ne == 3711) && (a == 5833));
    if (debug) { (void)ulist_verify(S, TRUE); }
    /* bool_t debug = ((sz == 1000) && (ct == 751)); */
    /* Remove {a} from {S}: */
    S->ct--;
    /* Now the items in occupied entries that follow {ix[k]} may need reinserting: */
    if (debug) { fprintf(stderr, "-- start actual deletion --\n"); }
    if (debug) { fprintf(stderr, "  vacated ix[%u]\n", k); }
    if (debug) { (void)ulist_verify(S, FALSE); }
    /* Get the hash table size {nh}: */
    ulist_hash_size_t nh = S->ix.ne;
    if (debug) { fprintf(stderr, "  nh = %u\n", nh); }
    uint32_t r = (k + 1) % nh; /* Scans entries after {S[k]}. */
    uint32_t d = 0; /* Always equal to {imod(r - (k+1), nh)}. */ 
    while(TRUE)
      { /* Now entries {ix[k+1..r]\ix[r]}} are OK and occupied. */
        if (debug) { fprintf(stderr, "  k = %u  r = %u  d = %u\n", k, r, d); }
        ulist_index_t j = S->ix.e[r];
        if ((j >= ct) || (S->hx.e[j] != r))
          { if (debug) { fprintf(stderr, "  found vacant ix[%u]\n", r); }
            /* If our logic is correct, there are no more inconsistencies. */
            break;
          }
        /* Slot {r} is occupied by item {it[j]}. */
        if (debug) { fprintf(stderr, "  ix[%u] occupied by %u", r, j); }
        /* At this point the slot {ix[r]} is distinct from {ix[k]}. */
        /* Check whether {it[j]} needs to be reinserted: */
        ulist_hash_val_t hj = S->hash(S->it.e[j], nh);
        demand (hj < nh, "bad hash function");
        uint32_t dhj = (uint32_t)imod(hj - (k+1), nh);
        if (debug) { fprintf(stderr, "  it[%u] hashes to %u  dhj = %u  d = %u", j, hj, dhj, d); }
        if (dhj > d) 
          { /* Item {it[j]} is not OK, must remove and reinsert it. */
            if (debug) { fprintf(stderr, "  presumed BAD\n"); }
            /* Store the index {j} into {ix[k]}: */
            S->ix.e[k] = j; S->hx.e[j] = k; 
            if (debug) { fprintf(stderr, "  filled ix[%u] with ix[%u] = %u\n", k, r, j); }
            /* Now the hot hole is {ix[r]}: */
            if (debug) { fprintf(stderr, "  ix[%u] is now vacant\n", r); }
            k = r; d = 0; 
            if (debug) { (void)ulist_verify(S, FALSE); }
          }
        else
          { if (debug) { fprintf(stderr, "  presumed OK\n"); } }
        /* Go on to the next slot of {ix[]}: */
        d++; r++; if (r == nh) { r = 0; }
      }
    #if (ulist_PARANOIA)
    (void)ulist_verify(S, TRUE);
    #endif /* (ulist_PARANOIA) */
    if (debug) { fprintf(stderr, "-- end actual deletion --\n"); }
    if (debug) { (void)ulist_verify(S, TRUE); }
    /* !!! implement auto-resize if less than 1/2 full. !!! */
    return a;
  }

ulist_item_t ulist_set(ulist_t *S, ulist_index_t i, ulist_item_t a)
  { /* Get the list count {ct}: */
    ulist_count_t ct = S->ct;
    /* Check index validity: */
    demand((i >= 0) && (i < ct), "invalid item index");
    /* Swap item number {i} with last item: */
    ulist_swap(S, i, ct-1);
    /* Delete the item number {i}, now in last position, and save it in {b}: */
    ulist_item_t b = ulist_delete_last(S); 
    /* Append the new item: */
    ulist_insert_last(S, a);
    /* Swap it back to position {i}: */
    ulist_swap(S, i, ct-1);
    /* We are done: */
    return b;
  }
