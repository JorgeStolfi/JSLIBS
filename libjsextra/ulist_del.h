#ifndef ulist_del_H
#define ulist_del_H

/* Item delete and and item replace operations for {ulist.h}. */
/* Last edited on 2009-03-06 20:37:08 by stolfi */

#include <stdint.h>

#include <vec.h>
#include <bool.h>

#include <ulist.h>

ulist_item_t ulist_delete_last(ulist_t *S);
  /* Deletes the last item of the list {S}, and returns it.
    The other elements remain in {S}, with their current indices.
    
    Fails if {S} is empty. 
    Cost: {K} expected, {K*ct(S)} worst-case. */    

ulist_item_t ulist_set(ulist_t *S, ulist_index_t i, ulist_item_t a);
  /* Sets item {i} of the list {S} to {a}.  Returns the original
    item with index {i}.
    
    Fails if the index {i} is not in the range {0..ct(S)-1},
    or if {S} contains some other item equivalent to {a}.
    Cost: {K} expected, {K*ct(S)} worst-case. */    

#endif
