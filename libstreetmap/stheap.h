/* stheap.h -- Heap of street map vertices sorted by cost */
/* Last edited on 2022-10-20 07:34:21 by stolfi */ 

#ifndef stheap_H
#define stheap_H

#include <stdint.h>

typedef struct st_Heap /* A heap of edge IDs */ 
  { int32_t n;     /* Number of items in heap. */
    int32_t sz;    /* Current max capacity of heap. */
    int32_t *ei;   /* The heaped edge IDs. */
  } st_Heap;

st_Heap *st_heap_new(int32_t sz);
  /* Allocates a new heap with space for {sz} edge IDs. */

void st_heap_insert(st_Heap *h, int32_t ei, float *c);
  /* Inserts edge ID {ei} in the heap {h}, sorting it by {c[ei]}. */

int32_t st_heap_pop(st_Heap *h, float *c);
  /* Removes from {h} the edge {ei} that has minimal {c[ei]}, and returns it. */

void st_heap_discard(st_Heap *h);
  /* Frees the storage allocated for {h}. */

#endif
