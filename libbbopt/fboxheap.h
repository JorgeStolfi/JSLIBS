/* fboxheap.h -- Heap of function boxes sorted by some criterion */
/* Last edited on 2023-02-20 06:38:37 by stolfi */

#ifndef fboxheap_h
#define fboxheap_h

#define _GNU_SOURCE
#include <stdint.h>

#include <fbox.h>

typedef struct FBoxHeap /* An f-box heap descriptor */ 
  { int32_t n;         /* Number of items in heap */
    int32_t sz;        /* Current capacity of heap. */
    FBox **b;      /* The heaped boxes are {b[0..n-1]}. */
    FBoxCmp cmp;   /* The box comparison criterion. */ 
  } FBoxHeap;

FBoxHeap *fbox_heap_new(int32_t sz, FBoxCmp cmp);
  /* Allocates a new heap with space for `sz' boxes. */

void fbox_heap_insert(FBoxHeap *h, FBox *b);
  /* Inserts box `b' in the heap `h', sorting it into its proper place. */

FBox *fbox_heap_pop(FBoxHeap *h);
  /* Removes from `h' a box that is minimal under `h->cmp', and reeturns it. */

#endif
