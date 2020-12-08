/* See fboxheap.h */
/* Last edited on 2005-06-05 15:22:58 by stolfi */

#include <fbox.h>
#include <fboxheap.h>
#include <stdlib.h>
#include <affirm.h>

void fbox_heap_resize(FBoxHeap *h, int sz);
  /* Resizes the heap to have space for exaclty `sz' boxes. */

FBoxHeap *fbox_heap_new(int sz, FBoxCmp cmp)
  { FBoxHeap *h;
    h = (FBoxHeap *)malloc(sizeof(FBoxHeap));
    affirm(h != NULL, "out of memory while allocating f-box heap header");
    h->n = 0;
    h->sz = sz;
    h->b = (FBox **)malloc(sz * sizeof(FBox *));
    affirm(h->b != NULL, "out of memory while allocating f-box heap storage");
    h->cmp = cmp;
    return h;
  }
  
void fbox_heap_resize(FBoxHeap *h, int sz)
  { FBox **b; int i; int n = h->n;
    b = (FBox **)malloc(sz * sizeof(FBox *));
    affirm(b != NULL, "out of memory while extending box heap");
    for (i = 0; i < n; i++) { b[i] = h->b[i]; }
    free(h->b);
    h->b = b; h->sz = sz;
  }

void fbox_heap_insert(FBoxHeap *h, FBox *b)
  { int i, j;
    i = h->n;
    if (i >= h->sz) { fbox_heap_resize(h, 2*h->n+1); }
    affirm(i < h->sz, "failed to extend box heap");
    /* Insert at bottom of heap: */
    h->b[i] = b; h->n++;
    /* Bubble it up to its proper place: */
    while (i > 0)
      { j = (i - 1) / 2;
        if (h->cmp(h->b[i], h->b[j]) < 0)
          { FBox *t = h->b[i];  h->b[i] = h->b[j]; h->b[j] = t; }
        i = j;
      } 
  }
  
FBox *fbox_heap_pop(FBoxHeap *h)
  { int n = h->n; int i, j, ia, ib; FBox *b;
    affirm(n > 0, "prog error: popping an empty heap");
    /* Save current root: */
    j = 0; b = h->b[j];
    n--;
    if (j < n) 
      { /* Fill the hole with the last element: */
        h->b[j] = h->b[n];
        /* Bubble it down to the base: */
        while ((ia = 2*j + 1) < n)
          { ib = ia + 1; 
            i = ((ib >= n) || (h->cmp(h->b[ia], h->b[ib]) < 0) ? ia : ib);
            if (h->cmp(h->b[i], h->b[j]) < 0)
              { FBox *t = h->b[i]; h->b[i] = h->b[j]; h->b[j] = t; j = i; }
            else
              { j = n; }
          }
      }
    h->n = n;
    if (h->n < h->sz/4) { fbox_heap_resize(h, 2*h->n+1); } 
    return b;
  }

    
