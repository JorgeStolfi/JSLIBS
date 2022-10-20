/* See stheap.h */
/* Last edited on 2022-10-20 05:59:22 by stolfi */

#include <stheap.h>
#include <stdint.h>
#include <stdlib.h>
#include <affirm.h>

void st_heap_resize(st_Heap *h, int32_t sz);
  /* Resizes the heap to have space for exactly {sz} boxes. */

st_Heap *st_heap_new(int32_t sz)
  { st_Heap *h;
    h = (st_Heap *)malloc(sizeof(st_Heap));
    affirm(h != NULL, "out of memory while allocating edge heap header");
    h->n = 0;
    h->sz = sz;
    h->ei = (int32_t *)malloc(sz * sizeof(int32_t));
    affirm(h->ei != NULL, "out of memory while allocating edge heap storage");
    return h;
  }
  
void st_heap_resize(st_Heap *h, int32_t sz)
  { int32_t *ei; int32_t i; int32_t n = h->n;
    ei = (int32_t *)malloc(sz * sizeof(int32_t));
    affirm(ei != NULL, "out of memory while extending edge heap");
    for (i = 0; i < n; i++) { ei[i] = h->ei[i]; }
    free(h->ei);
    h->ei = ei; h->sz = sz;
  }

void st_heap_insert(st_Heap *h, int32_t ei, float *c)
  { int32_t i, j;
    i = h->n;
    if (h->n+1 > h->sz) { st_heap_resize(h, 2*h->n+1); }
    affirm(i < h->sz, "failed to extend edge heap");
    /* Insert at bottom of heap: */
    h->ei[i] = ei; h->n++;
    /* Bubble it up to its proper place: */
    while (i > 0)
      { j = (i - 1) / 2;
        if (c[h->ei[i]] < c[h->ei[j]])
          { int32_t t = h->ei[i];  h->ei[i] = h->ei[j]; h->ei[j] = t; }
        i = j;
      } 
  }
  
int32_t st_heap_pop(st_Heap *h, float *c)
  { int32_t n = h->n; int32_t i, j, ia, ib; int32_t ei;
    affirm(n > 0, "prog error: popping an empty heap");
    /* Save current root: */
    j = 0; ei = h->ei[j];
    n--;
    if (j < n) 
      { /* Fill the hole with the last element: */
        h->ei[j] = h->ei[n];
        /* Bubble it down to the base: */
        while ((ia = 2*j + 1) < n)
          { ib = ia + 1; 
            i = ((ib >= n) || (c[h->ei[ia]] < c[h->ei[ib]]) ? ia : ib);
            if (c[h->ei[i]] < c[h->ei[j]])
              { int32_t t = h->ei[i]; h->ei[i] = h->ei[j]; h->ei[j] = t; j = i; }
            else
              { j = n; }
          }
      }
    h->n = n;
    /* if (h->n < h->sz/4) { st_heap_resize(h, 2*h->n+1); }  */
    return ei;
  }

void st_heap_discard(st_Heap *h)
  { if (h != NULL)
      { if (h->ei != NULL) { free(h->ei); }
        free(h);
      }
  }
