/* See intheap.h */
/* Last edited on 2023-03-18 11:27:34 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>

#include <bool.h>
#include <affirm.h>

#include <intheap.h>

void ihp_heap_insert(int32_t *h, int32_t *n, int32_t v, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { /* Create a vacancy for {v} at bottom of heap: */
    int32_t i = (*n), j;
    (*n)++;
    /* Bubble {v} rootwards to its proper place {h[i]}: */
    while ((i > 0) && (sgn*cmp(v,h[j=(i-1)/2]) < 0)) { h[i] = h[j]; i = j; }
    h[i] = v;
  }
  
int32_t ihp_heap_pop(int32_t *h, int32_t *n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { 
    if ((*n) <= 0)
      { affirm(FALSE, "prog error: popping an empty heap"); return 0; }
    else
      { int32_t m = (*n); 
        /* Save current root: */
        int32_t w = h[0];
        /* Declare the root as a vacant slot: */
        int32_t i = 0;  /* {h[i]} is a vacant slot. */
        /* Promote children into vacancy {h[i]} until it reaches the fringe: */
        int32_t ja = 1; /* {h[ja]} is the first child of {h[i]}. */
        while (ja < m)
          { /* Find largest child {h[j]} of {h[i]}: */
            int32_t jb = ja + 1; /*  {h[jb]} is the second child of {h[i]}. */
            int32_t j = ((jb < m) && (sgn*cmp(h[ja],h[jb]) > 0) ? jb : ja);
            /* Promote largest child into hole: */
            h[i] = h[j]; i = j; ja = 2*i + 1;
          }
        /* One less element in heap: */
        m--;
        if (i < m)
          { /* The vacancy {h[i]} did not end up at {h[m]}, so fill it with {h[m]}: */        
            int32_t v = h[m], j; 
            /* Bubble it up to the proper place: */
            while ((i > 0) && (sgn*cmp(v,h[j=(i-1)/2]) < 0)) { h[i] = h[j]; i = j; }
            h[i] = v;
          }
        *n = m;
        /* Return the removed element: */
        return w;
      }
  }
