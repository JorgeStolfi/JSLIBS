/* See intheap.h */
/* Last edited on 2024-11-16 12:18:43 by stolfi */

#include <stdint.h>
#include <stdlib.h>

#include <bool.h>
#include <affirm.h>

#include <intheap.h>

void ihp_heap_insert(int32_t *h, uint32_t *n, int32_t v, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { /* Create a vacancy for {v} at bottom of heap: */
    uint32_t i = (*n), j;
    (*n)++;
    /* Bubble {v} rootwards to its proper place {h[i]}: */
    while ((i > 0) && (sgn*cmp(v,h[j=(i-1)/2]) < 0)) { h[i] = h[j]; i = j; }
    h[i] = v;
  }
  
int32_t ihp_heap_pop(int32_t *h, uint32_t *n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { 
    if ((*n) <= 0)
      { affirm(FALSE, "prog error: popping an empty heap"); return 0; }
    else
      { uint32_t m = (*n); 
        /* Save current root: */
        int32_t w = h[0];
        /* Declare the root as a vacant slot: */
        uint32_t i = 0;  /* {h[i]} is a vacant slot. */
        /* Promote children into vacancy {h[i]} until it reaches the fringe: */
        uint32_t ja = 1; /* {h[ja]} is the first child of {h[i]}. */
        while (ja < m)
          { /* Find largest child {h[j]} of {h[i]}: */
            uint32_t jb = ja + 1; /*  {h[jb]} is the second child of {h[i]}. */
            uint32_t j = ((jb < m) && (sgn*cmp(h[ja],h[jb]) > 0) ? jb : ja);
            /* Promote largest child into hole: */
            h[i] = h[j]; i = j; ja = 2*i + 1;
          }
        /* One less element in heap: */
        m--;
        if (i < m)
          { /* The vacancy {h[i]} did not end up at {h[m]}, so fill it with {h[m]}: */        
            int32_t v = h[m]; 
            /* Bubble it up to the proper place: */
            while (i > 0)
              { uint32_t j = (uint32_t)((i-1)/2);
                if (sgn*cmp(v,h[j]) >= 0) { break; }
                h[i] = h[j]; i = j;
              }
            h[i] = v;
          }
        *n = m;
        /* Return the removed element: */
        return w;
      }
  }
