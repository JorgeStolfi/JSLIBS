/* See intheap.h */
/* Last edited on 2014-05-15 23:53:07 by stolfilocal */

#include <intheap.h>
#include <bool.h>
#include <affirm.h>
#include <stdlib.h>

void ihp_heap_insert(int *h, int *n, int v, int cmp(int x, int y), int sgn)
  { /* Create a vacancy for {v} at bottom of heap: */
    int i = (*n), j;
    (*n)++;
    /* Bubble {v} rootwards to its proper place {h[i]}: */
    while ((i > 0) && (sgn*cmp(v,h[j=(i-1)/2]) < 0)) { h[i] = h[j]; i = j; }
    h[i] = v;
  }
  
int ihp_heap_pop(int *h, int *n, int cmp(int x, int y), int sgn)
  { 
    if ((*n) <= 0)
      { affirm(FALSE, "prog error: popping an empty heap"); return 0; }
    else
      { int m = (*n); 
        /* Save current root: */
        int w = h[0];
        /* Declare the root as a vacant slot: */
        int i = 0;  /* {h[i]} is a vacant slot. */
        /* Promote children into vacancy {h[i]} until it reaches the fringe: */
        int ja = 1; /* {h[ja]} is the first child of {h[i]}. */
        while (ja < m)
          { /* Find largest child {h[j]} of {h[i]}: */
            int jb = ja + 1; /*  {h[jb]} is the second child of {h[i]}. */
            int j = ((jb < m) && (sgn*cmp(h[ja],h[jb]) > 0) ? jb : ja);
            /* Promote largest child into hole: */
            h[i] = h[j]; i = j; ja = 2*i + 1;
          }
        /* One less element in heap: */
        m--;
        if (i < m)
          { /* The vacancy {h[i]} did not end up at {h[m]}, so fill it with {h[m]}: */        
            int v = h[m], j; 
            /* Bubble it up to the proper place: */
            while ((i > 0) && (sgn*cmp(v,h[j=(i-1)/2]) < 0)) { h[i] = h[j]; i = j; }
            h[i] = v;
          }
        *n = m;
        /* Return the removed element: */
        return w;
      }
  }
