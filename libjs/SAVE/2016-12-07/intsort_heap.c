/* See intsort.h */
/* Last edited on 2005-06-05 13:33:18 by stolfi */

#include <intsort.h>

#include <stdlib.h>
#include <intheap.h>
#include <affirm.h>

void isrt_heapsort(int *h, int n, int cmp(int x, int y), int sgn)
  { /* Arrange elements into a sorted heap, in reverse order: */
    int m = 0;
    while (m < n) { ihp_heap_insert(h, &m, h[m], cmp, -sgn); }
    /* Pop them out in reverse order and store them at the end: */
    while (m > 0) { int v = ihp_heap_pop(h, &m, cmp, -sgn); h[m] = v; }
  }
