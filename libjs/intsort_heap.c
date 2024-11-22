/* See intsort.h */
/* Last edited on 2024-11-17 15:50:11 by stolfi */

#include <stdint.h>
#include <stdlib.h>

#include <intheap.h>
#include <affirm.h>

#include <intsort.h>

void isrt_heapsort(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { if (n <= 1) { return; }
    /* Arrange elements into a sorted heap, in reverse order: */
    uint32_t m = 0;
    while (m < n) { ihp_heap_insert(h, &m, h[m], cmp, -sgn); }
    /* Pop them out in reverse order and store them at the end: */
    while (m > 0) { int32_t v = ihp_heap_pop(h, &m, cmp, -sgn); h[m] = v; }
  }
