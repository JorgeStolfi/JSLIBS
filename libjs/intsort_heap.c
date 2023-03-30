/* See intsort.h */
/* Last edited on 2023-03-18 11:23:54 by stolfi */

#include <intsort.h>

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <intheap.h>
#include <affirm.h>

void isrt_heapsort(int32_t *h, int32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { /* Arrange elements into a sorted heap, in reverse order: */
    int32_t m = 0;
    while (m < n) { ihp_heap_insert(h, &m, h[m], cmp, -sgn); }
    /* Pop them out in reverse order and store them at the end: */
    while (m > 0) { int32_t v = ihp_heap_pop(h, &m, cmp, -sgn); h[m] = v; }
  }
