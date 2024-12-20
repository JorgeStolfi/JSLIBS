/* See {btc_bubble_parms_copy.h} */
/* Last edited on 2024-12-05 10:23:09 by stolfi */

#include <btc_bubble_t.h>
#include <btc_bubble_parms_copy.h>

void btc_bubble_parms_copy(int nb, btc_bubble_t bp_a[], btc_bubble_t bp_b[])
  { 
    int ib;
    for (ib = 0; ib < nb; ib++) { bp_b[ib] = bp_a[ib]; }
  }

