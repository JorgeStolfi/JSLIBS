/* See set32.h */
/* Last edited on 2013-10-25 01:21:15 by stolfilocal */

#include <set32.h>
#include <stdint.h>

int set32_count(set32_t A)
  { int n = 0;
    while (A != 0) { if (A & 1) { n++; }  A >>= 1; }
    return n;
  }
    
set32_elem_t set32_elem_from_index(set32_index_t j, set32_t A)
  { set32_elem_t el = 0;
    while (A != 0) 
      { if ((A & 1) == 1) { if (j == 0) { return el; } j--; } 
        A >>= 1; el++;
      }
    return set32_NO_ELEM;
  }

set32_index_t set32_index_from_elem(set32_elem_t el, set32_t A)
  { set32_index_t j = 0;
    while (el > 0) { if ((A & 1) == 1) { j++; } A >>= 1; el--; }
    return (set32_index_t)((A & 1) == 1 ? j : set32_NO_INDEX);
  }
