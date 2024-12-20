/* See {btc_date_lookup.h} */
/* Last edited on 2024-12-05 10:23:26 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <bool.h>

#include <btc_date_lookup.h>

int btc_date_lookup(int nd, char* dt[], char* dtx)
  { /* Binary search: */
    int i = 0, j = nd-1;
    while (i <= j)
      { int m = (i + j)/2;
        int c = strcmp(dtx, dt[m]);
        if (c < 0) 
          { j = m - 1; }
        else if (c > 0) 
          { i = m + 1; }
        else
          { return m; }
      }
    fprintf(stderr, "** date %s not found\n", dtx);
    assert(FALSE);
  }
