/* See fbox.h */
/* Last edited on 2024-12-05 10:21:46 by stolfi */

#include <stdint.h>
#include <stdlib.h>

#include <ia.h>

#include <fbox.h>

FBox *fbox_make(int32_t d, int32_t depth, Interval *xr, Interval fr)
  { FBox *b = (FBox *)malloc(sizeof(FBox) + (d-1)*sizeof(Interval));
    int32_t i;
    b->d = d;
    b->depth = depth;
    b->fr = fr;
    for (i = 0; i < d; i++) { b->xr[i] = xr[i]; }
    return b;
  }
    
void fbox_discard(FBox *b)
  { free(b); }

void fbox_print(FILE *f, FBox *b)
  { int32_t i;
    fprintf(f, "(%02d)", b->depth);
    ia_print(f, b->xr[0]);
    for (i = 0; i < b->d; i++)
      { fprintf(f, " × ");
        ia_print(f, b->xr[i]); 
      }
    fprintf(f, " f = ");
    ia_print(f, b->fr); 
    fflush(f);
  }
