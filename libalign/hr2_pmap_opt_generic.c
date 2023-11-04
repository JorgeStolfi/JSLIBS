/* See {hr2_pmap_opt_generic.h}. */
/* Last edited on 2023-10-20 18:46:58 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <affirm.h>

#include <hr2_pmap_opt_generic.h>

void hr2_pmap_opt_generic_encode(hr2_pmap_t *M, r3x3_t *R, double y[])
  { int32_t k = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { double Rij = R->c[i][j];
            demand(Rij >= 0, "element variation must be non-negative");
            if (Rij > 0) 
              { y[k] = (M->dir.c[i][j] - (i == j ? 1 : 0)) / Rij;  k++; }
          }
      }
    assert(k <= 8);
  }

void hr2_pmap_opt_generic_decode(double y[], r3x3_t *R, hr2_pmap_t *M)
  { int32_t k = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { double Rij = R->c[i][j];
            demand(Rij >= 0, "element variation must be non-negative");
            if (Rij > 0) 
              { M->dir.c[i][j] = y[k]*Rij + (i == j ? 1 : 0); k++; }
            else
              { M->dir.c[i][j] = (i == j ? 1 : 0); }
          }
      }
    assert(k <= 8);
    r3x3_inv(&(M->dir), &(M->inv));
  }
