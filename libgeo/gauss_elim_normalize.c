/* See gauss_elim.h */
/* Last edited on 2024-11-25 01:21:33 by stolfi */

#include <stdint.h>
#include <assert.h>

#include <gauss_elim_normalize.h>

void gauss_elim_normalize(uint32_t m, uint32_t n, double M[])
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m) && (j < n))
      { if ((*Mij) != 0.0)
          { /* Scale row {i} by {1/M[i,j]}: */
            double s = (*Mij);
            (*Mij) = 1.0;
            double *Mir = Mij+1;
            for (int32_t r = j+1; r < n; r++) { (*Mir) /= s;  Mir++; }
            i++; Mij += n;
          }
        j++; Mij++;
      }
  }

