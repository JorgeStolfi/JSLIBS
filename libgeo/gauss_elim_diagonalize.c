/* See gauss_elim_diagonalize.h */
/* Last edited on 2024-11-25 01:20:14 by stolfi */

#include <stdint.h>
#include <assert.h>

#include <gauss_elim_diagonalize.h>

void gauss_elim_diagonalize(uint32_t m, uint32_t n, double M[])
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m) && (j < n))
      { if ((*Mij) != 0.0)
          { /* Clear elements {M[k][j]} with {k < i} */
            double *Mkj = Mij-n;
            for (int32_t k = (int32_t)i-1; k >= 0; k--)
              { /* Sub from row {k} a multiple of row {i} that cancels {M[k,j]}: */
                if ((*Mkj) != 0.0)
                  { double s = (*Mkj)/(*Mij);
                    (*Mkj) = 0.0;
                    double *Mkr = Mkj+1; 
                    double *Mir = Mij+1;
                    for (int32_t r = j+1; r < n; r++)
                      { (*Mkr) -= s*(*Mir); Mkr++; Mir++; }
                  }
                Mkj -= n;
              }
            i++; Mij += n;
          }
        j++; Mij++;
      }
  }
