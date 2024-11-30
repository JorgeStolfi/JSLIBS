/* See gauss_elim_triangularize.h */
/* Last edited on 2024-11-25 03:25:12 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>

#include <gauss_elim_triangularize.h>

void gauss_elim_triangularize
  ( uint32_t m,
    uint32_t n,
    double M[],
    bool_t total,
    double tiny
  )
  { int32_t i = 0;
    int32_t j = 0;
    double *Mij = &(M[0]);
    while ((i < m) && (j < n))
      { double *Mkj;
        /* Elements {M[r][s]} with {r >= i} and {s < j} are all zero. */
        assert(isfinite(*Mij));
        if (i < m-1)
          { /* Pivoting: */
            /* Find row {kmax} in {i..m-1} that maximizes {|M[kmax,j]|} */
            int32_t kmax = i;
            double Mmax = (*Mij);
            Mkj = Mij + n;
            for (int32_t k = i+1; k < m; k++)
              { assert(isfinite(*Mkj));
                if (fabs(*Mkj) > fabs(Mmax)) { kmax = k; Mmax = (*Mkj); } 
                Mkj += n;
              }
            if (kmax != i)
              { /* Swap rows {i} and {kmax} negating one of them: */
                Mkj = &(M[kmax*(int32_t)n + j]);
                double *Mkr = Mkj; 
                double *Mir = Mij;
                for (int32_t r = j; r < n; r++)
                  { double t = (*Mkr); (*Mkr) = -(*Mir); (*Mir) = t;
                    Mkr++; Mir++;
                  }
              }
          }
        if ((*Mij) != 0.0)
          { /* Clear elements {M[k][j]} with {k > i} */
            Mkj = Mij + n;
            for (int32_t k = i+1; k < m; k++)
              { if ((*Mkj) != 0.0)
                  { /* Subtract from row {k} a multiple of row {i} that cancels {M[k,j]}: */
                    double s = (*Mkj)/(*Mij);
                    assert(isfinite(s));
                    (*Mkj) = 0.0;
                    double *Mkr = Mkj+1;
                    double *Mir = Mij+1;
                    for (int32_t r = j+1; r < n; r++)
                      { double old = (*Mkr);
                        (*Mkr) = old - s*(*Mir);
                        /* Feeble attempt to clean out entries that are just roundoff error: */
                        if (fabs(*Mkr) < tiny*fabs(old)) { (*Mkr) = 0.0; }
                        Mkr++; Mir++;
                      }
                  }
                Mkj += n;
              }
            /* Advance to next row: */
            i++; Mij += n;
          }
        else
          { if (! total) { /* Advance to next row anyway: */ i++; Mij += n; } }
        j++; Mij++;
      }
  }
