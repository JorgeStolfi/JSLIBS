/* See sym_eigen_tred2.h */
/* Last edited on 2024-12-05 18:43:46 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>

#include <sym_eigen_tred2.h>

void sym_eigen_tred2(uint32_t n, double A[], double d[], double e[], double Z[])
  { 
    auto void t2s300(int32_t i);
    auto void t2s500(int32_t i);
      /* Accumulate of transformation matrix. */
    
    if (n == 0)
      { /* Void case: */ }
    else if (n == 1)
      { /* Trivial case: */
        d[0] = A[0];
        Z[0] = 1.0;
        e[0] = 0;
      }
    else if (n == 2)
      { /* Already tridiagonal: */
        d[0] = A[0]; d[1] = A[3];
        Z[0] = Z[3] = 1.0; Z[1] = Z[2] = 0;
        e[0] = 0; e[1] = A[2];
      }
    else
      { /* Copy last column of {A} to {d}: */
        for (uint32_t i = 0; i < n; i++) 
          { d[i] = A[i*n + n-1]; }

        /* Copy half of {A} aabove diagonal to {Z}: */
        for (uint32_t i = 0; i < n; i++)
          { for (uint32_t j = i; j < n; j++) 
              { Z[i*n + j] = A[i*n + j]; }
          }

        for (int32_t i = (int32_t)n-1; i >= 1; i--) { t2s300(i); }
        
        for (int32_t i = 1; i < n; i++) { t2s500(i); }
        
        for (uint32_t i = 0; i < n; i++)
          { d[i] = Z[i*n + n-1];
            Z[i*n + n-1] = 0.0;
          }

        Z[(n-1)*n + n-1] = 1.0;
        e[0] = 0.0;
      }
    return;

    void t2s300(int32_t i)
      {
        int32_t r = i - 1;
        /* .......... scale row (algol tol then not needed) .......... */
        double scale = 0.0;
        for (int32_t k = 0; k <= r; k++) { scale = scale + fabs(d[k]); }

        double h = 0.0;
        if (scale == 0)
          { e[i] = d[r];
            for (int32_t j = 0; j <= r; j++)
              { d[j] = Z[j*(int32_t)n + r];
                Z[j*(int32_t)n + i] = 0.0;
                Z[i*(int32_t)n + j] = 0.0;
              }
          }
        else
          { for (int32_t k = 0; k <= r; k++)
              { d[k] = d[k] / scale;
                h = h + d[k] * d[k];
              }
            double dr = d[r];
            double g = -copysign(sqrt(h),dr);
            e[i] = scale * g;
            h = h - dr * g;
            d[r] = dr - g;
            
            /* .......... form {A}*u .......... */
            for (int32_t j = 0; j <= r; j++) { e[j] = 0.0; }
            for (int32_t j = 0; j <= r; j++)
              { double dj = d[j];
                Z[i*(int32_t)n + j] = dj;
                double gg = e[j] + Z[j*(int32_t)n + j] * dj;
                if (r > j)
                  { for (int32_t k = j+1; k <= r; k++)
                      { gg = gg + Z[j*(int32_t)n + k] * d[k];
                        e[k] = e[k] + Z[j*(int32_t)n + k] * dj;
                      } 
                  }
                e[j] = gg;
              }
              
            /* .......... form p,q .......... */
            double sum = 0.0;
            for (int32_t j = 0; j <= r; j++)
              { e[j] = e[j] / h;
                sum = sum + e[j] * d[j];
              }
            double hh = sum / (h + h);
            for (int32_t j = 0; j <= r; j++) { e[j] = e[j] - hh * d[j]; }

            /* .......... form reduced {A} .......... */
            for (int32_t j = 0; j <= r; j++) /* 280 */
              { double dj = d[j];
                double ej = e[j];
                for (int32_t k = j; k <= r; k++)
                  { Z[j*(int32_t)n + k] = Z[j*(int32_t)n + k] - dj * e[k] - ej * d[k]; }
                d[j] = Z[j*(int32_t)n + r];
                Z[j*(int32_t)n + i] = 0.0;
              }
          }
        d[i] = h;
      }
   
     void t2s500(int32_t i)
       { 
         assert(i > 0);
         uint32_t r = (uint32_t)(i - 1);
         Z[r*n + n-1] = Z[r*n + r];
         Z[r*n + r] = 1.0;
         double di = d[i];
         if (di != 0)
           { for (int32_t k = 0; k <= r; k++) { d[k] = Z[i*(int32_t)n + k] / di; }
             for (int32_t j = 0; j <= r; j++)
               { double g = 0.0;
                 for (int32_t k = 0; k <= r; k++)
                   { g = g + Z[i*(int32_t)n + k] * Z[j*(int32_t)n + k]; }
                 for (int32_t k = 0; k <= r; k++)
                   { Z[j*(int32_t)n + k] = Z[j*(int32_t)n + k] - g * d[k]; }
               }
            }
         for (int32_t k = 0; k <= r; k++) { Z[i*(int32_t)n + k] = 0.0; }
      }
  }
