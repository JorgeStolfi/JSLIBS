/* See sym_eigen_tred1.h */
/* Last edited on 2024-12-05 18:43:30 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>

#include <sym_eigen_tred1.h>

void sym_eigen_tred1(uint32_t n, double A[], double d[], double e[], double e2[])
  {
    
    auto void t1s300(int32_t i);
    auto void t1s280(int32_t r, double g, double h);
    
    if (n == 0)
      { /* Void case: */ }
    else if (n == 1)
      { /* Trivial case: */
        d[0] = A[0];
        e2[0] = 0; e[0] = 0;
      }
    else if (n == 2)
      { /* Already tridiagonal: */
        d[0] = A[0]; d[1] = A[3];
        /* Order matters because {e2} may be the same as {e}: */
        e2[0] = 0; e[0] = 0; 
        e2[1] = A[2]*A[2]; e[1] = A[2];
      }
    else
      { /* Copy last column of {A} to {d}: */
        for (int32_t i = 0; i < n; i++) { d[i] = A[i*(int32_t)n + (int32_t)n-1]; }
        /* Copy diagonal of {A} to last col of {A}: */
        for (int32_t i = 0; i < n; i++) { A[i*(int32_t)n + (int32_t)n-1] = A[i*(int32_t)n + i]; }
        /* Do the work: */
        for (int32_t i = (int32_t)n-1; i >= 0; i--) { t1s300(i); }
      }
    return;
        
    void t1s300(int32_t i)
      { 
        int32_t r = i - 1;

        if (i == 0)  { e[i] = 0.0; e2[i] = 0.0; return; } 

        /* .......... scale column (algol tol then not needed) .......... */
        double scale = 0.0;
        for (int32_t k = 0; k <= r; k++) { scale = scale + fabs(d[k]); }

        if (scale == 0.0) 
          { /* Copy {A[0..r,r] -> d[0..r]}, {A[0..r,i] -> A[0..r,r]}, clear {A[0..r,i]}: */                for (int32_t j = 0; j <= r; j++) /* 125 */
              { d[j] = A[j*(int32_t)n + r];
                A[j*(int32_t)n + r] = A[j*(int32_t)n + i];
                A[j*(int32_t)n + i] = 0.0;
              }
            e[i] = 0.0; e2[i] = 0.0;
            return;
          }

        double h = 0.0; 
        for (int32_t k = 0; k <= r; k++) /* 150 */
          { d[k] = d[k] / scale;
            h = h + d[k] * d[k];
          }

        e2[i] = scale * scale * h;
        double f = d[r];
        double g = - copysign(sqrt(h),f);
        e[i] = scale * g;
        h = h - f * g;
        d[r] = f - g;
        if (r > 0) { t1s280(r, g, h); }
        
        for (int32_t j = 0; j <= r; j++) /* 290 */
          { double dj = d[j];
            d[j] = A[j*(int32_t)n + r];
            A[j*(int32_t)n + r] = A[j*(int32_t)n + i];
            A[j*(int32_t)n + i] = dj * scale;
          }
        return;
      }

    void t1s280(int32_t r, double g, double h)
      {
        /* .......... form {A}*u .......... */
        for (int32_t j = 0; j <= r; j++) { e[j] = 0.0; }

        for (int32_t j = 0; j <= r; j++) /* 240 */
          { double dj = d[j];
            double g = e[j] + A[j*(int32_t)n + j] * dj;
            if (r > j)
              { for (int32_t k = j+1; k <= r; k++)
                  { g = g + A[j*(int32_t)n + k] * d[k];
                    e[k] = e[k] + A[j*(int32_t)n + k] * dj;
                  }
              }
            e[j] = g;
          }
          
        /* .......... form p .......... */
        double sum = 0.0;
        for (int32_t j = 0; j <= r; j++) 
          { e[j] = e[j] / h;
            sum = sum + e[j] * d[j];
          }

        double hh = sum / (h + h);
        
        /* .......... form q .......... */
        for (int32_t j = 0; j <= r; j++) { e[j] = e[j] - hh * d[j]; }
        
        /* .......... form reduced {A} .......... */
        for (int32_t j = 0; j <= r; j++) 
          { double dj = d[j], ej = e[j];
            for (int32_t k = j; k <= r; k++)
              { A[j*(int32_t)n + k] = A[j*(int32_t)n + k] - dj*e[k] - ej*d[k];  }
          }
      }
  }
