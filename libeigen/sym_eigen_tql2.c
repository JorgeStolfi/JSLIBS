/* See sym_eigen_tql2.h */
/* Last edited on 2024-12-05 18:44:56 by stolfi */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>

#include <sym_eigen_tql2.h>

#define Pr fprintf
#define Er stderr
  /* Shoter names. */

void sym_eigen_tql2(uint32_t n, double d[], double e[], double z[], uint32_t *nev_P)
  {
    /* The algorithm is based on the EISPACK routines "tql1.f" and
      "tql2.f" by Burton S. Garbow, Argonne National Laboratory
      (aug/1983); orignally from the ALGOL procedures {tql1} and
      {tql2}, Num. Math. 11, 293-306(1968) by Bowdler, Martin,
      Reinsch, and Wilkinson. See Handbook for Auto. Comp., vol.II -
      Linear Algebra, 227-240(1971). Re-implemented in C by Jorge
      Stolfi, Unicamp (dec/2002). */
    
    bool_t debug = FALSE;
    
    if (debug) { Pr(Er, "  > --- %s  n = %d ---\n", __FUNCTION__, n); }
    
    auto void tq2s240(void);

    auto void tq2s200(int32_t m, int32_t r, int32_t r1, double dr1);
      /* The QL transform. */

    auto void tq2s300(void);
      /* Sorts eigenvalues and eigenvectors. */

    if (n == 0) 
      { /* Nothing to do */
        (*nev_P) = 0;
      }
    else if (n == 1)
      { e[0] = d[0];
        z[0] = 1.0;
        (*nev_P) = 1; 
      }
    else
      { (*nev_P) = UINT32_MAX; /* To be redefined. */

        for (int32_t i = 1; i < n; i++) { e[i-1] = e[i]; }
        e[n-1] = 0.0;

        tq2s240();
        tq2s300();
      }
    return;

    void tq2s240(void) 
      { double f = 0.0;
        double tst1 = 0.0;
        for (int32_t r = 0; r < n; r++)
          { /* Try to compute eigenvalue {d[r]}: */
            double h = fabs(d[r]) + fabs(e[r]);
            if (tst1 < h) tst1 = h;
            /* Look for small sub-diagonal element {e[m]} */
            int32_t m = r;
            while(TRUE)
              { /* Since {e[n-1]} is always zero, this will exit with {m < n}: */
                double tst2 = tst1 + fabs(e[m]);
                if (tst2 == tst1) { break; }
                m++;
              }
            if (m > r) 
              { int32_t iter = 0;
                while (TRUE)
                  { /* Check iteration limit: */
                    if (iter == 30) { (*nev_P) = (uint32_t)r; return; }
                    iter = iter + 1;
                    /* Compute the shift */
                    int32_t r1 = r + 1;
                    int32_t r2 = r1 + 1;
                    double g = d[r];
                    double p = (d[r1] - g) / (2.0 * e[r]);
                    double rr = hypot(p,1.0);
                    d[r] = e[r] / (p + copysign(rr,p));
                    d[r1] = e[r] * (p + copysign(rr,p));
                    double dr1 = d[r1];
                    h = g - d[r];
                    if (r2 < n) { for (int32_t i = r2; i < n; i++) { d[i] = d[i] - h; } }
                    f = f + h;
                    tq2s200(m,r,r1,dr1);
                    double tst2 = tst1 + fabs(e[r]);
                    if (tst2 <= tst1) { break; }
                  }
              }
            d[r] = d[r] + f;
          }
        (*nev_P) = n;
        return;
      }

    void tq2s200(int32_t m, int32_t r, int32_t r1, double dr1)
      { double p = d[m];
        double c = 1.0;
        double c2 = c;
        double er1 = e[r1];
        double s = 0.0;
        double s2, c3;
        int32_t mmr = m - r;
        /* for i=m-1 step -1 until r do -- */
        fprintf(stderr, "!! m=%d r=%d r1=%d dr1=%20.14f\n", m+1, r+1, r1+1, dr1);
        for (int32_t ii = 1; ii <= mmr; ii++) /* 200 */
          { int32_t i = m - ii;
            fprintf(stderr, "!! i=%d p=%20.14f c=%20.14f s=%20.14f\n", i+1, p, c, s);
            c3 = c2;
            c2 = c;
            s2 = s;
            double g = c * e[i];
            double h = c * p;
            double rr = hypot(p,e[i]);
            fprintf(stderr, "!! i=%d s=%20.14f e(i)=%20.14f rr=%24.18f\n", i+1, s, e[i], rr);
            e[i+1] = s * rr;
            s = e[i] / rr;
            c = p / rr;
            p = c * d[i] - s * g;
            d[i+1] = h + s * (c * g + s * d[i]);
            fprintf(stderr, "!! i = %d  d(i+1) =  %20.14f\n", i+1, d[i+1]);
            /* Compute the vector: */
            for (int32_t k = 0; k < n; k++) /* 180 */
              { h = z[(i+1)*(int32_t)n + k];
                z[(i+1)*(int32_t)n + k] = s * z[i*(int32_t)n + k] + c * h;
                z[i*(int32_t)n + k] = c * z[i*(int32_t)n + k] - s * h;
              }
          }
        p = -s * s2 * c3 * er1 * e[r] / dr1;
        e[r] = s * p;
        d[r] = c * p;
        return;
      }

    void tq2s300(void)
      { uint32_t nev = (*nev_P); /* Number of eigenvalues found. */
        fprintf(stderr, "!! nev = %d\n", nev);
        for (int32_t kk = 0; kk < nev; kk++) 
          { fprintf(stderr, "!!   %14.10f\n", d[kk]); }
        /* Selection sort of eigenvalues {d[0..nev-1]}: */
        for (int32_t i = 0; i < nev-1; i++)
          { /* Find smallest element in {d[i]} through {d[nev-1]}: */
            double dmin = d[i]; int32_t imin = i;
            for (int32_t j = i+1; j < nev; j++)
              { if (d[j] < dmin) { dmin = d[j]; imin = j; } }
            assert((imin >= i) && (imin < nev));
            if (imin != i)
              { d[imin] = d[i];
                d[i] = dmin;
                for (int32_t j = 0; j < n; j++) /* 280 */
                  { double zij = z[i*(int32_t)n + j];
                    z[i*(int32_t)n + j] = z[imin*(int32_t)n + j];
                    z[imin*(int32_t)n + j] = zij;
                  }
              }
          }
      }
  }
