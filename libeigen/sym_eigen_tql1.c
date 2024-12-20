/* See sym_eigen_tql1.h */
/* Last edited on 2024-12-05 18:44:20 by stolfi */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>

#include <sym_eigen_tql1.h>

#define Pr fprintf
#define Er stderr
  /* Shoter names. */

void sym_eigen_tql1(uint32_t n, double d[], double e[], uint32_t *nev_P)
  {
    /* The algorithm is based on the EISPACK routine "tql1.f" 
      by Burton S. Garbow, Argonne National Laboratory
      (aug/1983); orignally from the ALGOL procedures {tql1} and
      {tql2}, Num. Math. 11, 293-306(1968) by Bowdler, Martin,
      Reinsch, and Wilkinson. See Handbook for Auto. Comp., vol.II -
      Linear Algebra, 227-240(1971). Re-implemented in C by Jorge
      Stolfi, Unicamp (dec/2002). */
      
    auto void tq1s290(void);
    
    auto void tq1s200(int32_t r, int32_t r1, int32_t m, double dr1);
      /* The QL transformation. */
      
    auto void tq1s230(double dr, int32_t r);
      /* Inserts {dr} in the proper place in {d[0]} thtough 
        {d[r-2]}. Does not increment {r} */
        
    if (n == 0)
      { /* Void: */
        (*nev_P) = 0;
      }
    else if (n == 1)
      { /* Trivial, {d[0]} is the eigenvalue: */
        (*nev_P) = 1;
      }
    else
      { (*nev_P) = UINT32_MAX; /* To be redefined. */
        /* Shift the {e} vector down, set {e[n-1] = 0}: */
        for (int32_t i = 1; i < n; i++) { e[i-1] = e[i]; }
        e[n-1] = 0.0;
        /* Do all the work: */
        tq1s290();
      }
    return;

    void tq1s290(void)
      { 
        double f = 0.0;
        double tst1 = 0.0;
 
        for (int32_t r = 0; r < n; r++) /* 290 */
          { /* Try to compute eigenvalue {d[r]}: */
            double h = fabs(d[r]) + fabs(e[r]);
            if (tst1 < h) { tst1 = h; }
            /* Look for small sub-diagonal element {e[m]}: */
            int32_t m = r;
            while (TRUE)
              { /* Since {e[(n-1)] == 0}, this will exit with {m < n}: */
                double tst2 = tst1 + fabs(e[m]);
                if (tst2 == tst1) { break; }
                m++;
              }
            assert(m < n);
            if (m > r)
              { int32_t iter = 0;
                while (TRUE)
                  { /* Check for iteration limit: */
                    if (iter >= 30) { (*nev_P) = (uint32_t)r; return; }
                    iter = iter + 1;
                    /* Compute shift: */
                    int32_t r1 = r + 1;
                    int32_t r2 = r1 + 1;
                    double g = d[r];
                    double p = (d[r1] - g) / (2.0 * e[r]);
                    double rr = hypot(p, 1.0);
                    d[r] = e[r] / (p + copysign(rr, p));
                    d[r1] = e[r] * (p + copysign(rr, p));
                    double dr1 = d[r1];
                    h = g - d[r];
                    if (r2 < n) { for (int32_t i = r2; i < n; i++) { d[i] = d[i] - h; } }
                    f = f + h;
                    tq1s200(r,r1,m,dr1);
                    double tst2 = tst1 + fabs(e[r]);
                    if (tst2 <= tst1) { break; }
                  }
              }
            double dr = d[r] + f;
            tq1s230(dr,r);/* 290 */
         }
         (*nev_P) = n;
         return;
      }

    void tq1s200(int32_t r, int32_t r1, int32_t m, double dr1)
      {
        double p = d[m];
        double c = 1.0;
        double c2 = c;
        double s2, c3;
        double er1 = e[r1];
        double s = 0.0;
        int32_t mmr = m - r;
        /* for i=m-1 step -1 until r do -- */
        for (int32_t ii = 1; ii <= mmr; ii++)
          { int32_t i = m - ii;
            c3 = c2;
            c2 = c;
            s2 = s;
            double g = c * e[i];
            double h = c * p;
            double rr = hypot(p,e[i]);
            e[i+1] = s * rr;
            s = e[i] / rr;
            c = p / rr;
            p = c * d[i] - s * g;
            d[i+1] = h + s * (c * g + s * d[i]);
          }

        p = -s * s2 * c3 * er1 * e[r] / dr1;
        e[r] = s * p;
        d[r] = c * p;
        return;
      }

    void tq1s230(double dr, int32_t r)
      { /* Bubble {dr} among {d[0..r-1]}: */
        assert((r >= 0) && (r < n));
        int32_t i = r;
        while ((i >= 1) && (dr < d[i-1])) { d[i] = d[i-1]; i--; }
        assert((i >= 0) && (i <= r));
        d[i] = dr;
        return;
      }
  }
