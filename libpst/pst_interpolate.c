/* See {pst_interpolate.h}  */

/* Created on 2011-06-17 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-03-15 23:50:15 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <math.h>
#include <assert.h>
#include <stdint.h>

#include <affirm.h>
#include <float_image.h>
#include <pst_interpolate.h>

void pst_interpolate_values(uint32_t n, double vS[], double wS[], uint32_t m, double *vR_P, double *wR_P)
  {
    double vR = NAN, wR = 0.0;
    if (n > 0)
      { 
        /* Check data: */
        for (int32_t k = 0; k < n; k++)
          { demand(isfinite(wS[k]) && (wS[k] > 0), "invalid or non-posive weight {wS[k]}");
            demand(isfinite(vS[k]), "invalid value {vS[k]}");
          }
        double a[n]; /* Extrapolation/interpolation coefficients */
        switch (n)
          { case 0:
              break;
            case 1: /* Degree 0 extrapolation: */
              a[0] = 1.0; 
              break;
            case 2: /* Affine extrapolation/interpolation: */
              switch (m)
                { case 0: a[0] = +1.5; a[1] = -0.5; break;
                  case 1: a[0] = +0.5; a[1] = +0.5; break;
                  case 2: a[0] = -0.5; a[1] = +1.5; break;
                  default: assert(FALSE);
                }
              break;
            case 3: /* Quadratic extrapolation/interpolation: */
              switch (m)
                { case 0: a[0] = +1.875; a[1] = -1.250; a[2] = +0.375; break;
                  case 1: a[0] = +0.375; a[1] = +0.750; a[2] = -0.125; break;
                  case 2: a[0] = -0.125; a[1] = +0.750; a[2] = +0.375; break;
                  case 3: a[0] = +0.375; a[1] = -1.250; a[2] = +1.875; break;
                  default: assert(FALSE);
                }
              break;
            case 4: /* Cubic extrapolation/interpolation: */
              switch (m)
                { case 0: a[0] = +2.1875; a[1] = -2.1875; a[2] = +1.3125; a[3] = -0.3125; break;
                  case 1: a[0] = +0.3125; a[1] = +0.9375; a[2] = -0.3125; a[3] = +0.0625; break;
                  case 2: a[0] = -0.0625; a[1] = +0.5625; a[2] = +0.5625; a[3] = -0.0625; break;
                  case 3: a[0] = +0.0625; a[1] = -0.3125; a[2] = +0.9375; a[3] = +0.3125; break;
                  case 4: a[0] = -0.3125; a[1] = +1.3125; a[2] = -2.1875; a[3] = +2.1875; break;
                  default: assert(FALSE);
                }
              break;
            default: assert(FALSE);
          }
        double sum_av = 0;
        double sum_a2m = 0;
        for (int32_t k = 0; k < n; k++)
          { sum_av += a[k]*vS[k];
            sum_a2m += a[k]*a[k]*(1.0/wS[k]);
          }
        vR = sum_av;
        wR = 1.0/sum_a2m;

        /* Check for overflow/underflow: */
        if ((! isfinite(vR)) || (! isfinite(wR)) || (wR == 0))
          { vR = NAN; wR = 0.0; }
      }
    (*vR_P) = vR;
    (*wR_P) = wR; 
  }

void pst_interpolate_select_data
  ( int32_t j0,
    pst_iterpolate_get_data_func_t func,
    int32_t *ja_P,
    int32_t *jb_P,
    int32_t *n_P,
    int32_t *m_P
  )
  { /* Find widest {jLO,jHI} preferably centered so that {wP[jLO+1..jHI-1] > 0}: */
    int32_t jLO = j0;
    int32_t jHI = j0 + 1;
    while (jHI - jLO - 1 < 4)
      { double vLO, wLO, vHI, wHI;
        func(jLO, &vLO, &wLO);
        func(jHI, &vHI, &wHI);
        if ((wLO == 0) && (wHI == 0)) { break; }
        if (wLO > 0) { jLO--; }
        if (wHI > 0) { jHI++; }
      }

    int32_t n = jHI - jLO - 1; 
    assert((n >= 0) && (n <= 4));
    int32_t m = j0 - jLO;
    assert((m >= 0) && (m <= n));
      
    (*ja_P) = jLO + 1;
    (*jb_P) = jHI - 1;
    (*n_P) = n;
    (*m_P) = m;
  }
