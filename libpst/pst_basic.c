/* See pst_basic.h */
/* Last edited on 2025-03-03 19:35:26 by stolfi */

#include <math.h>
#include <assert.h>
#include <values.h>
#include <stdlib.h>
#include <stdint.h>

#include <float_image.h>
#include <vec.h>
#include <r2.h>
#include <r3.h>
#include <argparser.h>

#include <pst_basic.h>

void pst_perturb_normal(r3_t *nrm, double amt)
  { /* Generate a Gaussian random vector with mean 0, deviation {1} in each coordinate: */
    r3_t rnd; 
    r3_throw_normal(&rnd);
    /* Scale {rnd} by {1/sqrt(2)} so that a small {amt} gives
      an rms displacement of {amt}, then mix it into {nrm}
      with weights {amt}:{1-amt}: */
    double crnd = amt/M_SQRT2;
    double cnrm = 1 - amt;
    r3_mix(crnd, &rnd, cnrm, nrm, nrm);
    /* Re-normalize: */
    double len = r3_dir(nrm, nrm);
    /* If we were so unlucky as to get a near-zero mix, try again: */
    while (len < 0.00001) { r3_throw_normal(nrm); len = r3_dir(nrm, nrm); }
  }

r3_t pst_normal_from_slope(r2_t *grd)
  { double dZdX = grd->c[0]; 
    double dZdY = grd->c[1]; 
    double m = sqrt(1.0 + dZdX*dZdX + dZdY*dZdY);
    double nx = dZdX/m; 
    double ny = dZdY/m; 
    double nz = 1.0/m; 
    return (r3_t){{ nx, ny, nz }};
  }
  
r2_t pst_slope_from_normal(r3_t *nrm, double maxSlope)
  { double nx = nrm->c[0], ny = nrm->c[1], nz = nrm->c[2];
    double nr = hypot(nx, ny);
    double nzmin = nr / maxSlope;
    if (nz < nzmin) { nz = nzmin; }
    double dZdX = -nx / nz;
    double dZdY = -ny / nz;
    return (r2_t) {{ dZdX, dZdY }};
  }

void pst_ensure_pixel_consistency(int32_t NC, int32_t wch, bool_t bad, float v[])
  { 
    float w = ((wch >= 0) && (wch < NC) ? v[wch] : 1.0f);
    demand(isfinite(w) && (w >= 0), "invalid pixel weight");
    bad = bad | (w == 0.0);
    for (int32_t c = 0; c < NC; c++) { if (! isfinite(v[c])) { bad = TRUE; } }
    for (int32_t c = 0; c < NC; c++) 
      { if (bad)
          { v[c] = (c == wch ? 0.0 : NAN); }
        else
          { assert(isfinite(v[c]) && ((c != wch) || (v[c] > 0))); }
      }
  }
