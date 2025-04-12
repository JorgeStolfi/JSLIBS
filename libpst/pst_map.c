/* See pst_map.h */
/* Last edited on 2025-03-15 23:47:13 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <vec.h>
#include <float_image.h>

#include <pst_basic.h>
#include <pst_interpolate.h>

#include <pst_map.h>
    
void pst_map_ensure_pixel_consistency(float_image_t *A, int32_t wch)
  { int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);

    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float vA[NC];
            float_image_get_pixel(A, x, y, vA);
            pst_ensure_pixel_consistency(NC, wch, FALSE, vA);
            float_image_set_pixel(A, x, y, vA);
          }
      }
  }

void pst_map_interpolate_samples
  ( float_image_t *A,
    int32_t c,
    int32_t wch,
    int32_t x0, int32_t y0,
    int32_t x1, int32_t y1,
    bool_t extrapolate,
    double *vR_P, double *wR_P
  )
  {
    int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);
    demand((c >= 0) || (c < NC), "invalid channel index {c}");

    int32_t dx = x1 - x0, dy = y1 - y0;
    demand(((dx == 0) && (dy == 1)) || ((dx == 1) && (dy == 0)), "pixels are not consecutive");

    auto void fetch(int32_t j, double *vR_P, double *wR_P);
      /* Fetches the value {v = A[c,x,y]} and the weight {w=A[2,x,y]},
        where {[x,y]} is {j} pixels ahead or behind {[x0,y0]} in the
        direction of {[x1,y1]}. Namely, {[x,y]} is {[x0,y0]} if {j=0},
        {[x1,y1]} if {j=1}, and generally {[x0+j*dx,y0+j*dy]}. However
        if pixel {x,y} does not exist, or {v} is not finite, of {w} is
        zero, sets {v=NAN} and {w=0}. Returns {v} in {*vR_P} and {w} in
        {*wR_P}. */

    int32_t ja, jb, n, m;
    pst_interpolate_select_data(0, fetch, &ja, &jb, &n, &m);
    
    /* Got enough data? */
    if (((m >= 1) && (n - m >= 1)) || (extrapolate && (n >= 3)))
      { /* Collect values to interpolate/extrapolate: */
        double vS[n], wS[n];
        for (int32_t k = 0; k < n; k++)
          { fetch(ja + k, &(vS[k]), &(wS[k]));
            assert(isfinite(vS[k]));
            assert(isfinite(wS[k]) && (wS[k] > 0));
          }
        /* Interpolate/extrapolate them: */
        pst_interpolate_values((uint32_t)n, vS, wS, (uint32_t)m, vR_P, wR_P);
      }
    else
      { /* Data is not enough or not well placed: */
        (*vR_P) = NAN; (*wR_P) = 0;
      }

    return;

    void fetch(int32_t j, double *vR_P, double *wR_P)
      { double vR = NAN, wR = 0.0;
        int32_t x = x0 + j*dx;
        int32_t y = y0 + j*dy; 
        if ((x >= 0) && (y >= 0) && (x < NX) && (y < NY))
          { wR = ((wch >= 0) && (wch < NC) ? float_image_get_sample(A,wch,x,y) : 1.0); 
            assert(isfinite(wR) && (wR >= 0));
            if (wR > 0)
              { vR = float_image_get_sample(A,c,x,y);
                if (! isfinite(vR)) { vR = NAN; wR = 0; }
              }
          }
        (*vR_P) = vR;
        (*wR_P) = wR;
      }
  }

vec_typeimpl(pst_map_vec_t,pst_map_vec,float_image_t *);
