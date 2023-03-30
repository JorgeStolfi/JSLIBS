/* See {float_image_hartley_spectrum.h}. */
/* Last edited on 2023-03-18 17:06:18 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>

#include <fftw3.h>
 
#include <bool.h>
#include <affirm.h>
#include <float_image.h>

#include <float_image_hartley_spectrum.h>

void float_image_hartley_spectrum(float_image_t *H, float_image_t *P, bool_t center)
  {
    int32_t NC, NX, NY;
    float_image_get_size(H, &NC, &NX, &NY);
    int32_t HX = NX/2;
    int32_t HY = NY/2;
    
    for (int32_t ic = 0; ic < NC; ic++)
      { for (int32_t fy = 0; fy < NY; fy++)
          { int32_t gy = (NY - fy) % NY;
            /* Need to scan only half of Fourier transform: */
            for (int32_t fx = 0; fx <= HX; fx++)
              { int32_t gx = (NX - fx) % NX;
                /* Compute power at frequency {fx,xy} and {gx,gy}: */
                double cfxy = float_image_get_sample(H, ic, fx, fy);
                double pxy = cfxy*cfxy;
                if ((fx != gx) || (fy != gy))
                  { double cgxy = float_image_get_sample(H, ic, gx, gy);
                    pxy += cgxy*cgxy;
                  }
                /* Save in power spectrum image, with requested shift: */
                int32_t rx, ry, sx, sy; /* Indices where to put {fx,fy,gx,gy}: */
                if (center)
                  { rx = (fx + HX) % NX;
                    ry = (fy + HY) % NY;
                    sx = (gx + HX) % NX;
                    sy = (gy + HY) % NY;
                  }
                else
                  { rx = fx; ry = fy; sx = gx; sy = gy; }
                if ((fx == gx) && (fy == gy))
                  { float_image_set_sample(P, ic, rx, ry, (float)pxy); }
                else
                  { float_image_set_sample(P, ic, rx, ry, (float)pxy/2);
                    float_image_set_sample(P, ic, sx, sy, (float)pxy/2);
                  }
              }
          }
       }
  }
