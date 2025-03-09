/* See pst_map_shift_values.h */
/* Last edited on 2025-02-27 14:37:17 by stolfi */

#include <math.h>
#include <assert.h>
#include <stdint.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <float_image.h>
#include <pst_basic.h>

#include <pst_map_shift_values.h>

void pst_map_shift_values(float_image_t *A, int32_t wch, sign_t dir, double shift[])
  { 
    int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);
    
    for (int32_t y = 0;  y < NY; y++)
      { for (int32_t x = 0;  x < NX; x++)
          { float vA[NC];
            float_image_get_pixel(A, x, y, vA);
            double wA = ((wch >= 0) && (wch < NC) ? vA[wch] : 1.0);
            demand(isfinite(wA) && (wA >= 0), "invalid weight value");
            for (int32_t c = 0; c < NC; c++)
              { if (c != wch)
                  { if (! isfinite(vA[c])) 
                      { wA = 0; }
                    else 
                      { vA[c] = (float)(vA[c] + dir*shift[c]);
                        if (! isfinite(vA[c])) { wA = 0; }
                      }
                  }
              }
            pst_ensure_pixel_consistency(NC, wch, (wA == 0), vA);
            float_image_set_pixel(A, x, y, vA);
          }
      }
  }    
