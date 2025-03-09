/* See pst_map.h */
/* Last edited on 2025-03-01 19:55:47 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <vec.h>
#include <float_image.h>

#include <pst_basic.h>

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

vec_typeimpl(pst_map_vec_t,pst_map_vec,float_image_t *);
