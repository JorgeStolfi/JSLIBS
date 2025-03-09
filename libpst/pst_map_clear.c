/* See pst_map_clear.h */
/* Last edited on 2025-03-06 13:04:23 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>

#include <pst_map_clear.h>

void pst_map_clear
  ( float_image_t *A,
    int32_t wch,
    int32_t xmin,
    int32_t xmax,
    int32_t ymin,
    int32_t ymax
  )
  { int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY); 
    
    if (xmin < 0) { xmin = 0; }
    if (xmax >= NX) { xmax = NX-1; }
    if (ymin < 0) { ymin = 0; }
    if (ymax >= NY) { ymax = NY-1; }
    
    float v[NC];
    for (int32_t c = 0; c < NC; c++) { v[c] = (c == wch ? 0.0 : NAN); }
    float_image_fill_rectangle_pixels(A, xmin, xmax, ymin, ymax, v);
  }
