#ifndef float_image_wfilter_H
#define float_image_wfilter_H

/* Tools for image filters based on weighted windows. */
/* Last edited on 2025-02-25 18:32:20 by stolfi */

#include <stdint.h>
#include <float_image.h>

float_image_t* float_image_wfilter_hann(float_image_t *A, int32_t nw);
  /* Creates a copy of image {A}, with same channel, col, and row count,
    and applies to it a 2D smoothing filter with an {nw × nw} Hann
    weight window. */
    
#endif
