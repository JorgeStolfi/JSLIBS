/* uint16_image_check_dither.h - Checks whether a digital image is a dither matrix. */
/* Last edited on 2024-12-05 10:31:04 by stolfi */

#ifndef uint16_image_check_dither_H
#define uint16_image_check_dither_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <uint16_image.h>
     
bool_t uint16_image_check_dither(uint16_image_t *img, bool_t die);
  /* Checks whether {img} is a valid a dither matrix. Namely,
    the {maxval} must be the number of pixels minus one,
    and, in each channel, each sample value in {0..maxval}
    must occur exactly once.
    
    If the image is a valid dither matrix, the procedure returns true.
    Otherwise, the procedure prints an error message to {stderr},
    and either returns false, if {die} is false, or aborts, if {die} 
    is true. */

#endif
