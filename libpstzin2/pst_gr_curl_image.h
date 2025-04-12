/* Last edited on 2025-03-14 19:07:12 by stolfi */
/* Fill an image with the curl of a height difference graph. */

#ifndef pst_gr_curl_image_H
#define pst_gr_curl_image_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <float_image.h>

#include <pst_gr.h>

void pst_gr_curl_image_fill(pst_gr_t *gr, float_image_t* U);
  /* Stores into channel 0 of the image {U} the curl of the faces of {gr}.
    
    Specifically, the image {U} should have one column and one row less
    than the height map associated with {gr}. Each sample {U[0,x,y]}
    is set to the curl of the face that contains the center of the pixel. */

#endif
