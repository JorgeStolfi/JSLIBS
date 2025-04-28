/* Last edited on 2025-04-27 02:56:33 by stolfi */
/* Fill an image with the curl of a height difference graph. */

#ifndef psgr_curl_image_H
#define psgr_curl_image_H

/* Created by Rafael F. V. Saracchini */

#include <stdint.h>

#include <float_image.h>

#include <psgr_types.h>
#include <psgr.h>

void psgr_curl_image_fill(psgr_t *gr, float_image_t* U);
  /* Stores into channel 0 of the image {U} the curl of the faces of {gr}.
    
    Specifically, the image {U} should have one column and one row less
    than the height map associated with {gr}. Each sample {U[0,x,y]}
    is set to the curl of the face that contains the center of the pixel. */

#endif
