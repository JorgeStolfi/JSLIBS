/* Stack of images with limited depth of focus. */
/* Last edited on 2024-08-02 15:42:20 by stolfi */

#ifndef multifok_result_H
#define multifok_result_H

#define _GNU_SOURCE
#include <stdint.h>

#include <float_image.h>
   
typedef struct mfok_result_t
  { float_image_t *timg; 
    float_image_t *azimg;
    float_image_t *dzimg;
    double noise;
  } mfok_result_t;
  /* The result of multifocus stereo: the composite sharp image {timg},
    the estimated height map {azimg}, and the estimated height uncertainly 
    map {dzimg}. */
  

#endif