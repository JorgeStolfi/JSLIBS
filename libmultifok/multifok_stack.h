/* Stack of images with limited depth of focus. */
/* Last edited on 2024-05-02 10:43:24 by stolfi */

#ifndef multifok_stack_H
#define multifok_stack_H

#define _GNU_SOURCE
#include <stdint.h>

#include <float_image.h>
  
typedef struct multifok_image_stack_t 
  { int32_t NI;              /* Number of images in stack. */
    int32_t NX, NY;          /* Image dimensions of all images. */
    int32_t NC;              /* Number of channels in the scene images {csimg[ki]}. */
    float_image_t *csimg[];  /* Images of scene. */
    float_image_t *shimg[];  /* Sharpness maps for the images. */
    double zFoc[];           /* Assumed height of focus plane in each image. */
  } multifok_image_stack_t;
  /* A stack of input scene images {csimg[ki]} for {ki} in {0..NI-1},
    the corresponding sharpness maps {shimg[ki]}, and the corresponding
    focus plane heights {zFoc[ki]}. All images must have the 
    same number of columns {NX} and of rows {NY}.  All scene images {csimg[ki]}
    must have the same number of channels {NC}. */
     
multifok_image_stack_t* mfmi_read_stack
  ( char *inPrefix,
    int32_t NI,
    double *zFoc
  );
  /* Reads a stack of scene images {csimg[0..NI-1]} and the corresponding 
    sharpness score maps {shimg[0..NI-1]}, from files 
    "{inprefix}{NNNNNNN}", and returns them in a {multifok_image_stack_t}
    object {stack}.  Also sets {stac.zFoc} to the vector {zFoc},
    assumed to have {NI} elements. */

#endif
