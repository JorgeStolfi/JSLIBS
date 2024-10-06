/* Stack of images with limited depth of focus. */
/* Last edited on 2024-09-19 12:57:30 by stolfi */

#ifndef multifok_image_stack_H
#define multifok_image_stack_H

#define _GNU_SOURCE
#include <stdint.h>

#include <float_image.h>
  
typedef struct multifok_image_stack_t 
  { int32_t NI;              /* Number of images in stack. */
    int32_t NC;              /* Number of channels in the scene images {csimg[ki]}. */
    int32_t NX, NY;          /* Image dimensions of all images. */
    float_image_t **csimg;   /* Images of scene. */
    float_image_t **shimg;   /* Sharpness maps for the images. */
    double *zFoc;            /* Assumed height of focus plane in each image. */
  } multifok_image_stack_t;
  /* A stack of input scene images {csimg[ki]} for {ki} in {0..NI-1},
    the corresponding sharpness maps {shimg[ki]}, and the corresponding
    focus plane heights {zFoc[ki]}. All images must have the 
    same number of columns {NX} and of rows {NY}.  All scene images {csimg[ki]}
    must have the same number of channels {NC}. */
     
multifok_image_stack_t* multifok_image_stack_read
  ( char *inPrefix,
    double zDep,
    int32_t NI,
    double zFoc[]
  );
  /* Reads a stack of scene images {csimg[ki]} and the corresponding
    sharpness score maps {shimg[ki]} from files
    "{inprefix}-{tag}-cs.png" and "{inprefix}-{tag}-sh.png",
    respectively; where {ki} ranges in {0..NI-1}, {tag} is
    "fd{DDD.DDDD}-z{FFF.FFFF}", {DDD.DDDD} is the value of {Zdep}, and
    {FFF.FFFF} is the value of {zFoc[ki]}.
    
    Returns those images in a {multifok_image_stack_t} object
    {stack}. */

#endif
