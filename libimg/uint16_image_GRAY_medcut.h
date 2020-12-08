/* uint16_image_GRAY_medcut.h - choose a representative set of N grays, from gray histogram */
/* Last edited on 2006-11-10 19:25:29 by stolfi */ 

#ifndef uint16_image_GRAY_medcut_H
#define uint16_image_GRAY_medcut_H

#include <jspnm.h>

typedef uint16_t *pgm_pixel_vector;

extern pgm_pixel_vector
uint16_image_GRAY_median_cut(
    long *gh,              /* Pixel count indexed by gray level */
    uint16_t maxval,   /* maxval of pixel values, also size of "gh" */
    int *newgraysp         /* In: desired number of grays, out: number chosen */
  );
  /*
    Chooses a good set of gray values for image quantization,
    given the image's histogram.
    
    Based on Paul Heckbert's median cut algorithm, as described in
    "Color Image Quantization for Frame Buffer Display", SIGGRAPH '82
    Proceedings, page 297.
  */

#endif
