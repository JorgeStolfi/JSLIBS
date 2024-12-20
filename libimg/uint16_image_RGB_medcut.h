/* uint16_image_RGB_medcut.h - choose a representative set of N colors, from color histogram */
/* Last edited on 2024-12-05 07:29:27 by stolfi */ 

#ifndef uint16_image_RGB_medcut_H
#define uint16_image_RGB_medcut_H

#include <jspnm.h>
#include <uint16_image_RGB_hist.h>

#define MAXCOLORS 32767

uint16_image_RGB_hist_vector uint16_image_RGB_median_cut
  ( uint16_image_RGB_hist_vector ch,  /* Color histogram of image. */
    int32_t colors,           /* Size of {ch}. */
    uint16_t maxval,  /* Max valid sample value. */
    uint32_t *newcolorsp       /* In: desired number of colors, out: number chosen */
  );
  /* Chooses a good set of pixel values for color image quantization,
    given the image's histogram.
    
    Based on Paul Heckbert's median cut algorithm, as described in
    "Color Image Quantization for Frame Buffer Display", SIGGRAPH '82
    Proceedings, page 297. */

#endif
