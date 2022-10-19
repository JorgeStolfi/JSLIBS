#ifndef neuromat_image_png_H
#define neuromat_image_png_H

/* NeuroMat generic PNG image tools. */
/* Last edited on 2021-08-29 01:05:19 by stolfi */

#define _GNU_SOURCE
#include <float_image.h>
  
void neuromat_image_png_write(char *dir, char *name, float_image_t *fim, float vlo, float vhi, double gamma);
  /* Writes the image {fim} to a file named "{dir}/{name}.png" in PNG
    format. The number of channels must be 1 (grayscale), 2 (gray+opacity), 3 (RGB), or 4 (RGB + Opacity).
    Samples (other than opacity samples) are affinely mapped from {[vlo__vhi]} to {[0__1]},]
    encoded with the given {gamma}, and quantized. */

#endif
