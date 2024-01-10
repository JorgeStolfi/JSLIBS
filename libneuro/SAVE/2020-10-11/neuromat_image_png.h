#ifndef neuromat_image_png_H
#define neuromat_image_png_H

/* NeuroMat generic PNG image tools. */
/* Last edited on 2013-12-06 04:51:38 by stolfilocal */

#define _GNU_SOURCE
#include <float_image.h>
  
void neuromat_image_png_write(char *prefix, char *tag, float_image_t *fim, double vlo, double vhi);
  /* Writes the image {fim} to a file named "{prefix}_{tag}.png" in PNG
    format. Samples are affinely mapped from {vlo} to 0 and from {vhi]}
    to 1, then gamma-mapped (gamma=2.5) and quantized. */

#endif
