/* uint16_image_RGB_hist.h - color histogram tools */
/* Last edited on 2009-01-07 01:38:23 by stolfi */ 

/* Sub-clone of Jef Poskanzer's ppmcmap.h */

#ifndef uint16_image_RGB_hist_H
#define uint16_image_RGB_hist_H

#include <jspnm.h>

typedef struct ppm_pixel_t { uint16_t c[3]; } ppm_pixel_t;
  /* An RGB pixel. */

typedef struct uint16_image_RGB_hist_item { ppm_pixel_t color; int value; } uint16_image_RGB_hist_item;
  /* An entry of an RGB histogram. */

typedef struct uint16_image_RGB_hist_item* uint16_image_RGB_hist_vector;

uint16_image_RGB_hist_vector uint16_image_RGB_hist_build 
  ( uint16_t** samples, 
    int chns,
    int cols, 
    int rows,
    int maxcolors, 
    int* colorsP
  );
  /* Returns a {uint16_image_RGB_hist_vector} with space allocated for {maxcolors} entries,
    where the first {*colorsP} entries are taken from the given samples.
    Assumes that each row consists of {cols} pixels, each with {chns} samples, 
    all consecutive. */

int ppm_equal(ppm_pixel_t *a, ppm_pixel_t *b);
  /* TRUE if the two pixels {a,b} are identical. */

void uint16_image_RGB_hist_add
  ( uint16_image_RGB_hist_vector chv, 
    int* colorsP, 
    int maxcolors, 
    ppm_pixel_t* colorP,
    int value, 
    int position
  );

void uint16_image_RGB_hist_free (uint16_image_RGB_hist_vector chv);

#endif
