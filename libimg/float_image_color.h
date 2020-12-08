#ifndef float_image_color_H
#define float_image_color_H

/* Tools specific for color images. */
/* Last edited on 2008-07-29 20:35:32 by stolfi */ 

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <float_image.h>
#include <frgb.h>

frgb_t fic_get_frgb_pixel(float_image_t *A, int cR, int cG, int cB, int x, int y);
/* Extracts the samples of channels {cR}, {cB} and {cG} of the pixel
  at column {x} and row {y} of {A} as an RGB color triple.
  
  Fails if {x} or {y} are outside their valid ranges or if {cR} is
  outside the range {0..A->NC-1}. If {cG} is outside that range, the
  procedure assumes {cB == cR}. Ditto for {cB}. In particular, if
  {cR=0,cG=1,cB=2}, but {A} is a monochromatic image, the sample at
  {x,y} is interpreted as a gray value and converted to an RGB
  triple. */

void fic_set_frgb_pixel(float_image_t *A, int cR, int cG, int cB, int x, int y, frgb_t *p);
/* Sets the samples of channels {cR}, {cB} and {cG} of the pixel at
  column {x} and row {y} of {A} to the components of the RGB triple {p}.

  Fails if {x} or {y} are outside their valid ranges, or if any of {cR},
  {cG}, or {cB} is outside the range {0..A->NC-1}. */

void fic_normalize_colors(float_image_t *A, int cR, int cG, int cB);
  /* Maps the colors of all pixels to fit in the [0..1] cube, by scaling 
    uniformly the three components and and adding or subtracting a constant
    amount of white.  This transformation preserves hue, and
    rescales brightness and chroma independently of each other.
    So grays remain gray, and equally bright colors remain equally 
    bright. */

#endif
