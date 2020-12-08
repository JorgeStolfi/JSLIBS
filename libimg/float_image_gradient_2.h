#ifndef float_image_gradient_2_H
#define float_image_gradient_2_H

/* Extra tools for image gradients. */
/* Last edited on 2009-07-03 17:02:24 by stolfi */ 

#include <bool.h>
#include <r2.h>
#include <gauss_table.h>
#include <float_image.h>

void float_image_gradient_sqr_relative_2
  ( float_image_t *A,
    int cA,
    double noise, 
    float_image_t *G,
    int cG    
  );
  /* Computes the relativized gradient squared of image {A}, stores it into channel {cG}
    of image {G} (which must have the same dimensions as {A}).
  
    If {cA} is non-negative, the sample of {G} at pixel {p} and channel
    {cG} is {g2/(v2+noise^2)}, where {g2} is the squared {3×3} Sobel
    gradient of channel {cA} of {A} at {p}, and {v2} is the local variance of 
    channel {cA} of {A} within a {5×5} Gaussian-like window.
    
    If {cA} is negative, the result is the arithmetic average of
    {g2/(v2 + noise^2)} over all channels of {A}. */

#endif
