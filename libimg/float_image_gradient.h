#ifndef float_image_gradient_H
#define float_image_gradient_H

/* Tools for image gradients. */
/* Last edited on 2012-01-10 04:12:23 by stolfilocal */ 

#include <bool.h>
#include <r2.h>
#include <gauss_table.h>
#include <float_image.h>

void float_image_gradient_sobel(float_image_t *A, int cA, float_image_t *DX, int cX, float_image_t *DY, int cY);
  /* Computes the gradient image of channel {cA} image {A}. If {DX} is
    not null and {cX} is non-negative, stores the X derivative into
    channel {cX} of image {DX} If {DY} is not null and {cY} is
    non-negative, stores the {Y} derivative into channel {cY} of image
    {DY}. The derivatives are computed using the {3×3} Sobel
    operator. */

void float_image_gradient_sqr_sobel(float_image_t *A, int cA, float_image_t *G, int cG);
  /* Computes the squared gradient image of channel {cA} image {A},
    stores it into channel {cG} of image {G}.
  
    More precisely, stores into each pixel {p} of channel {cG} of {G}
    the value {fx^2+fy^2}, where {fx} is the X derivative of channel
    {cA} of {A} at {p} and {fy} is the Y derivative at {p}. The
    derivatives are computed using the {3×3} Sobel operator. */

float_image_t *float_image_gradient_sqr_relative
  ( float_image_t *A, 
    double noise, 
    bool_t mono
  );
  /* Computes the relativized gradient squared {R} of image {A}.
  
    If {mono} is FALSE, the image {R} has the same dimensions and
    channel count as {A}. The sample of {R} at pixel {p} and channel
    {c} is {g2/(v2+noise^2)}, where {g2} is the squared {3×3} Sobel
    gradient of {A} at {p}, and {v2} is the local variance of {A}
    within a {5×5} Gaussian-like window.
    
    If {mono} is TRUE, the image {R} has a single channel, whose
    value is the arithmetic average of {g2/(v2 + noise^2)} 
    over all channels of {A}. */

#endif
