#ifndef float_image_overlay_H
#define float_image_overlay_H

/* Tools for overlaying float images according to opacity channels. */
/* Last edited on 2021-08-28 03:46:13 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <float_image.h>

/* OPACITY CHANNEL

  The procedures in this interface assume that a specified channel
  {icop} of an argument image is an /opacilty channel/ that determines
  the value other channels when the image is /ovrlaid/ on different
  backgrounds.
  
  Let {oA} be the sample on the opacity channel {icop} of a pixel of
  image {A}. The sample {oA} must be between 0 and 1. These procedures
  assume that a random point {p} inside that pixel has probability {oA}
  of being covered by {A}, so that the samples in other channels are
  determined by {A}, and probability {1-oA} of being not covered at all
  by {A}, so that the samples in other channels are determined by
  whatever in behind {A}.
  
  Thus, suppose {A} is overlaid on an background {Z} that is assumed to
  be fully opaque.  If {vA} is the sample in any other channel of that pixel of {A},
  and {vZ} is the sample of that same channel and pixel of {Z},
  the average value of the combination, in that channel and over that pixel,
  is by definition {oA*vA + (1-oA)*vZ}.
  
  It follows that the value {vA} is assumed irrelevant if  {oA} is zero,
  that is, if {A} is fully transparent in that pixel.  
  
  Moreover, two images {A,B} with opacity channels are overlaid, the
  events "point {p} is covered by {A}" and "point {p} is colvered by
  {B}" are assumed to be independent (even if {A} and {B} are copies of
  the same image, or were computed by the same algorithm). */
   
void float_image_overlay(float_image_t *A, float_image_t *B, int32_t icop, int32_t xlo, int32_t ylo);
  /* Overlays image {B}, displaced by {xlo} columns and {ylo} rows,  on top of image {A},
    assuming that channel {icop} of both is the opacity channel.
 
    The two images must have the same number of channels, and {icop}
    must be a valid channel index. The procedure ignored any pixels of
    {B} that, after displacement, lie outside the domain of {A}. It does
    not change any pixels of {A} that are not within the displaced
    domain of {B}
    
    For each pixel, let {oA,oB} be the sample values of {A} and {B} from
    that pixel channel {icop}, and {vA,vB} be the samples from that
    pixel but any other channel. The procedure sets {vA} to
    be {(oB*vB + (1-oB)*oA*vA)/oC} where {oC = 1 - (1-oA)*(1-oB) = oB +
    (1-oB)*oA}; then sets {oA} to {oC}.
    
    If {oC} is zero then the value of {vA} is not changed, and {oA} is
    still zero. This occurs when both images are effectively fully
    transparent at that pixel {oA=oB=0}, but may also occur also because
    of floating-point roundoff when {oA} and {oB} are too small. */

#endif
