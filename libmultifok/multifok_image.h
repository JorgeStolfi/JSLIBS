/* Test tools for {multifok_focus_op} and related funcs. */
/* Last edited on 2025-04-11 08:40:03 by stolfi */

#ifndef multifok_image_H
#define multifok_image_H

#include <stdint.h>

#include <bool.h>
#include <i2.h>
#include <float_image.h>

#define multifok_image_gamma 1.000
  /* Assumed encoding gamma of input and output PNG images. */

#define multifok_image_bias 0.0327
  /* Assumed encoding bias of input and output PNG images (irrelevant if gamma is 1). */

float_image_t *multifok_image_set_weight_channel(float_image_t *img, float_image_t *wht, int32_t chns, bool_t clear);
  /*  The image {img} must have {chns} or {chns+1} channels, where channels {0..chns-1}
    are data and channel {chns}, if it exists, gives the reliability weight of that
    data.  
    
    The procedure returns a newly allocated image {timg} with same size
    as {img} and exactly {chns+1} channels. The data channels
    {0..chns-1} of {timg} are copied from {img}. If {wht} is not null,
    the weight channel {chns} of {timg} is taken from channel 0 of
    {wht}. If {wht} is null, channel {chns} of {timg} is copied from
    channel {chns} of {img} if it exists, otherwise set to all ones.
    
    In any case, if {wht} is not null, must have the same col and row
    counts as {img} and a single channel.
    
    If {clear} is true, the procedure will also set to {NAN} the data samples
    of each channel that ends up with zero weight. */
 
void multifok_image_draw_crosses(float_image_t *img, int32_t ch, uint32_t NQ, i2_t pix[], float val);
  /* Draws open crosses into image {img} at the positions listed in {pix[0..NQ-1]}
    The crosses are drawn with value {val} into channel {ch}, and with value 0 in every other
    channels. */

#endif
