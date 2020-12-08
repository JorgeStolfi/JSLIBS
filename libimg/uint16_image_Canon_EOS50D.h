/* uint16_image_Canon_EOS50D.h - read AVC PIKE F-100 raw images. */
/* Last edited on 2017-06-20 20:51:56 by stolfilocal */

#ifndef uint16_image_Canon_EOS50D_H
#define uint16_image_Canon_EOS50D_H

#include <jspnm.h>
#include <uint16_image.h>


uint16_image_t *pnm_Canon_EOS50D_raw16_debayer(uint16_image_t *img,bool_t squeeze, bool_t verbose) ;
  /* Applies plain debayering to a Raw16 Canon EOS50D image.
  
    The image {img} must have {img.chns==1}. The procedure assumes the
    samples of {img} are R,G,B values interleved in the Bayer pattern
    ((R0,G0),(G1,B1)). It outputs a three-channel image file {omg} where each
    sampleof {img} has been copied into the appropriate channel.
    
    If {squeeze} is false, the output image has the same size as the original,
    and each sample is copied to the appropriate channel without
    changing its row or column. The other samples in {omg} are set to zero.
    Thus channels 0 and 2 will have only one significant sample in every 2x2 
    block of pixels, while channel 1 will have two significant samples per  
    block, in a checkerboard fashion.
    
    If {squeeze} is true, the output has half as many columns and rows
    as {img}. The sample in column {col} and row {row} of {img} is
    copied to the appropriate channel of column {col/2} and row
    {row/2}, rounded down. The two green Bayer channels G0 and G1 are
    averaged to form channel 1 of the output. The average is rounded
    down or up, in checkerboard fashion. */


#endif
