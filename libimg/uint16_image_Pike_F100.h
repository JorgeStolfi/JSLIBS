/* uint16_image_Pike_F100.h - read AVC PIKE F-100 raw images. */
/* Last edited on 2017-06-20 20:52:21 by stolfilocal */

#ifndef uint16_image_Pike_F100_H
#define uint16_image_Pike_F100_H

#include <jspnm.h>
#include <uint16_image.h>

uint16_image_t *uint16_image_Pike_F100_read(char *name, bool_t verbose);
  /* Reads file called "{name}", assumed to contain an AVC PIke F-100
    image in RAW16 mode. Returns it as a {uint16_image_t} with 
    {cols = rows = 1000}, {chns = 1}, {maxval = 2^16-1}. 
    If {verbose} is TRUE, prints a notice to {stderr}.
    
    Also reads (and prints, if {verbose} is true) the 
    2944 byte trailer of the image file.
    
    The single channel 0 contain R, G, and B pixels
    interleaved in the 2x2 pattern ((G0,R0),(B0,G1)). */

uint16_image_t *uint16_image_Pike_F100_fread (FILE *rd);
  /* Same as {uint16_image_Pike_F100_read} but reads from the open file {rd}.
    Does not read the 2944 byte trailer of the image file. */

void uint16_image_Pike_F100_read_trailer(FILE *rd, bool_t verbose);
  /* Reads (and prints, if {verbose} is true) the 2944 byte 
    trailer of the image file. */

int uint16_image_Pike_F100_sample_channel(int col, int row);
  /* Returns the channel index (0=R, 1=G, 2=B) of the pixel in column
    {col} and row {row}, according to the Pike F-100 Bayer screen. */
  
uint16_image_t *uint16_image_Pike_F100_debayer(uint16_image_t *img, bool_t squeeze, bool_t verbose);
  /* Applies plain debayering to a Raw16 Pike F-100 image.
  
    The image {img} must have {img.chns==1}. The procedure assumes the
    samples of {img} are R,G,B values interleved in the Bayer pattern
    above. It outputs a three-channel image file {omg} where each
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
  
uint16_image_t *uint16_image_Pike_F100_extract_bayer_channel(uint16_image_t *img, int bayer_col, int bayer_row, bool_t verbose);
  /* Extracts a selected Bayer channel from a Raw16 Pike F-100 image.
    The image {img} must have {img.chns==1}. The procedure assumes the
    samples of {img} are R,G,B values interleved in the Bayer pattern
    above. It outputs a single-channel image file {omg} containing one sample
    from each instance of the Bayer pattern; namely the sample in column
    {bayer_col} and row {bayer_row}, both being either 0 or 1.  The result image
    therefore has half as many rows and columns as the original. */
  
void uint16_image_Pike_F100_color_balance(uint16_image_t *img, float bal[], bool_t verbose);
  /* Multiplies by {bal[k]} the samples in each channel {k} of {img}.
    If {img} has three channels, assumes it has been de-bayerized.
    If {img} has only one channel, assumes the three channels are
    interleaved in the Bayer pattern. */
   
void uint16_image_Pike_F100_output_amp_balance(uint16_image_t *img, int split_col, double gain0, double gain1, bool_t verbose);
  /* Adjusts the samples of {img} to compensate for different gains in
    the output amplifiers of the KAI-1020 sensor, which result in a
    gain discontinuity around column 500. Samples in columns
    {0..split_col-1} are multiplied by {gain0}, samples in {split_col..NX-1} are
    multiplied by {gain1}. Can be used with single- or three-channel
    images, and with images produced by {uint16_image_Pike_F100_extract_bayer_channel}. */
 
uint16_image_t *pnm_Pike_F100_bayer_channel_white_mask(uint16_image_t *img, bool_t verbose);
  /* Returns a white field mask suitable for pixel sensitivity
    correction of a Bayer channel. The image {img} must be a photo of
    a smooth white surface taken with the lens out of focus. It must
    have been produced by {uint16_image_Pike_F100_extract_bayer_channel} and
    must have had its output amplifiers equalized. */

void uint16_image_Pike_F100_interpolate_bayer_pixel(uint16_image_t *img, int col, int row);
  /* Requires {img.chns == 3} and assumes that it is the result of
    {uint16_image_Pike_F100_debayer}. Computes the missing samples of
    pixel in column {col} and row {row} by interpolation of nearby
    ones. */

void uint16_image_Pike_F100_interpolate_bayer(uint16_image_t *img, bool_t verbose);
  /* Applies {uint16_image_Pike_F100_interpolate_bayer_pixel} to all pixels of {img}. */

#endif
