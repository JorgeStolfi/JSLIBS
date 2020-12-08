/* uint16_image_read_jpeg.h - routines to read images of {uint16_t} samples from JPEG files. */
/* Last edited on 2017-06-23 01:44:24 by stolfilocal */

/* Created by R. Minetto (IC-UNICAMP) sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#ifndef uint16_image_read_jpeg_H
#define uint16_image_read_jpeg_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#undef FALSE
#undef TRUE
#include <jpeglib.h>
#undef FALSE
#undef TRUE

#include <jsjpeg.h>

#include <bool.h>
#include <uint16_image.h>

uint16_image_t *uint16_image_read_jpeg_named(char *name, bool_t verbose, int32_t *spaceP);
  /* Reads an {uint16_image_t} image from the named JPEG file. 
  
    If {*spaceP} is not NULL, sets it to the image color space
    as in the {J_COLOR_SPACE} codes of {jpeglib.h}. 

    This procedure can read only 8-bits-per-pixel JPEG files.
    The parameter {BITS_IN_JSAMPLE} in the {jpeglib.h} file must be 8, and
    the JPEG library must have been compiled with that option.

    The parameters of the output image are described in
    {uint16_image_read_jpeg_INFO}. Note that For the color space
    {JCS_RGB565}, even though {img.maxval} is 255, the actual maximum
    samples in each channel are 31, 63, and 51. See
    {uint16_image_scale_channel} and {sample_conv_choose_maxval}.
    
    If the {name} is "-", reads from {stdin}. If {verbose} is TRUE,
    prints a notice to {stderr}. */
    
uint16_image_t *uint16_image_read_jpeg_file(FILE *rd, bool_t verbose, int32_t *spaceP);
  /* Same as {uint16_image_read_jpeg_named}, but from a previously opened file handle {rd}. */
    
#define uint16_image_read_jpeg_INFO \
  "The returned image will have {maxval = 255}, and may have 1 to 4 channels.  The" \
  " number, order, and meaning of the channels depends on the color space:\n" \
  "\n" \
  "    {JCS_GRAYSCALE}: 1 channel, intensity.\n" \
  "\n" \
  "    {JCS_RGB}, {JCS_EXT_RGB}, {JCS_RGB565}: 3 channels, red/green/blue.\n" \
  "\n" \
  "    {JCS_EXT_BGR}: 3 channels, blue/green/red.\n" \
  "\n" \
  "    {JCS_EXT_RGBX}: 4 channels, red/green/blue/junk.\n" \
  "\n" \
  "    {JCS_EXT_BGRX}: 4 channels, blue/green/red/junk.\n" \
  "\n" \
  "    {JCS_EXT_XBGR}: 4 channels, junk/blue/green/red.\n" \
  "\n" \
  "    {JCS_EXT_XRGB}: 4 channels, junk/red/green/blue.\n" \
  "\n" \
  "    {JCS_YCbCr}: 3 channels, Y/Cb/Cr (also known as YUV).\n" \
  "\n" \
  "    {JCS_CMYK}: 4 channels, cyan/magenta/yellow/black\n" \
  "\n" \
  "    {JCS_YCCK}: 4 channels, Y/Cb/Cr/black\n" \
  "\n" \
  "    {JCS_EXT_RGBA}: 4 channels, red/green/blue/alpha\n" \
  "\n" \
  "    {JCS_EXT_BGRA}: 4 channels, blue/green/red/alpha\n" \
  "\n" \
  "    {JCS_EXT_ABGR}: 4 channels, alpha/blue/green/red\n" \
  "\n" \
  "    {JCS_EXT_ARGB}: 4 channels, alpha/red/green/blue."

#endif
