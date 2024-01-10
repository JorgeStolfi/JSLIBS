/* uint16_image_io_jpeg.h - routines to read and write JPEG files. */
/* Last edited on 2017-06-20 20:50:45 by stolfilocal */

/* Created by R. Minetto (IC-UNICAMP) sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#ifndef uint16_image_io_jpeg_H
#define uint16_image_io_jpeg_H

#define _GNU_SOURCE
#include <stdio.h>

#undef FALSE
#undef TRUE
#include <jpeglib.h>
#undef FALSE
#undef TRUE

#include <bool.h>
#include <uint16_image.h>

/* These procedures read/write a JPEG file into/from the {uint16_image_t} 
  memory format defined in {uint16_image.h}.
  
  These procedures can read and write only 8-bits-per-pixel JPEG files.
  The parameter {BITS_IN_JSAMPLE} in the {jpeglib.h} file must be 8, and
  the JPEG library must have been compiled with that option.

  In theory, some JPEG files may use 12 bits per sample. However, in
  order to read those files it would be necessary to recompile the
  library with {BITS_IN_JSAMPLE=12}; and then it would not be able to
  read 8-bit JPEG files. (WHAT A CROCK!). Fortunately no one seems to
  generate 12-bit JPEGs, so we will not care about that case.
  */
  
#define jsjpeg_BITS_PER_SAMPLE (BITS_IN_JSAMPLE)
  /* The number of bits per sample in the JPEG file as defined by {jpeglib.h} (8 or 12). */

#define jsjpeg_MAXVAL (MAXJSAMPLE)
  /* The max sample value as as defined by {jpeglib.h} (255 or 4095). */

uint16_image_t *uint16_image_io_jpeg_read (char *name, bool_t verbose, int32_t *spaceP);
  /* Reads an image from the named JPEG file. 
  
    If {*spaceP} is not NULL, sets it to the image color space
    as in the {J_COLOR_SPACE} codes of {jpeglib.h}. 
    
    The parameters of the output image are described in {uint16_image_io_jpeg_read_INFO}.
    Note that For the color space {JCS_RGB565}, even though {img.maxval} is 255, 
    the actual maximum samples in each channel are 31, 63, and 51. 
    See {uint16_image_scale_channel} and {sample_conv_choose_maxval}.
    
    If the {name} is "-", reads from {stdin}. If {verbose} is TRUE,
    prints a notice to {stderr}. */
    
#define uint16_image_io_jpeg_read_INFO \
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
    
uint16_image_t *uint16_image_io_jpeg_fread (FILE *rd, bool_t verbose, int32_t *spaceP);
  /* Same as {uint16_image_io_jpeg_read}, but from a previously opened file handle {rd}. */

void uint16_image_io_jpeg_write (char *name, uint16_image_t *img, int32_t quality, bool_t verbose);
  /* Writes image {img} to the specified file in JPEG format.
  
    If {img.chns} is 1 uses colorspace {JCS_GRAYSCALE}, if 3 assumes
    {JCS_RGB}, if 4 assumes {JCS_RGBA}. Other values are invalid.
    
    The {quality} parameter is an integer in {1..100} that specifies
    the fidelity of the compressed image to the original: 100 means no
    compression; typical values range in {80..98}.
    
    The routine can write only 8-bits-per-sample JPEG files, so the
    maximum output sample value {omaxval} is 255. The maximum input
    image sample value {imaxval} can be arbitrary.
    
    If {imaxval} is equal to {omaxval} (255), the samples are sent to
    the JPEG encoder without any change. Otherwise samples of {img} are
    linearly scaled from the range {0..imaxval} to {0..omaxval} by the
    fractional factor {s = omaxval/imaxval}. If {imaxval>omaxval},
    information will be lost. If {imaxval} is less than {omaxval}, bt is
    not a divisor of it, the conversion will be neither exact not
    uniform; namely, input samples that differ by 1 will be mapped to
    outout samples that differ by {floor(s)} or {ceil(s)}.
    
    If {name} is "-", writes to {stdout}. If {verbose} is TRUE, prints
    a notice to {stderr}. */

void uint16_image_io_jpeg_fwrite (FILE *wr, uint16_image_t *img, int32_t quality, bool_t verbose);
  /* Same as {uint16_image_io_jpeg_write}, but to a previously opened file handle {wr}. */

#endif
