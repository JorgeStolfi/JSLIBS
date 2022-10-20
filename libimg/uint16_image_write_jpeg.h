/* uint16_image_write_jpeg.h - routines to write images of {uint16_t} samples to JPEG files. */
/* Last edited on 2017-06-23 01:44:11 by stolfilocal */

/* Created by R. Minetto (IC-UNICAMP) sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#ifndef uint16_image_write_jpeg_H
#define uint16_image_write_jpeg_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#undef FALSE
#undef TRUE
#include <jpeglib.h>
#undef FALSE
#undef TRUE

#include <jsjpeg.h>

#include <bool.h>
#include <uint16_image.h>
  
void uint16_image_write_jpeg_named(char *name, uint16_image_t *img, int32_t quality, bool_t verbose);
  /* Writes image {img} to the specified file in JPEG format.
  
    If {img.chns} is 1 uses colorspace {JCS_GRAYSCALE}, if 3 assumes
    {JCS_RGB}, if 4 assumes {JCS_RGBA}. Other values are invalid.
    
    The {quality} parameter is an integer in {1..100} that specifies
    the fidelity of the compressed image to the original: 100 means no
    compression; typical values range in {80..98}.
    
    The routine can write only 8-bits-per-sample JPEG files.  The
    parameter {BITS_IN_JSAMPLE} in the {jpeglib.h} file must be 8, and
    the JPEG library must have been compiled with that option. The
    maximum output sample value {omaxval} will be 255. The maximum input
    image sample value {imaxval} can be arbitrary.
    
    If {imaxval} is equal to {omaxval} (255), the samples are sent to
    the JPEG encoder without any change. Otherwise samples of {img} are
    linearly scaled from the range {0..imaxval} to {0..omaxval} by the
    fractional factor {s = omaxval/imaxval}. If {imaxval>omaxval},
    information will be lost. If {imaxval} is less than {omaxval}, but is
    not a divisor of it, the conversion will be neither exact not
    uniform; namely, input samples that differ by 1 will be mapped to
    outout samples that differ by {floor(s)} or {ceil(s)}.
    
    If {name} is "-", writes to {stdout}. If {verbose} is TRUE, prints
    a notice to {stderr}. */

void uint16_image_write_jpeg_file(FILE *wr, uint16_image_t *img, int32_t quality, bool_t verbose);
  /* Same as {uint16_image_write_jpeg_named}, but to a previously opened file handle {wr}. */

#endif
