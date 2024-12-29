/* uint16_image_read_png.h - routines to read PNG image files as images with {uint16_t} samples. */
/* Last edited on 2024-12-26 11:54:16 by stolfi */

/* Created by R. Minetto (IC-UNICAMP) as {ipng.h} sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#ifndef uint16_image_read_png_H
#define uint16_image_read_png_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <uint16_image.h>

#define uint16_image_read_png_MAX_CHNS (4)
  /* Max channels in a PNG image that can be read with this interface. */

uint16_image_t *uint16_image_read_png_named(char *name, double *gammaP, uint32_t imaxval[], bool_t verbose);
  /* Reads an image from the named PNG file.  The  resulting {uint16_image_t}
    may have 1 to 4 channels, as described in {uint16_image_read_png_named_INFO}. 
    
    Also stores in {*gammaP} the {gamma} exponent to be used when
    interpreting the samples. If the file has an "sRGB" chunk, the
    {*gammaP} is set set to {NAN}, to indicate the standard sRGB
    intensity encoding/decoding. Otherwise, if the the file includes a
    "gAMA" chunk, it sets {*gammaP} to the value stored in that chunk.
    Otherwise sets {*gammaP} to {NAN}.
    
    The array {imaxval} must have at least 4 elements.  The parameter {gammaP} is
    set, but no gamma conversion is applied to the image. 
    
    If the {name} is \"-\", reads from {stdin}. If {verbose} is TRUE," \
    prints a notice to {stderr}. */

uint16_image_t *uint16_image_read_png_file(FILE *rd, double *gammaP, uint32_t imaxval[], bool_t verbose);
  /* Same as {uint16_image_read_png_named}, but from a previously opened file handle {rd}. */

#define uint16_image_read_png_named_INFO \
  "The image read from a PNG file may have 1 to 4 channels:\n" \
  "\n" \
  "    1: PNG type was grayscale image.\n" \
  "\n" \
  "    2: PNG type was grayscale + alpha.\n" \
  "\n" \
  "    3: PNG type was RGB.\n" \
  "\n" \
  "    4: PNG type was RGB + alpha.\n" \
  "\n" \
  "  If an alpha channel is present, it will be the last one.\n" \
  "\n" \
  "   Stores in {imaxval[i]} the original max sample value for each" \
  " channel {i}, namely {2^{buse[i]}-1}" \
  " where {buse[i]} is the number of bits actually used for the samples" \
  " in that channel.  This number is obtained from the \"sBIT\" chunk in" \
  " the file; or, if that chunk is not present, from the number" \
  " of bits per sample of te input file.  The" \
  " possible values of {buse[i]} are 1 to 16, implying" \
  " {imaxval[i]} 1 to 65535, respectively.\n" \
  "\n" \
  "   Also stores in {*gammaP} the" \
  " {gamma} exponent specified in the \"gAma\" chunk of the file, if" \
  " present; otherwise sets {*gammaP} to NAN."

#endif
