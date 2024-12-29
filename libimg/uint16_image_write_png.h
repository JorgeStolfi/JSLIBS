/* uint16_image_write_png.h - routines to write images with {uint16_t} samples as PNG files. */
/* Last edited on 2024-12-26 11:55:33 by stolfi */

/* Created by R. Minetto (IC-UNICAMP) as {ipng.h} sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#ifndef uint16_image_write_png_H
#define uint16_image_write_png_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <uint16_image.h>

void uint16_image_write_png_named(char *name, uint16_image_t *img, double gamma, bool_t verbose);
  /* Writes image {img} to the specified file in PNG format.  The 
    image may have 1 to 4 channels which are interpteted as specified in 
    {uint16_image_write_png_INFO}.
    
    The output number of bits per sample {ubits} is set according to the
    image's max sample value {imaxval = img.maxval}, as explained in
    {uint16_image_write_png_SAMPLE_CONV_INFO}. The image tries to avoid
    or mininmize rounding errors as much as possible.
  
    If {name} is "-", writes to {stdout}. If {verbose} is TRUE, prints
    a notice to {stderr}. */

void uint16_image_write_png_file(FILE *wr, uint16_image_t *img, double gamma, bool_t verbose);
  /* Same as {uint16_image_read_png_named} and {uint16_image_write_png_named}, 
    but to a previously opened file handle {wr}. */
    
#define uint16_image_write_png_INFO \
  "When writing to a PNG file, the image may have 1 to 4 channels," \
  " which are interpreted as follows:\n" \
  "\n" \
  "   1: PNG type will be grayscale image.\n" \
  "\n" \
  "   2: PNG type will be grayscale + alpha.\n" \
  "\n" \
  "   3: PNG type will be RGB.\n" \
  "\n" \
  "   4: PNG type will be RGB + alpha.\n" \
  "\n" \
  "   The file's bits per sample {fbits} will be set to one of" \
  " the PNG valid values -- namely, 1, 2, 4, 8, or 16 -- so" \
  " that the the image's max sample value {imaxval} does not" \
  " exceed {2^fbits-1}.  Note that the file is never color-mapped, so" \
  " the {fbits} values 1, 2, and 4 are valid only if the image has" \
  " a single channel.\n" \
  "\n" \
  "   The file may specify a true sample bit size {ubits} less" \
  " than {fbits}, via the 'sBIT' chunk.  Depending on {ubits} and {imaxval}, the" \
  " samples may also have to be linearly rescaled so that {imaxval} is" \
  " mapped to {omaxval=2^ubits-1}." \
  " If {imaxval} is not a divisor of {omaxval}, the conversion may" \
  " slightly distort the sample values.\n" \
  "\n" \
  "   The {gamma} parameter is written to the file as the 'gAMA' chunk." \
  " The samples are not modified.  If {gamma} is NAN, zero, or negative," \
  " the chunk is omitted."

#define uint16_image_write_png_SAMPLE_CONV_INFO \
  "The true bit size {ubits} of samples in the file is chosen" \
  " depending on the input image's max sample value" \
  " {imaxval = img.maxval} so as to avoid or minimize rounding errors" \
  " in sample conversion.\n" \
  "\n" \
  "  Namely, if {imaxval} is a divisor of {2^k-1}, for some {k} between" \
  " 1 and 16, then {ubits} will be the smallest such {k}.  If {imaxval} is" \
  " not a divisor of {2^k-1} for any {k} up to 16, then {ubits} will be" \
  " set to 16.\n" \
  "\n" \
  "  Every image sample {ismp} in {0..imaxval} will then be scaled to" \
  " a file sample {osmp} in {0..omaxval}, where {omaxval} is {2^ubits-1}, by the" \
  " scaling factor {s=omaxval/imaxval}.\n" \
  "\n" \
  "  If the fraction {s} is an" \
  " integer (that is, {omaxval} is a multiple of {imaxval}), the" \
  " conversion will not entail any rounding and will preserve the original" \
  " brightness values exactly.\n" \
  "\n" \
  "  If the fraction {s} is not an integer, the mapping from image samples" \
  " to file samples will have 'hiccups' (uneven steps) due to" \
  " rounding.  Input sample values" \
  " that differ by 1 will map to samples that differ by either {floor(s)} or {ceil(s)}.\n" \
  "\n" \
  "  In any case, input sample values {0} and {imaxval} will" \
  " always map to {0} and {omaxval}, respectively.\n" \
  "\n" \
  "  In any case, the file's nominal bit size {fbits} will be the" \
  " smallest valid value that is not less than {ubits}.  If {fbits} is" \
  " different than {ubits}, an 'sBIT' chunk will be written to" \
  " specify {ubits} as the true bits size.\n" \
  "\n" \
  "  More precisely, scaling will be unnecessary or exact for the" \
  " following values of {imaxval}:" \
  "\n" \
  "    1: {ubits=1,omaxval=1}\n" \
  "\n" \
  "    3: {ubits=2,omaxval=3}\n" \
  "\n" \
  "    7: {ubits=3,omaxval=7}\n" \
  "\n" \
  "    5,15: {ubits=4,omaxval=15}\n" \
  "\n" \
  "    31: {ubits=5,omaxval=31}\n" \
  "\n" \
  "    9,21,63: {ubits=6,omaxval=63}\n" \
  "\n" \
  "    127: {ubits=7,omaxval=127}\n" \
  "\n" \
  "    17,51,85,255: {ubits=8,omaxval=255}\n" \
  "\n" \
  "    73,511: {ubits=9,omaxval=511}\n" \
  "\n" \
  "    11,33,93,341,1023: {ubits=10,omaxval=1023}\n" \
  "\n" \
  "    23,89,2047: {ubits=11,omaxval=2047}\n" \
  "\n" \
  "    13,35,39,45,65,91,105,117,195,273,315,455,585,\n" \
  "    819,1365,4095: {ubits=12,omaxval=4095}\n" \
  "\n" \
  "    8191: {ubits=13,omaxval=8191}\n" \
  "\n" \
  "    43,129,381,5461,16383: {ubits=14,omaxval=16383}\n" \
  "\n" \
  "    31,151,217,1057,4681,32767: {ubits=15,omaxval=32767}\n" \
  "\n" \
  "    257,771,1285,3855,4639,13107,21845,65535: {ubits=16,omaxval=65535}."

#endif
