/* Reading images with {uint16_t} samples from image files in variable formats. */
/* Last edited on 2017-06-29 18:29:06 by stolfilocal */

#ifndef uint16_image_read_gen_H
#define uint16_image_read_gen_H

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <uint16_image.h>
#include <argparser.h>

#include <image_file_format.h>

#define uint16_image_read_gen_MAX_CHANS (4)
  /* Max number of channels in images produced by this interface. */

uint16_image_t *uint16_image_read_gen_named
  ( const char *fname,
    image_file_format_t ffmt,
    uint32_t imaxval[], /* (OUT) Max sample value in each channel. */
    double *gammaP,     /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,      /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose
  );
  /* Reads an image from file {fname}, assumed to be in format {ffmt};
    which includes PNG, JPEG, and PNM (PBM, PGM, or PPM).
    
    If the {name} is "-", reads from {stdin}. Otherwise, the file name
    extension must be included in {fname}, but does not affect the
    assumed file format.

    The contents of the file is restructured into an
    in-memory image {img}, in the common format {uint16_image_t}.
    See {uint16_image_read_gen_INFO} below.
    
    Some image files have different sample ranges in each channel.
    The input procedures below will save in {imaxval[i]} the
    actual max sample value in channel {i} of the file.  The array {imaxval}
    must have at least {uint16_image_read_gen_MAX_CHANS} elements. 

    The max sample value of the returned image {img.maxval} will be at
    least as large as all of the channel {imaxval}s.  Note that the samples
    are not scaled to match this max value; the samples of channel {i}
    remains in the ranges {0..imaxval[i]}.

    If gamma and bias decoding parameters are specified (explicitly or
    implicitly) in the input file, they are stored in {*gammaP} and
    {*biasP}; otherwise those variables are set to {NAN}. The pixels
    themselves are NOT gamma-decoded. 

    If the parameter {verbose} is TRUE, the procedure prints some
    information messages to {stderr}. */

uint16_image_t *uint16_image_read_gen_file
  ( FILE *rd,
    image_file_format_t ffmt,
    uint32_t *imaxval[], /* (OUT) Max sample value in each chanel. */
    double *gammaP,      /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,       /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose
  );
  /* Same as {uint16_image_read_gen_named}, but reads from the file handle {rd}
    instead of from a named file. The file {rd} must have been opened for reading,
    and is left open on return. */

uint16_image_t *uint16_image_read_gen_frame
  ( const char *fpat,
    int fnum,
    image_file_format_t ffmt,
    uint32_t *imaxval[], /* (OUT) Max sample value in each chanel. */
    double *gammaP,      /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,       /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose
  );
  /* Same as {uint16_image_read_gen_named}, but the file name is obtained by
    substituting the frame number {fnum} for the "%..d" specification in
    {fpat}. The file name extension must be included in
    {fpat} but does not affect the assumed file format. */

#define uint16_image_read_gen_INFO \
  "The input file image format can be PNM (NetPBM: PBM, PGM, PPM), JPEG," \
  " or PNG (Portable Network Graphics). In general, the file may be either RGB" \
  " color (three samples per pixel) or grayscale (one sample" \
  " per pixel), and may contain an extra channel for transparency (alpha).  Some files may contain data in other color spaces, such as the subtractive cyan/yellow/magenta or the CIE Y/Cb/Cr (both possibly with an extra black channel).\n" \
  "\n" \
  "  Reading from PNM (PBM, PGM, PPM) files:\n" \
  "\n" \
  "  The image will have" \
  " the same number of channels as the file, and samples" \
  " will be in the same order.\n" \
  "\n" \
  "  The parameters {*gammaP} and {*biasP} will be set" \
  " to values that provide a good approximation of the" \
  " default PNM encoding and decoding (BT709).  See {sample_conv.h} for details.\n" \
  "\n" \
  "  Reading from PNG files:\n" \
  "\n" \
  "  " uint16_image_read_png_INFO "\n" \
  "\n" \
  "  Reading from JPG (JPEG) files:\n" \
  "\n" \
  "  Only 8-bit-per-pixel JPEG images are supported.  The returned" \
  " image may have 1, 3, or 4 channels, depending on the JPEG color space" \
  " specified in the file; which also determines their meaning.  Namely:\n" \
  "\n" \
  "    JCS_GRAYSCALE: 1 channel, intensity.\n" \
  "\n" \
  "    JCS_RGB, JCS_EXT_RGB, JCS_RGB565, JCS_EXT_BGR,\n" \
  "    JCS_EXT_RGBX, JCS_EXT_BGRX, JCS_EXT_XBGR,\n" \
  "    JCS_EXT_XRGB: 3 channels, red/green/blue.\n" \
  "\n" \
  "    JCS_YCbCr: 3 channels, Y/Cb/Cr (also known as YUV).\n" \
  "\n" \
  "    JCS_CMYK: 4 channels, cyan/magenta/yellow/black.\n" \
  "\n" \
  "    JCS_YCCK: 4 channels, Y/Cb/Cr/black.\n" \
  "\n" \
  "    JCS_EXT_RGBA, JCS_EXT_BGRA, JCS_EXT_ABGR,\n" \
  "    JCS_EXT_ARGB: 4 channels, red/green/blue/alpha.\n" \
  "\n" \
  "  Note that the three primary color channels are reordered, if needed, to" \
  " the more common red/green/blue/alpha order, and any" \
  " filler channels (\"X\") are discarded."
  
#define uint16_image_read_gen_MAXVAL_INFO \
  "The red, green, and blue samples of images with colorspace {JCS_RGB565} will" \
  " range in {0..31}, {0..63}, and {0..31}, respectively, and the returned image" \
  " will have {img.maxval} set to 63.  For all other colorspaces, the returned image" \
  " will have {img.maxval} 255 and the nominal range of samples in all channels will be {0..255}."

#endif
