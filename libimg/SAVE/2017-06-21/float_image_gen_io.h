/* I/O of {float_image_t} images in generic file formats. */
/* Last edited on 2017-06-16 20:20:40 by stolfilocal */

#ifndef float_image_gen_io_H
#define float_image_gen_io_H

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <float_image.h>
#include <argparser.h>

#include <float_image_file_format.h>

#define float_image_gen_read_MAX_CHANS (4)

/* IMAGE INPUT

  Each procedure in this section reads an image file in a specified
  format, including PNG, JPG, PNM (PBM, PGM, or PPM), and FNI; and converts 
  it to a common format {float_image_t}.
    
  For all formats except FNI, integer samples read from the file are
  converted from their original range {0..maxval} to the specified range
  {[vmin _ vmax]} by an affine map.  Depending on the image format,
  the {maxval} may be different for each channel.  The 
  estimated quantization error in each sample of each 
  channel {i} is returned in {qterr[i]}.  The array {qterr}
  must have at least {float_image_gen_read_MAX_CHANS} elements. 
  
  For FNI images, there is no conversion, and {qterr[i]} is set to zero.
  However, a warning is printed if the sample values are not 
  in the range {[vmin _ vmax]}.
  
  If gamma and bias decoding parameters are specified (explicitly or
  implicitly) in the input file, they are stored in {*gammaP} and
  {*biasP}; otherwise those variables are set to {NAN}. The pixels
  themselves are NOT gamma-decoded; the client may want to use
  {float_image_apply_gamma} on the result. */
    
#define float_image_gen_read_INFO \
  "The input file image format can be PNM (PBM, PGM, PPM), JPG," \
  " PNG, or FNI. In general, the file may be either RGB" \
  " color (three samples per pixel) or grayscale (one sample" \
  " per pixel), and may contain an alpha channel.\n" \
  "\n" \
  "  PNM (PBM, PGM, PPM) FORMAT:\n" \
  "\n" \
  "  The image will have" \
  " the same number of channels as the file, and samples" \
  " will be in the same order.\n" \
  "\n" \
  "  The parameters {*gammaP} and {*biasP} will be set" \
  " to values that provide a good approximation of the" \
  " default PNM encoding and decoding (BT709).  See {sample_conv.h} for details.\n" \
  "\n" \
  "  PNG FORMAT:\n" \
  "\n" \
  "  " uint16_image_io_png_read_INFO "\n" \
  "\n" \
  "  JPG (JPEG) FORMAT:\n" \
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
  "  Note that the channels are reordered, if needed, to" \
  " the more common red/green/blue/alpha order, and the" \
  " filler channels (\"X\") are discarded.\n" \
  "\n" \
  "  The samples of images with colorspace {JCS_RGB565} are" \
  " first mapped as integers from {0..31} and {0..63} to {0..255} and then to float values.  This introduces some rounding errors because {255/31} and {255/63} are not exact.   The maximum absolute error "

#define float_image_gen_io_size_MAX (float_image_max_size)
  /* Max expected image size along any axis. */

float_image_t *float_image_gen_read
  ( const char *fname,
    float_image_file_format_t ffmt,
    double vmin,      /* Output sample value corresponding to file sample value 0. */
    double vmax,      /* Output sample value corresponding to max file sample value. */
    double qterr[],   /* (OUT) Estimated quantization error in each channel. */
    double *gammaP,   /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,    /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose
  );
  /* Reads an image from file {fname}, assumed to be in format {ffmt},
    and converts it to a float image with pixel values in {[vmin _ vmax]}.

    If the {name} is "-", reads from {stdin}. Otherwise, the file name
    extension must be included in {fname}, but does not affect the
    assumed file format.

    If the parameter {verbose} is TRUE, the procedure prints some
    information messages to {stderr}. */

float_image_t *float_image_gen_fread
  ( FILE *rd,
    float_image_file_format_t ffmt,
    double vmin,      /* Output sample value corresponding to file sample value 0. */
    double vmax,      /* Output sample value corresponding to max file sample value. */
    double qterr[],   /* (OUT) Estimated quantization error in each channel. */
    double *gammaP,   /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,    /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose
  );
  /* Same as {float_image_gen_read}, but reads from the file handle {rd}
    instead of from a named file. The file {rd} must have been opened for reading,
    and is left open on return. */

float_image_t *float_image_gen_read_frame
  ( const char *fpat,
    int fnum,
    float_image_file_format_t ffmt,
    double vmin,      /* Output sample value corresponding to file sample value 0. */
    double vmax,      /* Output sample value corresponding to max file sample value. */
    double qterr[],   /* (OUT) Estimated quantization error in each channel. */
    double *gammaP,   /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,    /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose
  );
  /* Same as {float_image_gen_read}, but the file name is obtained by
    substituting the frame number {fnum} for the "%..d" specification in
    {fpat}. The file name extension must be included in
    {fpat} but does not affect the assumed file format. */

#endif
