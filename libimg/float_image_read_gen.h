/* Reafing {float_image_t} images from variable file formats. */
/* Last edited on 2017-09-02 17:11:15 by stolfilocal */

#ifndef float_image_read_gen_H
#define float_image_read_gen_H

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <float_image.h>
#include <argparser.h>

#include <image_file_format.h>

float_image_t *float_image_read_gen_named
  ( const char *fname,
    image_file_format_t ffmt,
    float v0,            /* Output sample value corresponding to file sample value 0. */
    float vM,            /* Output sample value corresponding to max file sample value. */
    uint16_t **maxvalP,  /* (OUT) Discrete nominal max value in file for each channel. */
    double *gammaDecP,   /* (OUT) Exponent of gamma-decoding specified or implied in the input file. */
    double *biasP,       /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose       /* If true, prints some information about the file and conversion. */ 
  );
  /* Reads an image from file {fname}, and converts it to a float-valued
    image {fim}.
    
    The input file is assumed to be in format {ffmt}, including PNG,
    JPG, PNM (PBM, PGM, or PPM), and FNI.

    If the {name} is "-", the procedure reads from {stdin}. Otherwise, the
    file name extension must be included in {fname}, but does not affect
    the assumed file format.

    If the parameter {verbose} is TRUE, the procedure prints some
    information messages to {stderr}.

    Images read from an FNI file may have any number of channels, with
    arbitrary meanings. For other formats, the number of channels {NC}
    in the returned image may be 1, 2, 3, or 4. Usually, the image is
    grauscale if {NC} is 1 or 2, and RGB color if {NC} is 3 or 4. When
    {NC} is 2 or 4, the last channel is opacity ("alpha"). However, JPEG
    files may use other 3- or 4-channel color spaces, such as CYMK or
    Y/Cr/Cb.
    
    Sample values read from an FNI file are returned without any change.
    The parameters {v0} and {vM} are ignored. For all other file formats
    except FNI, integer samples read from the file are converted to
    float values as described next, and detailed in
    {float_image_read_gen_INFO} and {float_image_read_gen_CONV_INFO}
    below.
    
    Depending on the image format, the range of integer samples stored
    in the file may be different for each channel. If {maxvalP} is not
    {NULL}, The procedure may set {*maxvalP} to the address of a newly
    allocated array {maxval[0..NC-1]}, such that {maxval[c]} is the max
    integer value allowed in that channel. Note that it may not be the
    actual maximum integer sample read from the file. For FNI images,
    however, if {maxvalP} is not {NULL}, the procedure sets {*maxvalP}
    to {NULL}.
    
    Samples in image files are usually related to light intensity values
    by a non-linear /gamma encoding/ that depends on two parameters, an
    exponent {gammaDec} and an offset {bias}. If the file specifies or
    implies values for these parameters, they are returned in
    {*gammaDecP} and {*biasP}. If the file does not specify them, the
    procedure assumes the parameters of the ITU-R BT709 standard; see
    {sample_conv.h}. For FNI images, {*gammaDecP} and {*biasP} are set
    to {NAN}. */
    
float_image_t *float_image_read_gen_file
  ( FILE *rd,
    image_file_format_t ffmt,
    float v0,           /* Output sample value corresponding to file sample value 0. */
    float vM,           /* Output sample value corresponding to max file sample value. */
    uint16_t **maxvalP, /* (OUT) Discrete nominal max value in file for each channel. */
    double *gammaDecP,  /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,      /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose      /* If true, prints some information about the file and conversion. */ 
  );
  /* Same as {float_image_read_gen_named}, but reads from the file handle {rd}
    instead of from a named file. The file {rd} must have been opened for reading,
    and is left open on return. */

float_image_t *float_image_read_gen_frame
  ( const char *fpat,
    int fnum,
    image_file_format_t ffmt,
    float v0,           /* Output sample value corresponding to file sample value 0. */
    float vM,           /* Output sample value corresponding to max file sample value. */
    uint16_t **maxvalP, /* (OUT) Discrete nominal max value in file for each channel. */
    double *gammaDecP,  /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,      /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose      /* If true, prints some information about the file and conversion. */ 
  );
  /* Same as {float_image_read_gen_named}, but the file name is obtained by
    substituting the frame number {fnum} for the "%..d" specification in
    {fpat}. The file name extension must be included in
    {fpat} but does not affect the assumed file format. */
    
#define float_image_read_gen_INFO \
  "The input file image format can be PNM (PBM, PGM, PPM), JPG," \
  " PNG, or FNI. In general, the file may be either RGB" \
  " color (three samples per pixel) or grayscale (one sample" \
  " per pixel), and may contain an alpha channel.\n" \
  "\n" \
  "  FNI FORMAT:\n" \
  "\n" \
  "  The FNI image file format corresponds is a faithful and" \
  " complete dump of the {float_image_t} memory" \
  " representation.  Therefore, the image is read just" \
  " as it is in the file, without any conversion," \
  " scaling, or gamma correction.\n" \
  "\n" \
  "  PNM (PBM, PGM, PPM) FORMAT:\n" \
  "\n" \
  "  The image will have" \
  " the same number of channels as the file, and samples" \
  " will be in the same order.\n" \
  "\n" \
  "  The samples will be gamma-corrected using exponent and offset" \
  " parameters that provide a good approximation of the" \
  " default PNM decoding (ITU-R BT709).  See" \
  " {sample_conv.h} for details.\n" \
  "\n" \
  "  PNG FORMAT:\n" \
  "\n" \
  "  " float_image_read_png_INFO "\n" \
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
  " first mapped as integers from {0..31} and {0..63} to {0..255}" \
  " and then to float values.  This introduces some rounding errors" \
  " because {255/31} and {255/63} are not exact."
  
#define float_image_read_gen_CONV_INFO \
  "For all input formats except FNI, the integer samples read from" \
  " the file are converted to {float} values as follows.\n" \
  "\n" \
  "  First, the integer samples are mapped affinely from their original" \
  " discrete range {0..M} to some sub-interval of " \
  " {[-1 _ 1]}," \
  " depending of two given parameters {v0} and {vM}.  Specifically, the smalled discrete" \
  " sample value 0 is mapped to {E(v0/vd)} and the " \
  " maxmimum sample value {M} is mapped to {E(vM/vd)}, where" \
  " {vd = max{|v0|,|vM|}} and {E} is the \"gamma\" encoding function.   The {E} function uses the exponent" \
  " and offset specified or implied by the file, as per above.  See" \
  " {sample_conv_gamma} for details.\n" \
  "\n" \
  "  Then a \"gamma\" decoding mapping {v = D(v)} is applied to every" \
  " sample, using the same exponent and bias.  As a resut of these two" \
  " steps, the original sample values {0} and {M} are mapped to {v0/vd}" \
  " and {vM/vd}, respectively.\n" \
  "\n" \
  "  Then the sample values are multiplied by the scaling factor {vd} above.\n" \
  "\n" \
  "  Note that the final float samples will be in the range {v0} to {vM} only" \
  " of there is no gamma correction ({*gammaDecP} is set to 1) or if {{v0,vM}} is" \
  " a subset of {{-1,0,+1}}.  The two parameters must be distinct, but {vM] may" \
  " be less than {v0} to invert the brightness scale."

#endif
