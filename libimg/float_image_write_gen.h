/* Write a {float_image_t} to a file with variable format. */
/* Last edited on 2024-12-20 17:41:21 by stolfi */

#ifndef float_image_write_gen_H
#define float_image_write_gen_H

#include <stdio.h>

#include <bool.h>
#include <float_image.h>
#include <argparser.h>

#include <image_file_format.h>

#define float_image_write_gen_MAX_CHANS (4)

/* Each procedure in this section writes given {float_image_t} in-memory image
  to an image file in a specified
  format, including PNG, JPG, PNM (PBM, PGM, or PPM), and FNI; converting it 
  appropriately.
    
  For FNI images, there is no conversion; the samples are written
  out in a suitabl decimal floating-point format.
  
  For all other formats except FNI, {float} samples from the image
  are converted from their original assumed range {[v0 _ vM]} 
  to the range {0..maxval}, where {maxval} depends on the image file 
  format. 
  
  Image samples that are outside the range {[v0_vM]} are clipped to that
  range, with a warning. However, a warning is printed if the sample
  values are not in the range {[v0 _ vM]}.
  
  The scaled samples are then gamma-encoded with exponent {expoEnc} and
  bias {bias} before quantization. See {sample_conv_gamma}. If
  {expoEnc} is {NAN}, uses a the ITU-R BT709 encoding exponent
  {sample_conv_gamma_BT709_ENC_EXPO}. If {bias} is {NAN}, uses the
  BT709 value {sample_conv_gamma_BT709_BIAS}. To suppress gamma
  encoding, specify {expoEnc=1}. */
    
#define float_image_write_gen_INFO \
  "The output file image format can be PNM (PBM, PGM, PPM), JPG," \
  " PNG, or FNI.  The number {NC} of channels in the image to be written" \
  " must be between 1 and 4.  The specific variant of the specified file" \
  " format will depend on {NC}.  It may be either RGB color (three samples" \
  " per pixel) or grayscale (one sample" \
  " per pixel), and may contain an alpha channel, depending on the format.  The" \
  " samples will be quantized to integer in some range {0..maxval}, where" \
  " {maxval} depends on the the file format:\n" \
  "\n" \
  "  PNM (PBM, PGM, PPM) FORMAT:\n" \
  "\n" \
  "  The input image must not have an alpha channel.  The procedure will write" \
  " a raw PGM file if {C} is 1, and a raw PPM file if {NC} is 3.  The max sample" \
  " value will be 65535 (16 bits per sample).\n" \
  "\n" \
  "  The given {expoEnc} exponent will be used in the quantization but cannot be" \
  " saved in the file.  PNM files assume the ITU-R BT709 gamma" \
  " encoding.  See {sample_conv.h} for details.\n" \
  "\n" \
  "  PNG FORMAT:\n" \
  "\n" \
  "  The file will contain a true color (not colormapped) image in" \
  " the {GRAY}, {GRAY+ALPHA}, {RGB}, or {RGB+ALPHA} format if {NC} is 1, 2, 3 or 4," \
  " respectively. The max quantized sample value {maxval} in the file will be 65535 (16 bits" \
  " per sample).  The {expoEnc} exponent will be written as a \"gAMA\" chunk.\n" \
  "\n" \
  "  JPG (JPEG) FORMAT:\n" \
  "\n" \
  "  The image may have only 1, 3 or 4 channels, and the file's colorspace will" \
  " be {JCS_GRAYSCALE}, {JCS_RGB}, or {JCS_RGBA}, respectively.  The max quantized" \
  " sample value will be 255 (8 bits per sample)."

void float_image_write_gen_named
  ( const char *fname,         /* File name. */
    float_image_t *fimg,       /* Image to write out. */
    image_file_format_t ffmt,  /* File format. */
    double v0,                 /* Output sample value corresponding to file sample value 0. */
    double vM,                 /* Output sample value corresponding to max file sample value. */
    double expoEnc,           /* Exponent parameter for gamma encoding. */
    double bias,               /* Bias parameter for gamma encoding. */
    bool_t verbose
  );
  /* Writes the image {fimg} to file {fname} in format {ffmt},
    mapping {[v0 _ vM]} (with gamma encoding) to {0..maxval}
    as described in {float_image_write_gen_INFO}.

    If the {name} is "-", writes to {stdout}. Otherwise, the file name
    extension must be included in {fname}, but does not affect the
    assumed file format.

    If the parameter {verbose} is TRUE, the procedure prints some
    information messages to {stderr}. */

void float_image_write_gen_file
  ( FILE *wr,
    float_image_t *fimg,
    image_file_format_t ffmt,
    double v0,      /* Output sample value corresponding to file sample value 0. */
    double vM,      /* Output sample value corresponding to max file sample value. */
    double expoEnc,  /* Exponent parameter for gamma encoding. */
    double bias,      /* Bias parameter for gamma encoding. */
    bool_t verbose
  );
  /* Same as {float_image_write_gen_named}, but reads from the file handle {rd}
    instead of from a named file. The file {wr} must have been opened for writing,
    and is left open on return. */

void float_image_write_gen_frame
  ( const char *fpat,
    int32_t fnum,
    float_image_t *fimg,
    image_file_format_t ffmt,
    double v0,       /* Output sample value corresponding to file sample value 0. */
    double vM,       /* Output sample value corresponding to max file sample value. */
    double expoEnc,  /* Exponent parameter for gamma encoding. */
    double bias,     /* Bias parameter for gamma encoding. */
    bool_t verbose
  );
  /* Same as {float_image_write_gen_named}, but the file name is obtained by
    substituting the frame number {fnum} for the "%..d" specification in
    {fpat}. The file name extension must be included in
    {fpat} but does not affect the assumed file format. */

#endif
