/* See {float_image_write_gen.h} */
/* Last edited on 2020-11-06 00:14:07 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <argparser.h>

#include <uint16_image.h>
#include <uint16_image_write_jpeg.h>
#include <uint16_image_write_png.h>
#include <uint16_image_write_pnm.h>

#include <ix.h>
#include <sample_conv.h>
#include <float_image.h>
#include <float_image_to_uint16_image.h>
#include <image_file_format.h>

#include <float_image_write_gen.h>


void float_image_write_gen_frame
  ( const char *fpat,
    int fnum,
    float_image_t *fimg,
    image_file_format_t ffmt,
    double v0,        /* Output sample value corresponding to file sample value 0. */
    double vM,        /* Output sample value corresponding to max file sample value. */
    double gammaEnc,  /* Exponent parameter for gamma encoding. */
    double bias,      /* Bias parameter for gamma encoding. */
    bool_t verbose
  )
  {
    /* Insert the frame number in {fpat}: */
    char *fname = NULL;
    int nch = asprintf(&fname, fpat, fnum);
    demand(nch > 0, "invalid file name pattern");

    /* Write the file: */
    float_image_write_gen_named(fname, fimg, ffmt, v0, vM, gammaEnc, bias, verbose);

    free(fname);
  }

void float_image_write_gen_named
  ( const char *fname,
    float_image_t *fimg,
    image_file_format_t ffmt,
    double v0,        /* Output sample value corresponding to file sample value 0. */
    double vM,        /* Output sample value corresponding to max file sample value. */
    double gammaEnc,  /* Exponent parameter for gamma encoding. */
    double bias,      /* Bias parameter for gamma encoding. */
    bool_t verbose
  )
  { FILE *wr = open_write(fname, verbose);
    float_image_write_gen_file(wr, fimg, ffmt, v0, vM, gammaEnc, bias, verbose);
    fclose(wr);
  }

void float_image_write_gen_file
  ( FILE *wr,
    float_image_t *fimg,
    image_file_format_t ffmt,
    double v0,        /* Output sample value corresponding to file sample value 0. */
    double vM,        /* Output sample value corresponding to max file sample value. */
    double gammaEnc,  /* Exponent parameter for gamma encoding. */
    double bias,      /* Bias parameter for gamma encoding. */
    bool_t verbose
  )
  { 
    bool_t debug = FALSE;
    
    if (ffmt == image_file_format_FNI)
      { /* Write the image without any conversion: */
        float_image_write(wr, fimg); 
        return;
      }
 
    int32_t NC = (int32_t)fimg->sz[0]; /* Num channels. */
    
    /* Choose exponent and bias for gamma encoding: */
    if (isnan(gammaEnc)) { gammaEnc = sample_conv_BT709_ENC_GAMMA; }
    if (isnan(bias)) { bias = sample_conv_BT709_BIAS; }
    
    /* Copy {fimg} and scale so that the range fits snugly in {[-1_+1]}, then apply gamma: */
    double vd = fmax(fabs(v0), fabs(vM));
    float_image_t *gimg = float_image_copy(fimg);
    for (int c = 0; c < NC; c++) 
      { /* Rescale samples from {[-|vM| _ +|vM|]} to {[-1 _ +1]} or from {[0 _ vM]} to {[0 _ 1]}: */
        float_image_rescale_samples(gimg, c, 0, (float)vd, 0.0, 1.0);
        /* Apply gamma encoding: */
        if (gammaEnc != 1.0) { float_image_apply_gamma(gimg, c, gammaEnc, bias); }
      }

    /* Determine the range for the last encoding step, from {[eslo_eshi]} to {0..maxval}: */
    double eslo = sample_conv_gamma((float)(v0/vd), gammaEnc, bias);
    double eshi = sample_conv_gamma((float)(vM/vd), gammaEnc, bias);
    double slo[NC];      /* Low end of first mapping. */
    double shi[NC];      /* High end of first mapping. */
    for (int c = 0; c < NC; c++) { slo[c] = eslo; shi[c] = eshi; }
 
    /* Determine the max quantized sample value: */
    bool_t maxval = (ffmt == image_file_format_JPG ? 255 : 65535);
    
    /* Quantize samples of {gimg} to the {uint16_t} image {pimg}: */
    bool_t isMask  = FALSE;
    bool_t yup = FALSE;
    bool_t verbose_quant = FALSE;
    if (debug) { fprintf(stderr, "  quantizing image...\n"); }
    uint16_image_t *pimg = float_image_to_uint16_image(gimg, isMask, NC, slo, shi, NULL, maxval, yup, verbose_quant);
    float_image_free(gimg);
    
    /* Write {pimg} to file: */
    if (debug) { fprintf(stderr, "  writing PNG file...\n"); }
    switch (ffmt)
      {
        case image_file_format_JPG:
          { int32_t quality = 95;
            uint16_image_write_jpeg_file(wr, pimg, quality, verbose);
          }
          break;

        case image_file_format_PNG:
          { uint16_image_write_png_file (wr, pimg, gammaEnc, verbose);
          } 
          break;

        case image_file_format_PNM:
          { bool_t forceplain = FALSE;
            uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
          }
          break;

        default:
          demand(FALSE, "unimplemented image file format");
      }
    uint16_image_free(pimg);
  }
