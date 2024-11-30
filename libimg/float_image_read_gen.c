/* See {float_image_read_gen.h} */
/* Last edited on 2020-11-06 00:27:22 by jstolfi */

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
#include <uint16_image_read_pnm.h>
#include <uint16_image_read_jpeg.h>
#include <uint16_image_read_png.h>

#include <ix.h>
#include <sample_conv.h>
#include <float_image.h>
#include <float_image_from_uint16_image.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <image_file_format.h>

#include <float_image_read_gen.h>

double float_image_read_gen_pick_parameter(const char *p_name, double p_given, double p_file, double p_default, bool_t verbose);
  /* Chooses the gamma or bias to be used for decoding an input file.
    The result is the client-given value {p_given} if not {NAN},
    else the value {p_file} specified implicitly or implicitly in the 
    input file if not {NAN}, else the default value {p_default}.
    Prints the choice if verbose, using the string {p_name} as the 
    parameter's name.  */

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
  )
  {
    /* Insert the frame number in {fpat}: */
    int nch = char *fname = jsprintf(fpat, fnum);
    demand(nch > 0, "invalid file name pattern");
    /* Read the file: */
    float_image_t *fimg = float_image_read_gen_named(fname, ffmt, v0, vM, maxvalP, gammaDecP, biasP, verbose);
    free(fname);
    return fimg;
  }
  
float_image_t *float_image_read_gen_named
  ( const char *fname,
    image_file_format_t ffmt,
    float v0,           /* Output sample value corresponding to file sample value 0. */
    float vM,           /* Output sample value corresponding to max file sample value. */
    uint16_t **maxvalP, /* (OUT) Discrete nominal max value in file for each channel. */
    double *gammaDecP,  /* (OUT) Exponent of gamma-decoding specified or implied in the input file. */
    double *biasP,      /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose      /* If true, prints some information about the file and conversion. */ 
  )
  {
    FILE *rd = open_read(fname, verbose);
    float_image_t *fimg = float_image_read_gen_file(rd, ffmt, v0, vM, maxvalP, gammaDecP, biasP, verbose);
    fclose(rd);
    return fimg;
  }
  
float_image_t *float_image_read_gen_file
  ( FILE *rd,
    image_file_format_t ffmt,
    float v0,           /* Output sample value corresponding to file sample value 0. */
    float vM,           /* Output sample value corresponding to max file sample value. */
    uint16_t **maxvalP, /* (OUT) Discrete nominal max value in file for each channel. */
    double *gammaDecP,  /* (OUT) Gamma specified or implied in the input file. */
    double *biasP,      /* (OUT) Bias parameter for gamma conversion, idem. */
    bool_t verbose      /* If true, prints some information about the file and conversion. */ 
  )
  {   
    if (ffmt == image_file_format_FNI)
      { /* Read the image without any conversion: */
        float_image_t *fim = float_image_read(rd); 
        /* Check parameter ranges: */
        if (gammaDecP != NULL) { (*gammaDecP) = NAN; }
        if (biasP != NULL) { (*biasP) = NAN; }
        if (maxvalP != NULL) { (*maxvalP) = NULL; }
        return fim;
      }
    
    /* Read the input file as a {uint16_image_t}, without any sample conversion: */
    uint16_image_t *pimg = NULL; /* The quantized image, read below. */
    double gammaEnc_file = NAN; /* Encoding gamma specified by the file itself. */
    double bias_file = NAN; /* Encoding bias specified by the file itself. */    
    int32_t NC_MAX = 4; /* Max number of channels expected in image file, except FNI. */
    uint32_t maxval_file[NC_MAX]; 
    switch (ffmt)
      {
        case image_file_format_JPG:
          { int32_t space; /* JPEG color space code. */
            uint16_image_t *pimg = uint16_image_read_jpeg_file(rd, verbose, &space);
            /* At present, {uint16_image_read_jpeg_named} can suport only two kinds. */
            /* Revise the next line if {uint16_image_read_jpeg_named} is expanded to support other kinds. */
            assert((space == JCS_GRAYSCALE) || (space == JCS_RGB));
            /* Use the client-given gamma for decoding: */
            for (uint32_t c = 0;  c < pimg->chns; c++) { maxval_file[c] = 255; }
            gammaEnc_file = NAN;
            bias_file = NAN;
          }
          break;

        case image_file_format_PNG:
          { bool_t verbose_png = TRUE;
            pimg = uint16_image_read_png_file (rd, &gammaEnc_file, maxval_file, verbose_png);
            bias_file = NAN;
          } 
          break;

        case image_file_format_PNM:
          { pimg = uint16_image_read_pnm_file(rd);
            /* Parameters that approximate the standard PNM ITU-R BT.709 decoding: */
            for (uint32_t c = 0;  c < pimg->chns; c++) { maxval_file[c] = pimg->maxval; }
            gammaEnc_file = sample_conv_BT709_ENC_GAMMA; 
            bias_file = sample_conv_BT709_BIAS;
          }
          break;

        default:
          demand(FALSE, "unimplemented image file format");
      }

    /* Get and check image channel count: */
    int32_t NC = (int32_t)pimg->chns; /* Num channels. */
    if (verbose) { fprintf(stderr, "width = %d height  = %d  channels = %d\n", pimg->cols, pimg->rows, NC); }
    demand((NC >= 1) && (NC <= 4), "invalid channel count in frame");

    /* Decide which exponent and bias should be used for gamma decoding: */
    double gammaEnc = float_image_read_gen_pick_parameter("gamma", NAN, gammaEnc_file, sample_conv_BT709_ENC_GAMMA, verbose);
    double bias =  float_image_read_gen_pick_parameter("bias",  NAN, bias_file,  sample_conv_BT709_BIAS, verbose);

    /* Determine the range for the first decoding step, from {0..maxval} to {[eslo_eshi]}: */
    double vd = fmax(fabs(v0), fabs(vM));
    double eslo = sample_conv_gamma((float)(v0/vd), gammaEnc, bias);
    double eshi = sample_conv_gamma((float)(vM/vd), gammaEnc, bias);
    double slo[NC];      /* Low end of first mapping. */
    double shi[NC];      /* High end of first mapping. */
    for (int c = 0; c < NC; c++) { slo[c] = eslo; shi[c] = eshi; }
    
    /* Convert to float image in the range {[eslo_eshi]}: */
    bool_t isMask = FALSE;          /* TRUE for masks, FALSE for images. */
    bool_t yup = FALSE;             /* TRUE to reverse the indexing of rows. */
    bool_t verbose_float = verbose; /* TRUE to debug the conversion. */
    float_image_t *fimg = float_image_from_uint16_image(pimg, isMask, slo, shi, yup, verbose_float);
    
    /* Discards the pixel array and header of the PNM image. */
    uint16_image_free(pimg);

    /* Aplly gamma correction: */
    for (int c = 0; c < NC; c++) 
       { /* Apply gamma decoding: */
         if (gammaEnc != 1.0) { float_image_apply_gamma(fimg, c, 1/gammaEnc, bias); }
         /* Rescale samples from {[v0/vd _ vM/vd]} to {[v0 _ vM]}: */
         float_image_rescale_samples(fimg, c, 0.0f, 1.0f, 0.0f, (float)vd);
       }

    if (gammaDecP != NULL) { (*gammaDecP) = 1/gammaEnc; }
    if (biasP != NULL) { (*biasP) = bias; }

    return fimg;
  }

double float_image_read_gen_pick_parameter(const char *p_name, double p_given, double p_file, double p_default, bool_t verbose)
  {
    demand(! isnan(p_default), "invalid default");
    double p;
    if (! isnan(p_given))
      { p = p_given;
        if (verbose && (! isnan(p_file)) && (p_given != p_file)) 
          { fprintf(stderr, "file %s (= %.4f) overriden\n", p_name, p_file); }
      }
    else if (!isnan(p_file))
      { p = p_file;
        if (verbose) 
          { fprintf(stderr, "using file %s = %.4f\n", p_name, p_file); }
      }
    else
      { p = p_default;
        if (verbose) 
          { fprintf(stderr, "using default %s = %.4f\n", p_name, p_default); }
      }
    demand(p > 0.0, "invalid parameter value");
    return p;
  }

