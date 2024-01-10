/* See {dsm_image.h} */
/* Last edited on 2017-06-20 20:46:45 by stolfilocal */

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
#include <uint16_image_io_jpeg.h>
#include <uint16_image_io_png.h>

#include <indexing.h>
#include <sample_conv.h>
#include <float_image.h>
#include <float_image_from_to_uint16_image.h>
#include <frgb.h>
#include <frgb_ops.h>

#include <dsm_image_format.h>
#include <dsm_gamma.h>

#include <dsm_image.h>

float_image_t *dsm_image_read_frame
  ( const char *fpat, 
    int fnum, 
    dsm_image_format_t ffmt,
    double gamma,
    double vmax,
    bool_t gray
  )
  {
    /* Insert the frame number in {fpat}: */
    char *fname = NULL;
    int nch = asprintf(&fname, fpat, fnum);
    demand(nch > 0, "invalid file name pattern");
    /* Read the file: */
    float_image_t *fimg = dsm_image_read(fname, ffmt, gamma, vmax, gray);

    free(fname);
    return fimg;
  }
  
float_image_t *dsm_image_read
  ( const char *fname, 
    dsm_image_format_t ffmt,
    double gamma,
    double vmax,
    bool_t gray
  )
  {
    bool_t verbose = TRUE;
    FILE *rd = open_read(fname, verbose);
    float_image_t *fimg = dsm_image_fread(rd, ffmt, gamma, vmax, gray);
    fclose(rd);
    return fimg;
  }
  
float_image_t *dsm_image_fread
  ( FILE *rd,
    dsm_image_format_t ffmt,
    double gamma,
    double vmax,
    bool_t gray
  )
  {   
    bool_t verbose = TRUE;
    
    /* Read the input file as a {uint16_image_t}, without any sample conversion: */
    uint16_image_t *pimg; /* The quantized image. */
    double gamma_file = NAN; /* Encoding gamma specified by the file itself. */
    double bias_file = NAN; /* Encoding bias specified by the file itself. */    
    switch (ffmt)
      {
        case dsm_image_format_JPG:
          { int kind;
            pimg = uint16_image_io_jpeg_fread(rd, &kind);
            /* At present, {uint16_image_io_jpeg_read} can suport only two kinds. */
            /* Revise the next line if {uint16_image_io_jpeg_read} is expanded to support other kinds. */
            assert((kind == JCS_GRAYSCALE) || (kind == JCS_RGB));
            /* Use the client-given gamma for decoding: */
            gamma_file = NAN;
            bias_file = NAN;
          }
          break;

        case dsm_image_format_PNG:
          { bool_t verbose_png = TRUE;
            pimg = uint16_image_io_png_fread (rd, &gamma_file, verbose_png);
            gamma_file = 1/gamma_file; /* PNG "gAMA" gives the encoding gamma, not decoding. */
            bias_file = NAN;
          } 
          break;

        case dsm_image_format_PNM:
          { pimg = uint16_image_read_pnm_file(rd);
            /* Parameters that approximate the standard PNM ITU-R BT.709 decoding: */
            gamma_file = 1.0/0.450; 
            bias_file = 0.0327;
          }
          break;

        default:
          assert(FALSE);
      }

    /* Get and check image channel count: */
    int32_t NC = (int32_t)pimg->chns; /* Num channels. */
    if (verbose) { fprintf(stderr, "width = %d height  = %d  channels = %d\n", pimg->cols, pimg->rows, NC); }
    demand((NC == 1) || (NC == 3), "invalid channel count in frame");

    /* Determine the range for the first decoding step, from {0..maxval} to {[eslo_eshi]}: */
    double eslo = (vmax < 0 ? -1.0 : 0.0);
    double eshi = 1.0;
    double slo[NC];      /* Low end of first mapping, (0 or {-1}). */
    double shi[NC];      /* High end of first mapping ({+1}). */
    for (int c = 0; c < NC; c++) { slo[c] = eslo; shi[c] = eshi; }
    
    /* Convert to float image in the range {[-1 _ +1]} or {[0 _ 1]}: */
    bool_t isMask = FALSE;        /* TRUE for masks, FALSE for images. */
    bool_t yup = FALSE;           /* TRUE to reverse the indexing of rows. */
    bool_t verbose_float = FALSE; /* TRUE to debug the conversion. */
    float_image_t *fimg = float_image_from_uint16_image(pimg, isMask, slo, shi, yup, verbose_float);
    
    /* Discards the pixel array and header of the PNM image. */
    uint16_image_free(pimg);

    /* Decide which gamma and bias should be used for decoding: */
    double bias = NAN;
    gamma = dsm_gamma_pick_parameter("gamma", gamma, gamma_file, 1.0, verbose);
    bias =  dsm_gamma_pick_parameter("bias",  bias,  bias_file,  sample_conv_DEFAULT_BIAS, verbose);

    /* Aplly gamma correction: */
    for (int c = 0; c < NC; c++) 
       { /* Apply gamma correction: */
         if (gamma != 1.0) { float_image_apply_gamma(fimg, c, gamma, bias); }
         /* Rescale samples from {[-1 _ +1]} to {[vmax _ -vmax]} or from {[0 _ 1]} to {[0 _ vmax]}: */
         float_image_rescale_samples(fimg, c, 0, 1, 0, (float)fabs(vmax));
       }

    /* This step must be done after gamma correction: */
    if ((NC > 1) && gray)
      { /* Convert to grayscale: */
        assert(NC == 3);
        float_image_make_grayscale(fimg);
        /* Reduce to single channel: */
        fimg->sz[0] = 1; 
        if (verbose) { fprintf(stderr, "converted to grayscale: channels = %d\n", (int32_t)fimg->sz[0]); }
      }
      
    return fimg;
  }

void dsm_image_write_frame(const char *fpat, int fnum, float_image_t *fimg, double gamma, double bias)
  {
    /* Insert the frame number in {fpat}: */
    char *fname = NULL;
    int nch = asprintf(&fname, fpat, fnum);
    demand(nch > 0, "invalid file name pattern");

    /* Write the file: */
    dsm_image_write(fname, fimg, gamma, bias);

    free(fname);
  }

void dsm_image_write
  ( const char *fname, 
    float_image_t *fimg,
    double gamma,
    double vmax
  )
  { bool_t verbose = TRUE;
    FILE *wr = open_write(fname, verbose);
    dsm_image_fwrite(wr, fimg, gamma, vmax);
    fclose(wr);
  }

void dsm_image_fwrite
  ( FILE *wr, 
    float_image_t *fimg,
    double gamma,
    double vmax
  )
  { 
    bool_t debug = FALSE;
    
    int32_t NC = (int32_t)fimg->sz[0]; /* Num channels. */
    
    /* Copy {fimg} and scale so that the range fits snugly in {[-1_+1]}, then apply gamma: */
    float_image_t *gimg = float_image_copy(fimg);
    if (isnan(gamma)) { gamma = 1.0; }
    double bias = sample_conv_DEFAULT_BIAS;;
    for (int c = 0; c < NC; c++) 
      { /* Rescale samples from {[-|vmax| _ +|vmax|]} to {[-1 _ +1]} or from {[0 _ vmax]} to {[0 _ 1]}: */
        float_image_rescale_samples(gimg, c, 0, (float)fabs(vmax), 0.0, 1.0);
        /* apply inverse gamma correction: */
        if (gamma != 1.0) { float_image_apply_gamma(gimg, c, 1/gamma, bias); }
      }

    /* Determine the last scaling step, from {[0 _ 1]} or {[-1 _ +1]} to {0..maxval}: */
    double eslo = (vmax < 0 ? -1.0 : 0.0);
    double eshi = 1.0;
    double slo[NC];      /* Low end of first mapping, (0 or {-1}). */
    double shi[NC];      /* High end of first mapping ({+1}). */
    for (int c = 0; c < NC; c++) { slo[c] = eslo; shi[c] = eshi; }

    /* Quantize {gimg} to the PNM image {pimg}, with the specified gamma: */
    bool_t isMask  = FALSE;
    bool_t yup = FALSE;
    bool_t verbose_pnm = FALSE;
    bool_t maxval = 65535; /* For a 16-bit-per-sample PNM image. */
    if (debug) { fprintf(stderr, "  quantizing image...\n"); }
    uint16_image_t *pimg = float_image_to_uint16_image(gimg, isMask, NC, slo, shi, NULL, maxval, yup, verbose_pnm);
    float_image_free(gimg);
    
    /* Write {pimg} as a PNG file with gamma 1.0 (since we did gamma ourselves): */
    if (debug) { fprintf(stderr, "  writing PNG file...\n"); }
    bool_t verbose_png = FALSE;
    uint16_image_io_png_fwrite (wr, pimg, 1.0, verbose_png);
    uint16_image_free(pimg);
  }
