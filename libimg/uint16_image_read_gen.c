/* See {uint16_image_read_gen.h} */
/* Last edited on 2024-12-20 17:18:51 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <argparser.h>

#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_read_jpeg.h>
#include <uint16_image_read_png.h>

#include <sample_conv.h>
#include <sample_conv_gamma.h>
#include <image_file_format.h>

#include <uint16_image_read_gen.h>

uint16_image_t *uint16_image_read_gen_named
  ( const char *fname, 
    image_file_format_t ffmt,
    uint32_t imaxval[], /* (OUT) Max sample value in each channel. */
    double *expoP,      /* (OUT) Gamma correction exponent specified/implied by input file. */
    double *biasP,      /* (OUT) Gamma correction bias parameter, idem. */
    bool_t verbose
  )
  {
    FILE *rd = open_read(fname, verbose);
    uint16_image_t *fimg = uint16_image_read_gen_file(rd, ffmt, imaxval, expoP, biasP, verbose);
    fclose(rd);
    return fimg;
  }
  
uint16_image_t *uint16_image_read_gen_frame
  ( const char *fpat, 
    int32_t fnum, 
    image_file_format_t ffmt,
    uint32_t imaxval[],  /* (OUT) Max sample value in each chanel. */
    double *expoP,       /* (OUT) Gamma conversion exponent specified/implied by input file. */
    double *biasP,       /* (OUT) Gamma conversion bias, idem. */
    bool_t verbose
  )
  {
    /* Insert the frame number in {fpat}: */
    char *fname = jsprintf(fpat, fnum);
    demand(strlen(fname) != 0, "duh? empty file name");
    /* Read the file: */
    uint16_image_t *fimg = uint16_image_read_gen_named(fname, ffmt, imaxval, expoP, biasP, verbose);
    free(fname);
    return fimg;
  }
  
uint16_image_t *uint16_image_read_gen_file
  ( FILE *rd,
    image_file_format_t ffmt,
    uint32_t imaxval[],  /* (OUT) Max sample value in each chanel. */
    double *expoP,       /* (OUT) Gamma conversion exponent specified/implied by input file. */
    double *biasP,       /* (OUT) Gamma conversion bias parameter, idem. */
    bool_t verbose
  )
  {   
    /* Read the input file as a {uint16_image_t}, without any sample conversion: */
    uint16_image_t *pimg; /* The quantized image. */
    double expo_file = NAN; /* Encoding expo specified by the file itself. */
    double bias_file = NAN; /* Encoding bias specified by the file itself. */    
    switch (ffmt)
      {
        case image_file_format_JPG:
          { int32_t space;
            pimg = uint16_image_read_jpeg_file(rd, verbose, &space);
            /* At present, {uint16_image_read_jpeg_named} can suport only two kinds. */
            /* Revise the next line if {uint16_image_read_jpeg_named} is expanded to support other kinds. */
            /* assert((space == JCS_GRAYSCALE) || (space == JCS_RGB)); */
            /* Use the client-given gamma correction parms for decoding: */
            expo_file = NAN;
            bias_file = NAN;
            /* The nominal {maxval} should be 255 for any {space}: */
            assert(pimg->maxval == 255);
            if (space == JCS_RGB565)
              { /* Packed RGB with 5,6,5 bits per sample: */
                assert(pimg->chns == 3);
                imaxval[0] = 31;
                imaxval[1] = 63;
                imaxval[2] = 31;
              }
            else
              { for (int32_t ic = 0; ic < pimg->chns; ic++)  { imaxval[ic] = 255; } }
          }
          break;

        case image_file_format_PNG:
          { bool_t verbose_png = TRUE;
            pimg = uint16_image_read_png_file (rd, &expo_file, imaxval, verbose_png);
            if (isnan(expo_file)) 
              { expo_file = sample_conv_gamma_sRGB_DEC_EXPO;
                bias_file = sample_conv_gamma_sRGB_BIAS;
              }
            else
              { expo_file = 1/expo_file; /* PNG "gAMA" gives the encoding gamma, not decoding. */
                bias_file = sample_conv_gamma_sRGB_BIAS; /* Guess, better than nothing. */
              }
          } 
          break;

        case image_file_format_PNM:
          { pimg = uint16_image_read_pnm_file(rd);
            /* Parameters that approximate the standard PNM ITU-R BT.709 decoding: */
            expo_file = sample_conv_gamma_BT709_DEC_EXPO; 
            bias_file = sample_conv_gamma_BT709_BIAS;
            for (int32_t ic = 0; ic < pimg->chns; ic++)  { imaxval[ic] = pimg->maxval; }
          }
          break;

        default:
          assert(FALSE);
      }
    (*expoP) = expo_file;
    (*biasP) = bias_file;
    return pimg;
  }

