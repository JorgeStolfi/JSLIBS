/* See {uint16_image_write_png.h} */
/* Last edited on 2017-09-14 16:48:40 by jstolfi */

/* Created by R. Minetto (IC-UNICAMP) as {ipng.c} sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include <png.h>
#ifndef _SETJMP_H
  #include <setjmp.h>
#endif

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jspng.h>
#include <uint16_image.h>
#include <sample_conv.h>

#include <uint16_image_write_png.h>

void uint16_image_write_png_named(char *name, uint16_image_t *img, double gamma, bool_t verbose)
  {
    FILE *wr = open_write (name, verbose);
    uint16_image_write_png_file(wr, img, gamma, verbose);
    if ((wr != stdout) && (wr != stderr)) { fclose(wr); }
  }

void uint16_image_write_png_file(FILE *wr, uint16_image_t *img, double gamma, bool_t verbose)
  {
    bool_t debug = FALSE;
    bool_t debug_conv = FALSE;
    
    png_uint_32 chns = img->chns;
    png_uint_32 cols = img->cols;
    png_uint_32 rows = img->rows;
    
    uint16_t imaxval = img->maxval;  /* Image's {maxval}. */
    demand(imaxval + 0 <= 65535, "maxval is too big"); /* Paranoia, since {uint16_t} is 16 bits. */
    demand(imaxval > 0, "maxval is zero");

    /* Create {libpng} work records {pw, pi}. */
    /* Avoid {setjmp/longjmp} complexity for now, use default error handling. */
    png_struct *pw = notnull(png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL), "no mem");
    png_info *pi = notnull(png_create_info_struct(pw), "no mem");

    /* Assign {wr} as the output file for writes to {pw}: */
    png_init_io(pw, wr);

    /* Set the {zlib} compression level: */
    png_set_compression_level(pw, 7); /* Level 6 is default, level 9 is maximum. */

    /* Determine a suitable bit size {img_smp_bits} for the input samples: */
    /* Pick for the smallest {k} such that {imaxval} divides {2^k-1}; or 16 if none: */
    uint8_t img_smp_bits = 16;  /* . */
    uint16_t omaxval = 65535; 
    for (int32_t k = 1; k < 16; k++) 
      { uint32_t kval = (1u << k) - 1u;
        if ((kval % imaxval) == 0)
          { img_smp_bits = (uint8_t)k; 
            omaxval = (uint16_t)kval; 
            break;
          }
      }
    assert((omaxval >= imaxval) && (omaxval + 0 <= 65535));
    
    /* Select output image type and nominal bit size: */
    png_byte ctype;
    png_byte png_smp_bits;
    if (chns == 1)
      { ctype = PNG_COLOR_TYPE_GRAY;
        /* Bit size can be 1, 2, 4, 8, or 16: */
        png_smp_bits = (omaxval == 1 ? 1 : (omaxval <= 3 ? 2 : (omaxval <= 15 ? 4: (omaxval <= 255 ? 8 : 16))));
      }
    else
      { if (chns == 2)
          { ctype = PNG_COLOR_TYPE_GRAY_ALPHA; }
        else if (chns == 3)
          { ctype = PNG_COLOR_TYPE_RGB; }
        else if (chns == 4)
          { ctype = PNG_COLOR_TYPE_RGB_ALPHA; }
        else
          { demand(FALSE, "invalid number of channels"); }
        /* Bit size can be only 8 or 16: */
        png_smp_bits = (omaxval <= 255 ? 8 : 16);
      }
    assert(img_smp_bits <= png_smp_bits);
    
    if (img_smp_bits != png_smp_bits)
      { /* Must write an "sBit" chunk: */
        png_color_8 sBIT;
        sBIT.gray = img_smp_bits;
        sBIT.red = img_smp_bits;
        sBIT.green = img_smp_bits;
        sBIT.blue = img_smp_bits;
        sBIT.alpha = img_smp_bits;
        png_set_sBIT(pw, pi, &sBIT);
      }
    
    if (verbose) 
      { if (imaxval != omaxval) 
          { fprintf(stderr, "%s: ", __FUNCTION__);
            fprintf(stderr, "rescaling samples from {0..%d} to {0..%d}\n", imaxval, omaxval);
            if ((omaxval % imaxval) != 0) { fprintf(stderr, "!! warning: scaling will not be exact\n"); }
          }
      }
    
    png_set_IHDR
      ( pw, pi, cols, rows, png_smp_bits, ctype, 
        PNG_INTERLACE_NONE, 
        PNG_COMPRESSION_TYPE_DEFAULT, 
        PNG_FILTER_TYPE_DEFAULT
      );
    if (debug) { jspng_dump_info(stderr, __FUNCTION__, "before write_info", pw, pi); }

    /* Provide a "gAMA" chunk. */
    if (! isnan(gamma) && (gamma > 0)) 
      { png_set_gAMA(pw, pi, gamma);
        if (debug) { jspng_dump_info(stderr, __FUNCTION__, "after set_gAMA", pw, pi); }
      }

    /* Write the PNG header: */
    png_write_info (pw, pi);
    if (debug) { jspng_dump_info(stderr, __FUNCTION__, "after write_info", pw, pi); }
    fflush(wr);

    /* -- OUTPUT TRANSFORMATION REQUESTS -------------------------------- */

    /* Ask the PNG library to squeeze input samples into 1,2,4-bit images if they fit: */
    if (png_smp_bits < 8) { png_set_packing(pw); } 

    /* ------------------------------------------------------------------ */
    
    /* Decide how many bytes per sample are needed for {png_write_row}. */
    /* Note that samples with less than 8 bits must be given unpacked, one per byte. */
    assert(sizeof(png_byte) == 1);
    assert((png_smp_bits <= 8) || (png_smp_bits == 16));
    uint32_t bytes_per_png_sample = (png_smp_bits <= 8 ? 1 : 2); /* Bytes per pixel for {png_write_row} */
    
    /* Compute how many bytes are needed to store one row of PNG file samples: */
    /* Cannot use {png_get_rowbytes()} pecause it does not account for {png_set_packing} (CROCK!). */
    uint32_t samples_per_row = (uint32_t)(cols*chns); /* Samples per row in PNG file and image. */
    uint32_t bytes_per_png_row = samples_per_row * bytes_per_png_sample; /* Bytes per row for {png_write_row} */
    if (debug) { fprintf(stderr, "bytes_per_png_row = %u\n", bytes_per_png_row); }
        
    /* Allocate buffer for one row of PNG samples: */
    png_byte *png_rowP = notnull(malloc(bytes_per_png_row + 1), "no mem");

    /* Process rows: */
    assert((imaxval + 0 <= 65535) && (omaxval + 0 <= 65535)); /* So that there is no overflow in the conversion. */
    uint32_t iround = imaxval / 2u;  /* Bias for sample conversion rounding. */
    assert(iround < imaxval);
    for (int32_t row = 0; row < rows; ++row)
      {
        if (debug_conv) { fprintf(stderr, "    [ "); }
        png_bytep png_P = png_rowP;
        uint16_t *img_P = img->smp[row];
        int32_t k;
        for (k = 0; k < samples_per_row; k++)
          { uint32_t ismp = (*img_P); img_P++;
            /* Apply scale factor: */
            uint32_t osmp;
            if (imaxval == omaxval) 
              { osmp = ismp; }
            else
              { osmp = (uint32_t)((omaxval*ismp + iround)/imaxval); }
            /* Store in PNG row buffer: */
            if (debug_conv) { fprintf(stderr, " %u->%u", ismp, osmp); }
            if (png_smp_bits <= 8)
              { /* {png_write_row} wants one byte per sample: */
                (*png_P) = (png_byte)osmp; png_P++;
              }
            else
              { /* {png_write_row} wants two bytes per sample: */
                (*png_P) = (png_byte)(osmp >> 8); png_P++;
                (*png_P) = (png_byte)(osmp & 255); png_P++;
              }
          }
        png_write_row(pw, png_rowP);
        if (debug_conv) { fprintf(stderr, " ]\n"); }
      }    

    png_write_end(pw, pi);

    /* png_destroy_write_struct(&pw, &pi); */
    png_free_data(pw, pi, PNG_FREE_ALL, -1);
    free(png_rowP);
  }

