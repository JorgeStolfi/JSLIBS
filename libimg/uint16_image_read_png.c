/* See {uint16_image_read_png.h} */
/* Last edited on 2019-02-06 18:21:14 by jstolfi */

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

#include <uint16_image_read_png.h>

uint16_image_t *uint16_image_read_png_named(char *name, double *gammaP, sample_uint32_t imaxval[], bool_t verbose)
  {
    FILE *rd = open_read (name, verbose);
    uint16_image_t *img = uint16_image_read_png_file(rd, gammaP, imaxval, verbose);
    if (rd != stdin) { fclose(rd); }
    return img;
  }

uint16_image_t *uint16_image_read_png_file(FILE *rd, double *gammaP, uint32_t imaxval[], bool_t verbose)
  {
    bool_t debug = TRUE;
    
    /* Check the PNG signature: */
    unsigned char sig[jspng_MAGIC_BYTES];
    fread(sig, 1, jspng_MAGIC_BYTES, rd);
    demand(png_sig_cmp(sig, 0, jspng_MAGIC_BYTES) == 0, "not a PNG file");

    /* Allocate the reader structures: */
    png_structp pr = notnull(png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL), "no mem");
    png_infop pi = notnull(png_create_info_struct(pr), "no mem");
    png_infop pe = notnull(png_create_info_struct(pr), "no mem");

    /* Tell {libpng} to abotr in case of errors: */
    if (setjmp(png_jmpbuf(pr)))
      { demand(FALSE, "libpng error");
        /* png_destroy_read_struct(&pr, &pi, &pe); */
        /* fclose(rd); */
        return NULL;
      }
    
    /* Attach the file handle to the PNG reader, tell it we already read the sig: */
    png_init_io(pr, rd);
    png_set_sig_bytes(pr, jspng_MAGIC_BYTES);
    
    /* Read the header: */
    png_read_info(pr, pi);
    if (debug) { jspng_dump_info(stderr, __FUNCTION__, "before transforms", pr, pi); }

    /* Get image dimensions: */
    png_uint_32 cols = png_get_image_width(pr, pi);
    demand((cols > 0) && (cols <= jspng_MAX_SIZE), "invalid num of cols"); 
    
    png_uint_32 rows = png_get_image_height(pr, pi);
    demand((rows > 0) && (rows <= jspng_MAX_SIZE), "invalid num of rows"); 
     
    /* Get channel count before transformations: */
    png_uint_32 orig_chns = png_get_channels(pr, pi);
    demand((orig_chns >= 1) && (orig_chns <= 4), "invalid unmapped channel count"); 

    /* Get number of bits per pixel before all transformations: */
    uint8_t orig_bits = png_get_bit_depth(pr, pi); /* Bits per sample in file. */
    demand((orig_bits >= 1) && (orig_bits <= 16), "invalid unmapped bits per sample in file");

    /* Get color type before transformations: */
    uint8_t orig_color_type = png_get_color_type(pr, pi);
    
    /* Bits per sample in actual color data (not in palette index) before unpacking: */
    uint8_t default_smp_bits = (orig_color_type == PNG_COLOR_TYPE_PALETTE ? 8 : orig_bits);
    
    /* -- INPUT TRANSFORMATION REQUESTS -------------------------------- */

    /* Request expansion of 1,2,4-bit samples to one byte: */
    /* This will leave each sample in the LOWER {bits} bits of each byte. */
    if (orig_bits < 8) 
      { if (debug) { fprintf(stderr, "requesting unpacking of %u-bit samples\n", orig_bits); }
        png_set_packing(pr);
      } 

    /* Request conversion of palette colors to RGB or RGBA: */
    if (orig_color_type == PNG_COLOR_TYPE_PALETTE) 
      { if (debug) { fprintf(stderr, "requesting expansion of %u-bit palette indices to true RGB\n", orig_bits); }
        png_set_palette_to_rgb(pr); 
      }

    /* If palette with a transparent color, request conversion to RGBA: */
    if (png_get_valid(pr, pi, PNG_INFO_tRNS)) 
      { demand(orig_color_type != PNG_COLOR_TYPE_GRAY_ALPHA, "tRNS chunk in GRAY+ALPHA file");
        demand(orig_color_type != PNG_COLOR_TYPE_RGB_ALPHA, "tRNS chunk in RGB+ALPHA file");
        if (debug) { fprintf(stderr, "requesting conversion of 'tRNS' chunk to alpha channel\n"); }
        png_set_tRNS_to_alpha(pr);
        /* If there is "tRNS" chunk, it seems that small GRAY/RGB pixels are scaled to 8-bit: */
        if (default_smp_bits < 8) { default_smp_bits = 8; }
      }
    else
      { /* No tRNS chunk; default alpha bit size is same as other channels: */
        default_smp_bits = default_smp_bits;
      }
    
    /* If there is an 'sBIT' chunk, request that the original bit sizes per channel be preserved: */
    png_color_8p sBIT;
    bool_t has_sBIT = png_get_sBIT(pr, pi, &sBIT);
    if (has_sBIT) 
      { if (debug) 
          { fprintf(stderr, "requesting samples at true bit sizes");
            fprintf(stderr, " Y%d R%d G%d B%d A%d\n", sBIT->gray, sBIT->red, sBIT->green, sBIT->blue, sBIT->alpha);
          }
        png_set_shift(pr, sBIT);
      }

    /* Request de-interlacing: */
    (void)png_set_interlace_handling(pr);

    /* ------------------------------------------------------------------ */

    /* Update the derived info fields: */
    png_read_update_info(pr, pi);
    if (debug) { jspng_dump_info(stderr, __FUNCTION__, "after transform requests", pr, pi); }
     
    /* Get channel count after transformations: */
    png_uint_32 chns = png_get_channels(pr, pi);
    demand((chns >= 1) && (chns <= jspng_MAX_CHANS), "invalid channel count"); 
    
    /* Get number of bits per pixel after transformations, ensure that unpacking worked: */
    uint8_t png_smp_bits = png_get_bit_depth(pr, pi); /* Bits per sample in transformed image. */
    demand((png_smp_bits == 8) || (png_smp_bits == 16), "packed sample expansion ({png_set_packing}) failed");
    
    /* Get color type after transformations, ensure that palette was expanded: */
    uint8_t color_type = png_get_color_type(pr, pi);
    if (orig_color_type == PNG_COLOR_TYPE_PALETTE)
      { demand
          ( (color_type == PNG_COLOR_TYPE_RGB) || (color_type == PNG_COLOR_TYPE_RGB_ALPHA),
            "palette expansion ({png_set_palette_to_rgb}) yielded wrong color type"
          );
        demand(png_smp_bits == 8, "invalid sample size for palette colors");
      }
    
    /* Get the gamma exponent of the file (but do not apply gamma correction): */
    
    (*gammaP) = NAN;
    if (png_get_valid(pr, pi, PNG_INFO_sRGB))
      { /* The sRGB color space has its own gamma. Return {NAN} to signify it: */
        if (verbose) { fprintf(stderr, "uses sRGB intensity encoding\n"); }
        (*gammaP) = NAN;
      }
    else if (png_get_gAMA(pr, pi, gammaP)) 
      { if (verbose) { fprintf(stderr, "file gamma = %.8f\n", (*gammaP)); } }
    else
      { if (verbose) { fprintf(stderr, "file does not specify a gamma value\n"); }
        (*gammaP) = NAN;
      }
    
    /* Consistency of channel count and color space, and set {buse[0..chns-1]}: */
    uint8_t buse[chns]; /* Used bits per channel. */
    char *ctype;
    if (color_type == PNG_COLOR_TYPE_GRAY)
      { demand(chns == 1, "inconsistent channel count for GRAY");
        buse[0] = (has_sBIT ? sBIT->gray : default_smp_bits);
        ctype = "GRAY";
      }
    else if (color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
      { demand(chns == 2, "inconsistent channel count GRAY+ALPHA");
        buse[0] = (has_sBIT ? sBIT->gray : default_smp_bits);
        buse[1] = (has_sBIT ? sBIT->alpha : default_smp_bits);
        ctype = "GRAY+ALPHA";
      }
    else if (color_type == PNG_COLOR_TYPE_RGB)
      { demand(chns == 3, "inconsistent channel count for RGB");
        buse[0] = (has_sBIT ? sBIT->red : default_smp_bits);
        buse[1] = (has_sBIT ? sBIT->green : default_smp_bits);
        buse[2] = (has_sBIT ? sBIT->blue : default_smp_bits);
        ctype = "RGB";
      }
    else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA)
      { demand(chns == 4, "inconsistent channel count for RGB+ALPHA");
        buse[0] = (has_sBIT ? sBIT->red : default_smp_bits);
        buse[1] = (has_sBIT ? sBIT->green : default_smp_bits);
        buse[2] = (has_sBIT ? sBIT->blue : default_smp_bits);
        buse[3] = (has_sBIT ? sBIT->alpha : default_smp_bits);
        ctype = "RGB+ALPHA";
      }
    else if (color_type == PNG_COLOR_TYPE_PALETTE)
      { demand(FALSE, "lbpng failed to expand the palette"); }
    else 
      { demand(FALSE, "invalid color type"); }

    if (verbose) 
      { fprintf(stderr, "PNG file: color type = %s  chns = %u", ctype, chns);
        fprintf(stderr, "  bits per channel =");
        for (int chn = 0; chn < chns; chn++) { fprintf(stderr, " %u", buse[chn]); }
        fprintf(stderr, "  size = %u x %u\n", cols, rows);
      }
    
    /* Validate used bit counts, compute {imaxval[0..chns]} and choose output {omaxval}: */
    sample_uint32_t omaxval = 0;   /* Output maxval. */
    for (int chn = 0; chn < chns; chn++)
      { demand ((buse[chn] > 0) && (buse[chn] <= png_smp_bits), "invalid bit count in sBIT chunk"); 
        imaxval[chn] = (1u << buse[chn]) - 1u;
        demand(imaxval[chn] <= PNM_MAX_SAMPLE, "file sample size too big");
        if (imaxval[chn] > omaxval) { omaxval = imaxval[chn]; }
      }

    /* Allocate the image in memory: */
    uint16_image_t *img = uint16_image_new((int32_t)cols, (int32_t)rows, (int32_t)chns);
    img->maxval = (uint16_t)omaxval;
  
    if (verbose) 
      { fprintf(stderr, "image in memory: chns = %u", chns);
        fprintf(stderr, "  maxval = %u", img->maxval);
        fprintf(stderr, "  size = %u x %u\n", cols, rows);
      }
    
    /* Compute how many bytes per sample are needed for {png_read_row}. */
    /* Note that samples with less than 8 bits will be unpacked, one per byte. */
    assert(sizeof(png_byte) == 1);
    assert((png_smp_bits == 8) || (png_smp_bits == 16));
    uint32_t bytes_per_png_sample = png_smp_bits/8;  /* Bytes per sample after transformations. */
    
    /* Compute how many bytes are needed to store one row of PNG file samples: */
    uint32_t samples_per_row = (uint32_t)(cols*chns); /* Samples per row in PNG file and image. */
    uint32_t bytes_per_png_row = samples_per_row * bytes_per_png_sample; /* Bytes per row for {png_write_row} */
    if (debug) 
      { fprintf(stderr, "bytes_per_png_row = %u\n", bytes_per_png_row);
        fprintf(stderr, "png_get_rowbytes() = %u\n", (uint32_t)png_get_rowbytes(pr, pi));
      }
    assert(bytes_per_png_row == (uint32_t)png_get_rowbytes(pr, pi));

    /* Allocate buffer for the PNG image data: */
    png_bytep *png_dataP = notnull(malloc(rows*sizeof(png_bytep)), "no mem");
    for (uint32_t row = 0;  row < rows; row++)
      { png_dataP[row] = notnull(malloc(bytes_per_png_row), "no mem"); }

    /* Read the image samples: */
    if (debug) { fprintf(stderr, "loading image ...\n"); }
    png_read_image(pr, png_dataP);
    if (debug) { fprintf(stderr, "image loaded\n"); }

    /* Read and convert the pixels, row by row: */
    assert((png_smp_bits == 8) || (png_smp_bits == 16));
    for (uint32_t row = 0;  row < rows; row++)
      { 
        /* Process one row of the image: */
        png_bytep png_P = png_dataP[row];
        uint16_t *img_P = img->smp[row];
        for (uint32_t col = 0;  col < cols; col++)
          { for (uint32_t chn = 0;  chn < chns; chn++)
              {
                /* Get sample {smp} from PNG row buffer: */
                uint32_t smp = (*png_P); png_P++;
                if (png_smp_bits == 16)
                  { /* PNG file has two bytes per sample: */
                    smp = (smp << 8) | (*png_P); png_P++;
                  }
                
                /* Sample range consistency: */
                if (smp > imaxval[chn]) 
                  { fprintf(stderr, "invalid sample");
                    fprintf(stderr, "  col = %d  row = %d  chn = %d", col, row, chn);
                    fprintf(stderr, "  smp = %u max = %u", smp, imaxval[chn]);
                    /* Try to get the upper and lower {buse[chn]} bits: */
                    uint8_t bsz = buse[chn];
                    fprintf(stderr, "  smp(lo) = %u", smp & ((1u << bsz) - 1u));
                    fprintf(stderr, "  smp(hi) = %u", smp >> (png_smp_bits - bsz));
                    /* Try to downscale {smp} from {2^png_smp_bits-1} to {imaxval[chn]}: */
                    uint32_t full_maxval = (1u << png_smp_bits) - 1u;
                    fprintf(stderr, "  smp(sc) = %u", (smp * imaxval[chn])/full_maxval);
                    fprintf(stderr, "\n");
                  }
                assert(smp <= imaxval[chn]);
                
                /* Store into {img}: */
                (*img_P) = (uint16_t)smp;
                img_P++;
              }
          }
      }    

    /* Read any post-image chunks: */
    png_read_end(pr, pe);

    for (uint32_t row = 0;  row < rows; row++) { free(png_dataP[row]); }
    free(png_dataP);
    png_destroy_read_struct(&pr, &pi, &pe);
    /* png_free_data(pr, pi, PNG_FREE_ALL, -1); */
    return img; 
  }
  
