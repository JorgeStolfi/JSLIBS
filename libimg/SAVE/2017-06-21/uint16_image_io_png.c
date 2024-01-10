/* See {uint16_image_io_jpeg.h} */
/* Last edited on 2017-06-20 20:51:02 by stolfilocal */

/* Created by R. Minetto (IC-UNICAMP) as {ipng.c} sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <setjmp.h>

#include <png.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_io_png.h>
#include <sample_conv.h>

/* INTERNAL PROTOTYPES */

void jspng_dump_info(FILE *wr, const char *func, char *label, png_structp pr, png_infop pi);

/* IMPLEMENTATIONS */

uint16_image_t *uint16_image_io_png_read (char *name, double *gammaP, sample_uint32_t imaxval[], bool_t verbose)
  {
    FILE *rd = open_read (name, verbose);
    uint16_image_t *img = uint16_image_io_png_fread(rd, gammaP, imaxval, verbose);
    if (rd != stdin) { fclose(rd); }
    return img;
  }

void uint16_image_io_png_write (char *name, uint16_image_t *img, double gamma, bool_t verbose)
  {
    FILE *wr = open_write (name, verbose);
    uint16_image_io_png_fwrite(wr, img, gamma, verbose);
    if ((wr != stdout) && (wr != stderr)) { fclose(wr); }
  }

#define uint16_image_io_png_MAGIC_BYTES 8
/* Number of bytes in the PNG file signature. */

#define uint16_image_io_png_MAX_SIZE (2147483647u)
/* Max width and height of a PNG file that can be read by this interface. */

uint16_image_t *uint16_image_io_png_fread (FILE *rd, double *gammaP, uint32_t imaxval[], bool_t verbose)
  {
    bool_t debug = TRUE;
    
    /* Check the PNG signature: */
    unsigned char sig[uint16_image_io_png_MAGIC_BYTES];
    fread(sig, 1, uint16_image_io_png_MAGIC_BYTES, rd);
    demand(png_sig_cmp(sig, 0, uint16_image_io_png_MAGIC_BYTES) == 0, "not a PNG file");

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
    png_set_sig_bytes(pr, uint16_image_io_png_MAGIC_BYTES);
    
    /* Read the header: */
    png_read_info(pr, pi);
    if (debug) { jspng_dump_info(stderr, __FUNCTION__, "before transforms", pr, pi); }

    /* Get image dimensions: */
    png_uint_32 cols = png_get_image_width(pr, pi);
    demand((cols > 0) && (cols <= uint16_image_io_png_MAX_SIZE), "invalid num of cols"); 
    
    png_uint_32 rows = png_get_image_height(pr, pi);
    demand((rows > 0) && (rows <= uint16_image_io_png_MAX_SIZE), "invalid num of rows"); 
     
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
    demand((chns >= 1) && (chns <= uint16_image_io_png_MAX_CHNS), "invalid channel count"); 
    
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
    for (uint32_t row = 0; row < rows; row++)
      { png_dataP[row] = notnull(malloc(bytes_per_png_row), "no mem"); }

    /* Read the image samples: */
    if (debug) { fprintf(stderr, "loading image ...\n"); }
    png_read_image(pr, png_dataP);
    if (debug) { fprintf(stderr, "image loaded\n"); }

    /* Read and convert the pixels, row by row: */
    assert((png_smp_bits == 8) || (png_smp_bits == 16));
    for (int32_t row = 0; row < rows; row++)
      { 
        /* Process one row of the image: */
        png_bytep png_P = png_dataP[row];
        uint16_t *img_P = img->smp[row];
        for (int32_t col = 0; col < cols; col++)
          { for (int32_t chn = 0; chn < chns; chn++)
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

    for (uint32_t row = 0; row < rows; row++) { free(png_dataP[row]); }
    free(png_dataP);
    png_destroy_read_struct(&pr, &pi, &pe);
    /* png_free_data(pr, pi, PNG_FREE_ALL, -1); */
    return img; 
  }
  
void uint16_image_io_png_fwrite (FILE *wr, uint16_image_t *img, double gamma, bool_t verbose)
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
    for (uint32_t k = 1; k < 16; k++) 
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

void jspng_dump_info(FILE *wr, const char *func, char *label, png_structp pr, png_infop pi)
  { 
    fprintf(wr, "------------------------------------------------------------\n");
    fprintf(wr, "PNG file info - %s - %s\n", func, label);
    fprintf(wr, "image_width = %u\n", (uint32_t)png_get_image_width(pr, pi));
    fprintf(wr, "image_height = %u\n", (uint32_t)png_get_image_height(pr, pi));
    fprintf(wr, "channels = %u\n", (uint32_t)png_get_channels(pr, pi));
    fprintf(wr, "bit_depth = %d\n", png_get_bit_depth (pr, pi));
    /* fprintf(wr, "usr_bit_depth = %d\n", (int32_t)pr->usr_bit_depth; */
    fprintf(wr, "color_type = %d\n", png_get_color_type(pr, pi));
    fprintf(wr, "interlace_type = %d\n", png_get_interlace_type(pr, pi));
    fprintf(wr, "compression_type = %d\n", png_get_compression_type(pr, pi));
    fprintf(wr, "filter_type = %d\n", png_get_filter_type(pr, pi));
    
    png_color_8p sBIT; 
    if (png_get_sBIT(pr, pi, &sBIT)) 
      { fprintf(wr, "sBIT = %16p = ( %u %u %u %u %u )\n", (void *)sBIT, sBIT->gray, sBIT->red, sBIT->green, sBIT->blue, sBIT->alpha); }
    else
      { fprintf(wr, "sBIT not specified\n"); }
    
    int has_tRNS = png_get_valid(pr, pi, PNG_INFO_tRNS);
    fprintf(wr, "has tRNS = %d\n", has_tRNS);
    
    double gamma;
    if (png_get_gAMA(pr, pi, &gamma))
      { fprintf(wr, "gAMA = %25.16e\n", gamma); }
    else
      { fprintf(wr, "gAMA not specified\n"); }
    
    fprintf(wr, "channels = %d\n", png_get_channels(pr, pi));
    fprintf(wr, "rowbytes = %u\n", (uint32_t)png_get_rowbytes(pr, pi));
    png_bytep sg = (png_bytep)png_get_signature(pr, pi);
    if (sg != NULL)
      { fprintf(wr, "signature = %16p = %02x %02x %02x %02x %02x %02x %02x %02x\n", (void *)sg, sg[0], sg[1], sg[2], sg[3], sg[4], sg[5], sg[6], sg[7]); }
    else
      { fprintf(wr, "png_get_signature failed\n"); }
    /* fprintf(wr, " = %d\n", png_get_(pr, pi)); */
    fprintf(wr, "------------------------------------------------------------\n");

  }
