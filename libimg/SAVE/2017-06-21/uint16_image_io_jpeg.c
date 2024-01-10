/* See {uint16_image_io_jpeg.h} */
/* Last edited on 2017-06-20 20:50:01 by stolfilocal */

/* Created by R. Minetto (IC-UNICAMP) sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <jpeglib.h>
#undef FALSE
#undef TRUE

#include <affirm.h>
#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_io_jpeg.h>

uint16_image_t *uint16_image_io_jpeg_read (char *name, bool_t verbose, int32_t *spaceP)
  { FILE *rd = open_read(name, verbose);
    uint16_image_t *img = uint16_image_io_jpeg_fread(rd, verbose, spaceP);
    if (rd != stdin) { fclose(rd); }
    return img;
  }

void uint16_image_io_jpeg_write (char *name, uint16_image_t *img, int32_t quality, bool_t verbose)
  { FILE *wr = open_write(name, verbose);
    uint16_image_io_jpeg_fwrite(wr, img, quality, verbose);
    if (wr != stdout) { fclose(wr); }
  }

uint16_image_t *uint16_image_io_jpeg_fread (FILE *rd, bool_t verbose, int32_t *spaceP)
  {
    /* Make sure that the JPEG library was compiled with 8-bit samples: */
    assert(jsjpeg_BITS_PER_SAMPLE == 8);
    assert(jsjpeg_MAXVAL == 255);
    
    /* Create the decompressor structures: */
    struct jpeg_decompress_struct jdec;
    struct jpeg_error_mgr jerm;
    jdec.err = jpeg_std_error(&jerm);
    jpeg_create_decompress(&jdec);
    jpeg_stdio_src(&jdec, rd);
    
    /* Read the file header and start the decompressor: */
    demand(jpeg_read_header(&jdec, TRUE) == 1, "{jpeg_read_header} failed");
    (void) jpeg_start_decompress(&jdec);  /* !!! Should check the return value? !!! */
    
    #define jdemand(cond,msg) \
      if (! (cond)) { jpeg_abort_decompress(&jdec); demand(FALSE, (msg)); }
    
    /* Get the image size: */
    int32_t cols = jdec.output_width;
    int32_t rows = jdec.output_height;
    int32_t chns = jdec.output_components;
    jdemand((chns > 0) && (chns <= 4), "invalid channel count");

    /* Return the image colorspace: */
    (*spaceP) = jdec.jpeg_color_space;

    /* Alocate the JPEG decompression sample buffer: */    
    int32_t spr = cols * chns;   /* Samples per row in decompressor output. */
    JSAMPARRAY buffer = (*jdec.mem->alloc_sarray)((j_common_ptr)&jdec, JPOOL_IMAGE, spr, 1);

    /* Alocate the in-core image: */    
    uint16_image_t *img = uint16_image_new(cols, rows, chns);
    img->maxval = jsjpeg_MAXVAL;

    /* Read scanline by scanline: */
    while (jdec.output_scanline < rows) {
      int32_t r = jdec.output_scanline;
      uint32_t nread = jpeg_read_scanlines(&jdec, buffer, 1);
      jdemand(nread == 1, "{jpeg_read_scanlines} failed");
      uint16_t *img_row = img->smp[r];
      JSAMPLE *bu_row = buffer[0];
      int32_t k;
      for (k = 0; k < spr; k++) { img_row[k] = (int32_t) bu_row[k]; }
    }
    
    /* Cleanup: */
    (void) jpeg_finish_decompress(&jdec);  /* !!! Should check the return value? !!! */
    jpeg_destroy_decompress(&jdec);

    return img;
    
    #undef jdemand
  }

void uint16_image_io_jpeg_fwrite (FILE *wr, uint16_image_t *img, int32_t quality, bool_t verbose)
  {
    assert(jsjpeg_BITS_PER_SAMPLE == 8);
    assert(jsjpeg_MAXVAL == 255);

    /* !!! Extend to other types besides grayscale and RGB. !!! */

    /* Quality checking: */
    demand((quality > 0) && (quality <= 100), "invalid {quality}");
    
    /* Channel count checking: */
    int32_t chns = img->chns;
    int32_t cols = img->cols;
    int32_t rows = img->rows;
    demand((chns == 1) || (chns == 3), "invalid channel count");
    
    uint16_t imaxval = img->maxval;  /* Image's {maxval}. */
    demand(imaxval + 0 <= 65535, "maxval is too big"); /* Paranoia, since {uint16_t} is 16 bits. */
    demand(imaxval > 0, "maxval is zero");
    
    /* Determine the max JPEG sample value: */
    uint32_t omaxval = jsjpeg_MAXVAL;  /* Max output sample value. */

    /* Determine the sample scaling factor: */
    if (imaxval != omaxval) 
      { fprintf(stderr, "%s: ", __FUNCTION__);
        fprintf(stderr, "rescaling from {0..%d} to {0..%d}\n", imaxval, omaxval);
        if (imaxval > omaxval)
          { fprintf(stderr, "!! warning: samples will be downscaled with loss of precision.\n"); }
        else if ((omaxval % imaxval) != 0)
          { fprintf(stderr, "!! warning: scaling will not be exact.\n"); }
      }
    uint32_t iround = imaxval / 2u;  /* Bias for sample conversion rounding. */
    assert(iround < imaxval);

    /* Create the compressor structures: */
    struct jpeg_compress_struct jcmp;
    struct jpeg_error_mgr jerm;
    jcmp.err = jpeg_std_error(&jerm);
    jpeg_create_compress(&jcmp);
    jpeg_stdio_dest(&jcmp, wr);

    /* Start the compressor: */
    J_COLOR_SPACE space = (chns == 1 ? JCS_GRAYSCALE : JCS_RGB);
    jcmp.image_width = cols;
    jcmp.image_height = rows;
    jcmp.input_components = chns;
    jcmp.in_color_space = space;
    jpeg_set_defaults(&jcmp);
    /* jpeg_set_colorspace(&jcmp, space); */
    jpeg_set_quality(&jcmp, quality, TRUE);
    jpeg_start_compress(&jcmp, TRUE);
    
    #define jdemand(cond,msg) \
      if (! (cond)) { jpeg_abort_compress(&jcmp); demand(FALSE, (msg)); }

    /* Alocate the JPEG compression sample buffer: */
    int32_t spr = chns * cols; /* Samples per scanline. */
    JSAMPARRAY buffer = (*jcmp.mem->alloc_sarray)((j_common_ptr) &jcmp, JPOOL_IMAGE, spr, 1);

    /* Write scanline by scanline: */
    while (jcmp.next_scanline < rows) {
      int32_t r = jcmp.next_scanline;
      uint16_t *img_row = img->smp[r];
      JSAMPLE *bu_row = buffer[0];
      int32_t k;
      for (k = 0; k < spr; k++) 
        { uint32_t smp = img_row[k];
          /* Apply scale factor: */
          if (imaxval != omaxval) { smp = (uint32_t)((omaxval*(uint64_t)smp + iround)/imaxval); }
          bu_row[k] = (JSAMPLE)smp;
        }

      uint32_t nwrote = jpeg_write_scanlines(&jcmp, buffer, 1);
      jdemand(nwrote == 1, "{jpeg_read_scanlines} failed");
    }

    /* Cleanup: */
    jpeg_finish_compress(&jcmp);
    jpeg_destroy_compress(&jcmp);
    
    #undef jdemand
  }
