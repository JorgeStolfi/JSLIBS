/* See {uint16_image_read_jpeg.h} */
/* Last edited on 2017-06-23 01:44:52 by stolfilocal */

/* Created by R. Minetto (IC-UNICAMP) sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#undef FALSE
#undef TRUE
#include <jpeglib.h>
#undef FALSE
#undef TRUE

#include <affirm.h>
#include <jsfile.h>
#include <jsjpeg.h>
#include <uint16_image.h>

#include <uint16_image_read_jpeg.h>

uint16_image_t *uint16_image_read_jpeg_named(char *name, bool_t verbose, int32_t *spaceP)
  { FILE *rd = open_read(name, verbose);
    uint16_image_t *img = uint16_image_read_jpeg_file(rd, verbose, spaceP);
    if (rd != stdin) { fclose(rd); }
    return img;
  }

uint16_image_t *uint16_image_read_jpeg_file(FILE *rd, bool_t verbose, int32_t *spaceP)
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
