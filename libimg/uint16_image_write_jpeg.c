/* See {uint16_image_write_jpeg.h} */
/* Last edited on 2024-12-05 07:45:56 by stolfi */

/* Created by R. Minetto (IC-UNICAMP) sometime in 2008--2009. */
/* Adapted by J. Stolfi (IC-UNICMP) on 2011-05-14. */

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

#include <uint16_image_write_jpeg.h>

void uint16_image_write_jpeg_named (char *name, uint16_image_t *img, int32_t quality, bool_t verbose)
  { FILE *wr = open_write(name, verbose);
    uint16_image_write_jpeg_file(wr, img, quality, verbose);
    if (wr != stdout) { fclose(wr); }
  }

void uint16_image_write_jpeg_file (FILE *wr, uint16_image_t *img, int32_t quality, bool_t verbose)
  {
    assert(jsjpeg_BITS_PER_SAMPLE == 8);
    assert(jsjpeg_MAXVAL == 255);

    /* !!! Extend to other types besides grayscale and RGB. !!! */

    /* Quality checking: */
    demand((quality > 0) && (quality <= 100), "invalid {quality}");
    
    /* Channel count checking: */
    uint32_t chns = (uint32_t)img->chns;
    uint32_t cols = (uint32_t)img->cols;
    uint32_t rows = (uint32_t)img->rows;
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
    jcmp.input_components = (int32_t)chns;
    jcmp.in_color_space = space;
    jpeg_set_defaults(&jcmp);
    /* jpeg_set_colorspace(&jcmp, space); */
    jpeg_set_quality(&jcmp, quality, TRUE);
    jpeg_start_compress(&jcmp, TRUE);
    
    #define jdemand(cond,msg) \
      if (! (cond)) { jpeg_abort_compress(&jcmp); demand(FALSE, (msg)); }

    /* Alocate the JPEG compression sample buffer: */
    uint32_t spr = chns * cols; /* Samples per scanline. */
    JSAMPARRAY buffer = (*jcmp.mem->alloc_sarray)((j_common_ptr) &jcmp, JPOOL_IMAGE, spr, 1);

    /* Write scanline by scanline: */
    while (jcmp.next_scanline < rows) {
      uint32_t r = jcmp.next_scanline;
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
