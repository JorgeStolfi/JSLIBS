/* {jsjpeg.h} - definitions and tools for JPEG image files.  */
/* Last edited on 2017-06-23 01:43:57 by stolfilocal */

#ifndef jsjpeg_H
#define jsjpeg_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#undef FALSE
#undef TRUE
#include <jpeglib.h>
#undef FALSE
#undef TRUE

#include <bool.h>
  
#define jsjpeg_BITS_PER_SAMPLE (BITS_IN_JSAMPLE)
  /* The number of bits per sample in the JPEG file as defined by {jpeglib.h} (8 or 12).

    In theory, JPEG files may use either 8 or 12 bits per sample.
    However, in order to read those files it would be necessary to
    recompile the library with {BITS_IN_JSAMPLE=12}; and then it would
    not be able to read 8-bit JPEG files. (WHAT A CROCK!). Fortunately
    no one seems to generate 12-bit JPEGs, so it may be safe to ignore
    that possibility. */

#define jsjpeg_MAXVAL (MAXJSAMPLE)
  /* The max sample value as as defined by {jpeglib.h} (255 or 4095). */

#endif
