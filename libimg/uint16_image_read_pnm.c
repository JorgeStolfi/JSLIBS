/* See uint16_image_read.h */
/* Last edited on 2017-06-22 02:36:01 by stolfilocal */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>

#include <uint16_image_read_pnm.h>

uint16_image_t *uint16_image_read_pnm_named(char *name, bool_t verbose)
  { FILE *rd = open_read(name, verbose);
    uint16_image_t *img = uint16_image_read_pnm_file(rd);
    if (rd != stdin) { fclose(rd); }
    return img;
  }

uint16_image_t *uint16_image_read_pnm_file(FILE *rd)
  { /* We cannot rely on {pnm_readpnm} because PBMPLUS has different pixel layout. */
    /* Read file header: */
    int cols, rows, chns;
    uint16_t maxval;
    bool_t raw, bits;
    pnm_format_t format;
    pnm_read_header(rd, &cols, &rows, &chns, &maxval, &raw, &bits, &format);
    /* Allocate image and set maxval: */
    uint16_image_t *img = uint16_image_new(cols, rows, chns);
    img->maxval = maxval;
    /* Read pixels: */
    int row;
    for (row = 0; row < rows; row++)
      { pnm_read_pixels(rd, img->smp[row], cols, chns, maxval, raw, bits); }
    return(img);
  }
  
