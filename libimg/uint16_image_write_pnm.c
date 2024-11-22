/* See uint16_image_write.h */
/* Last edited on 2017-07-01 00:05:22 by stolfilocal */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>

#include <uint16_image_write_pnm.h>

void uint16_image_write_pnm_named(char *name, uint16_image_t *img, bool_t forceplain, bool_t verbose)
  { FILE *wr = open_write(name, verbose);
    uint16_image_write_pnm_file(wr, img, forceplain, verbose);
    if (wr != stdout) { fclose(wr); }
  }
    
void uint16_image_write_pnm_file(FILE *wr, uint16_image_t *img, bool_t forceplain, bool_t verbose)
  { /* We cannot rely on {pnm_writepnm} because PBMPLUS has different pixel layout. */
    if (img->maxval <= 0){ pnm_error("invalid maxval (%u)", img->maxval); }
    pnm_format_t format;
    bool_t raw;
    bool_t bits; 
    pnm_choose_output_format(img->maxval, img->chns, forceplain, &format, &raw, &bits);
    pnm_write_header(wr, img->cols, img->rows, img->maxval, format);
    int row;
    for (row = 0; row < img->rows; row++)
      { pnm_write_pixels(wr, img->smp[row], img->cols, img->chns, img->maxval, raw, bits); }
    fflush(wr);
  }

