/* See {float_pnm_read_stream.h}. */
/* Last edited on 2024-12-26 12:41:09 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <float_image_buffer.h>
#include <float_pnm_stream.h>

#include <float_pnm_read_stream.h>

/* Internal prototypes: */

float_pnm_stream_t *float_pnm_read_stream_new(FILE *rd, bool_t isMask, uint32_t badval, uint32_t bufrows)
  { /* Allocate top record: */
    float_pnm_stream_t *str = float_pnm_stream_new(isMask, badval);
    /* Read the input file header: */
    pnm_read_header
     ( rd, &(str->cols), &(str->rows), &(str->chns), &(str->maxval), 
       &(str->raw), &(str->bits), &(str->format)
     );
    str->smp = uint16_image_alloc_pixel_row(str->cols, str->chns);
    str->ftb = pnm_make_floatize_table(str->maxval, isMask, badval);
    str->buf = float_image_buffer_new((int32_t)str->chns, (int32_t)str->cols, (int32_t)str->rows, (int32_t)bufrows);
    return str;
  }

double *float_pnm_read_stream_get_row(FILE *rd, float_pnm_stream_t *str, int32_t y) 
  { /* Row index must be valid: */
    if ((y < 0) || (y >= str->rows)) { return NULL; }
    /* Roll buffer forward until row {y} is in buffer: */
    demand(float_image_buffer_row_pos(str->buf, y) >= 00, "row has been discarded and cannot be read again");
    while (float_image_buffer_row_pos(str->buf, y) == +1)
      { float_pnm_read_stream_load_next_row(rd, str); }
    assert(float_image_buffer_row_pos(str->buf, y) == 00);
    return float_image_buffer_get_row(str->buf, y);
  }

void float_pnm_read_stream_load_next_row(FILE *rd, float_pnm_stream_t *str) 
  { /* Get row in {ibuf} for converted row {str->ynext}: */
    int32_t y = str->buf->ylim;
    demand(y < str->rows, "no more rows to read");
    assert(float_image_buffer_row_pos(str->buf, y) == +1);
    float_image_buffer_advance(str->buf);
    assert(float_image_buffer_row_pos(str->buf, y) == 00);
    /* Read one row of sample values from input image into {str->smp}: */
    pnm_read_pixels(rd, str->smp, str->cols, str->chns, str->maxval, str->raw, str->bits);
    uint32_t nspr = str->chns * str->cols; /* Number of samples per row. */
    /* Convert samples {str->smp[0..nspr-1]} to [0_1] scale, save in {str->buf}: */
    uint16_t *sP = str->smp; /* Scans raw samples. */
    double *dP = float_image_buffer_get_row(str->buf, y); /* Start of row {yb} in {ibuf}. */
    assert(dP != NULL);
    int32_t k;
    for (k = 0; k < nspr; k++, dP++, sP++)
      { demand((*sP) <= str->maxval, "invalid pixel value");
        (*dP) = str->ftb[(*sP)];
      }
  }

