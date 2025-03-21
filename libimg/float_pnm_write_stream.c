/* See {float_pnm_write_stream.h}. */
/* Last edited on 2024-12-26 12:43:16 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <float_image_buffer.h>
#include <float_pnm_stream.h>

#include <float_pnm_write_stream.h>

float_pnm_stream_t *float_pnm_write_stream_new 
  ( FILE *wr, 
    uint16_t maxval, 
    uint32_t rows, 
    uint32_t cols, 
    uint32_t chns, 
    bool_t isMask,
    uint32_t badval, 
    bool_t forceplain,
    uint32_t bufrows
  )
  { /* Allocate top record: */
    float_pnm_stream_t *str = float_pnm_stream_new(isMask, badval);
    /* Select output format and write image header: */
    str->chns = (uint32_t)chns;
    str->rows = (uint32_t)rows;
    str->cols = (uint32_t)cols;
    str->maxval = maxval;
    pnm_choose_output_format
      ( maxval, chns, forceplain,
        &(str->format), &(str->raw), &(str->bits)
      );
    pnm_write_header(wr, cols, rows, str->maxval, str->format);
    str->smp = uint16_image_alloc_pixel_row(str->cols, str->chns);
    str->ftb = NULL;
    /* Allocate the floated-row buffer: */
    str->buf = float_image_buffer_new((int32_t)chns, (int32_t)cols, (int32_t)rows, (int32_t)bufrows);
    return str;
  }

double *float_pnm_write_stream_get_row(FILE *wr, float_pnm_stream_t *str, int32_t y)
  { /* Row index must be valid: */
    if ((y < 0) || (y >= str->rows)) { return NULL; }
    /* Roll buffer forward until row {y} is in buffer: */
    demand(float_image_buffer_row_pos(str->buf, y) >= 00, "row has been written out and cannot be created again");
    while (float_image_buffer_row_pos(str->buf, y) == +1) 
      { float_pnm_write_stream_dump_first_row(wr, str); }
    /* Did we succeed? */
    assert(float_image_buffer_row_pos(str->buf, y) == 00);
    return float_image_buffer_get_row(str->buf, y);
  }

void float_pnm_write_stream_dump_first_row(FILE *wr, float_pnm_stream_t *str) 
  { int32_t y = str->buf->yini;
    demand(y < str->rows, "no more rows to write out");
    assert(float_image_buffer_row_pos(str->buf, y) == 00);
    /* Convert first row, write it, clear it: */
    double *dP = float_image_buffer_get_row(str->buf, y); /* Start of row {yb} in {ibuf}. */
    uint16_t *sP = str->smp; /* Scans raw samples. */
    int32_t k;
    uint32_t nspr = str->chns * str->cols;
    for (k = 0; k < nspr; k++, dP++, sP++)
      { (*sP) = pnm_quantize((*dP), str->maxval, str->isMask, str->badval); (*dP) = 0.0; }
    pnm_write_pixels(wr, str->smp, str->cols, str->chns, str->maxval, str->raw, str->bits);
    /* Roll the buffer forward: */
    float_image_buffer_advance(str->buf);
    assert(float_image_buffer_row_pos(str->buf, y) == -1);
  }
