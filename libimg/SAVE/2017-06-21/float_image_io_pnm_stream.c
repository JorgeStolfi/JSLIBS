/* See {float_image_from_to_uint16_image_buffer.h}. */
/* Last edited on 2017-06-20 20:48:57 by stolfilocal */

#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <float_image_buffer.h>

#include <float_image_io_pnm_stream.h>

/* Internal prototypes: */

float_image_io_pnm_stream_t *float_image_read_pnm_stream_new(FILE *rd, bool_t isMask, uint32_t badval, int bufrows)
  { /* Allocate top record: */
    float_image_io_pnm_stream_t *str = (float_image_io_pnm_stream_t *)notnull(malloc(sizeof(float_image_io_pnm_stream_t)), "no mem");
    /* Read the input file header: */
    pnm_read_header
     ( rd, &(str->cols), &(str->rows), &(str->chns), &(str->maxval), 
       &(str->raw), &(str->bits), &(str->format)
     );
    str->smp = uint16_image_alloc_pixel_row(str->cols, str->chns);
    str->badval = badval;
    str->isMask = isMask;
    str->ftb = pnm_make_floatize_table(str->maxval, isMask, badval);
    str->buf = float_image_buffer_new(str->chns, str->cols, str->rows, bufrows);
    return str;
  }

void fpnm_stream_free(float_image_io_pnm_stream_t *str)
  { if (str != NULL)
      { if (str->smp != NULL) { free(str->smp); str->smp = NULL; }
        if (str->ftb != NULL) { free(str->ftb); str->ftb = NULL; }
        if (str->buf != NULL) { float_image_buffer_free(str->buf); }
        free(str); str = NULL; 
      }
  }

double *float_image_read_pnm_stream_get_row(FILE *rd, float_image_io_pnm_stream_t *str, int y) 
  { /* Row index must be valid: */
    if ((y < 0) || (y >= str->rows)) { return NULL; }
    /* Roll buffer forward until row {y} is in buffer: */
    demand(float_image_buffer_row_pos(str->buf, y) >= 00, "row has been discarded and cannot be read again");
    while (float_image_buffer_row_pos(str->buf, y) == +1)
      { float_image_read_pnm_stream_load_next_row(rd, str); }
    assert(float_image_buffer_row_pos(str->buf, y) == 00);
    return float_image_buffer_get_row(str->buf, y);
  }

void float_image_read_pnm_stream_load_next_row(FILE *rd, float_image_io_pnm_stream_t *str) 
  { /* Get row in {ibuf} for converted row {str->ynext}: */
    int y = str->buf->ylim;
    demand(y < str->rows, "no more rows to read");
    assert(float_image_buffer_row_pos(str->buf, y) == +1);
    float_image_buffer_advance(str->buf);
    assert(float_image_buffer_row_pos(str->buf, y) == 00);
    /* Read one row of sample values from input image into {str->smp}: */
    pnm_read_pixels(rd, str->smp, str->cols, str->chns, str->maxval, str->raw, str->bits);
    int nspr = str->chns * str->cols; /* Number of samples per row. */
    /* Convert samples {str->smp[0..nspr-1]} to [0_1] scale, save in {str->buf}: */
    uint16_t *sP = str->smp; /* Scans raw samples. */
    double *dP = float_image_buffer_get_row(str->buf, y); /* Start of row {yb} in {ibuf}. */
    assert(dP != NULL);
    int k;
    for (k = 0; k < nspr; k++, dP++, sP++)
      { demand((*sP) <= str->maxval, "invalid pixel value");
        (*dP) = str->ftb[(*sP)];
      }
  }

float_image_io_pnm_stream_t *float_image_write_pnm_stream_new 
  ( FILE *wr, 
    uint16_t maxval, 
    int rows, 
    int cols, 
    int chns, 
    bool_t isMask,
    uint32_t badval, 
    bool_t forceplain,
    int bufrows
  )
  { /* Allocate top record: */
    float_image_io_pnm_stream_t *str = (float_image_io_pnm_stream_t *)notnull(malloc(sizeof(float_image_io_pnm_stream_t)), "no mem");
    /* Select output format and write image header: */
    str->chns = chns;
    str->rows = rows;
    str->cols = cols;
    str->maxval = maxval;
    pnm_choose_output_format
      ( maxval, chns, forceplain,
        &(str->format), &(str->raw), &(str->bits)
      );
    pnm_write_header(wr, cols, rows, str->maxval, str->format);
    str->smp = uint16_image_alloc_pixel_row(str->cols, str->chns);
    str->isMask = isMask;
    str->badval = badval; 
    str->ftb = NULL;
    /* Allocate the floated-row buffer: */
    str->buf = float_image_buffer_new(chns, cols, rows, bufrows);
    return str;
  }

double *float_image_write_pnm_stream_get_row(FILE *wr, float_image_io_pnm_stream_t *str, int y)
  { /* Row index must be valid: */
    if ((y < 0) || (y >= str->rows)) { return NULL; }
    /* Roll buffer forward until row {y} is in buffer: */
    demand(float_image_buffer_row_pos(str->buf, y) >= 00, "row has been written out and cannot be created again");
    while (float_image_buffer_row_pos(str->buf, y) == +1) 
      { float_image_write_pnm_stream_dump_first_row(wr, str); }
    /* Did we succeed? */
    assert(float_image_buffer_row_pos(str->buf, y) == 00);
    return float_image_buffer_get_row(str->buf, y);
  }

void float_image_write_pnm_stream_dump_first_row(FILE *wr, float_image_io_pnm_stream_t *str) 
  { int y = str->buf->yini;
    demand(y < str->rows, "no more rows to write out");
    assert(float_image_buffer_row_pos(str->buf, y) == 00);
    /* Convert first row, write it, clear it: */
    double *dP = float_image_buffer_get_row(str->buf, y); /* Start of row {yb} in {ibuf}. */
    uint16_t *sP = str->smp; /* Scans raw samples. */
    int k;
    int nspr = str->chns * str->cols;
    for (k = 0; k < nspr; k++, dP++, sP++)
      { (*sP) = pnm_quantize((*dP), str->maxval, str->isMask, str->badval); (*dP) = 0.0; }
    pnm_write_pixels(wr, str->smp, str->cols, str->chns, str->maxval, str->raw, str->bits);
    /* Roll the buffer forward: */
    float_image_buffer_advance(str->buf);
    assert(float_image_buffer_row_pos(str->buf, y) == -1);
  }
