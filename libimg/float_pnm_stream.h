#ifndef float_pnm_stream_H
#define float_pnm_stream_H

/* I/O streams over a PBM/PGM/PPB image with multi-row float-format buffer. */
/* Last edited on 2017-06-21 23:48:32 by stolfilocal */

#include <jspnm.h>
#include <float_image_buffer.h>
#include <bool.h>

typedef struct float_pnm_stream_t
  { /* Image file data: */
    pnm_format_t format;       /* Image file format (PGM_FORMAT, RPGM_FORMAT, etc.). */
    int rows;                  /* Image width. */
    int cols;                  /* Image height. */
    int chns;                  /* Image channels. */
    uint16_t maxval;           /* Maximum pixel value */
    bool_t raw;                /* TRUE if `raw' format variant. */
    bool_t bits;               /* TRUE for PBM format (`raw' or `plain'). */
    /* Sample conversion parameters: */
    bool_t isMask;             /* Determines details of sample conversion. */
    uint32_t badval;           /* Integer sample value meaning `invalid', or {>maxval} if none. */
    double *ftb;               /* Table to convert from integer to float samples. */
    uint16_t *smp;             /* Buffer for one row of integer pixels. */
    float_image_buffer_t *buf; /* Buffer for one or more rows of floated pixels. */
  } float_pnm_stream_t;
  /* A {float_pnm_stream_t} is an entity that allows reading or writing PNM
    image files in a row-by-row fashion. It provides buffer storage
    for a row {smp[0..chns*cols-1]} of samples from the file, and one or more 
    rows of float samples {buf}. See {float_image_buffer.h} for details.
    
    When reading, {buf->ynext} is the next image row to read from the
    file, of 0 if no rows were read yet. When writing, {buf->yfrst} is
    the next row to be written to the image file, or {buf->rows} if the
    entire image was written. The stream is empty when
    {buf->yfrst==buf->ynext} .*/

float_pnm_stream_t *float_pnm_stream_new(bool_t isMask, uint32_t badval);
  /* Allocates a {float_pnm_stream_t} descriptor record and saves the
    {isMask} and {badval} fields in it.  The pointers {smp,ftb,buf}
    are set to {NULL}. The other fields are set to arbitary values. */

void float_pnm_stream_free(float_pnm_stream_t *str);
  /* Frees all storage associated with {str}, including 
    {*(str.smp)}, {*(str.ftb)}, {*(str.buf)}, and the 
    descriptor {*str} itself. */

#endif
