#ifndef float_image_io_pnm_stream_H
#define float_image_io_pnm_stream_H

/* I/O streams over a PBM/PGM/PPB image with multi-row float-format buffer. */
/* Last edited on 2017-06-20 21:29:42 by stolfilocal */

#include <jspnm.h>
#include <float_image_buffer.h>
#include <bool.h>

typedef struct float_image_io_pnm_stream_t
  { /* Image file data: */
    pnm_format_t format;       /* Image file format (PGM_FORMAT, RPGM_FORMAT, etc.). */
    int rows;                  /* Image width. */
    int cols;                  /* Image height. */
    int chns;                  /* Image channels. */
    uint16_t maxval;       /* Maximum pixel value */
    bool_t raw;                /* TRUE if `raw' format variant. */
    bool_t bits;               /* TRUE for PBM format (`raw' or `plain'). */
    /* Sample conversion parameters: */
    bool_t isMask;             /* Determines details of sample conversion. */
    uint32_t badval;           /* Integer sample value meaning `invalid', or {>maxval} if none. */
    double *ftb;               /* Table to convert from integer to float samples. */
    uint16_t *smp;         /* Buffer for one row of raw integer pixels. */
    float_image_buffer_t *buf; /* Buffer for one or more rows of floated pixels. */
  } float_image_io_pnm_stream_t;
  /* A {float_image_io_pnm_stream_t} is an entity that allows reading or writing PNM
    image files in a row-by-row fashion. It provides buffer storage
    for a horizontal band of samples from the input image. The
    samples, converted to float values, are stored in the
    sub-structure {buf}. See {float_image_buffer.h} for details.
    
    When reading, {buf->ynext} is the next image row to read from the
    file, of 0 if no rows were read yet. When writing, {buf->yfrst} is
    the next row to be written to the image file, or {buf->rows} if the
    entire image was written. The stream is empty when
    {buf->yfrst==buf->ynext} .*/

/* ********************************************************************** */
/* INPUT STREAMS */

float_image_io_pnm_stream_t *float_image_read_pnm_stream_new(FILE *rd, bool_t isMask, uint32_t badval, int bufrows);
  /* Creates a new stream structure {str} primed for reading from the
    open PBM/PGM/PPM image file {rd}.
    
    Initializes the fields {rows,cols,chns,maxval,raw,bits,format} of
    {str} by reading the image header from {rd}. Saves the {isMask,badval}.
    Allocates the integer sample buffer {str->smp} and the float
    sample buffer {str->buf}. The latter will have space for {bufrows}
    consecutive rows, and will be initially empty.
    
    Also creates the conversion table {ftb} filled with
    {pnm_floatize(ival,maxval,isMask,badval)} from {jspnm.h}, for
    {ival} in {0..maxval}. Note that if {badval} is in the range
    {0..maxval}, any input samples equal to {badval} will be read as
    {NAN}. Use {badval==PNM_NO_BADVAL} to avoid this. */
    
double *float_image_read_pnm_stream_get_row(FILE *rd, float_image_io_pnm_stream_t *str, int y);
  /* Makes sure that row {y} of pixels from file {rd} is loaded in the
    buffer {str->buf}, and returns its address there.
   
   If row {y} it is not already stored in the buffer, loads into {str}
   every row from row {str->ynext} up to row {y}, with
   {float_image_read_pnm_stream_load_next_row} (using {str->smp} as a temproary
   buffer). May cause some earlier rows to be deallocated.
   
   Returns NULL if {y} is not in {0..str->rows-1}. Fails if {y} is a
   row that has been already deallocated from the buffer, i.e. if
   {y<str->buf->yini}. */
   
void float_image_read_pnm_stream_load_next_row(FILE *rd, float_image_io_pnm_stream_t *str);
  /* Reads the next image row (namely row {str->ynext}) from file {rd}
    into {str->smp}, converts its samples to floating point, stores them
    into {*(str->rowptr[k])} for the appropriate index {k}. May cause
    row {str->yfrst} to disappear from the stream. 
    Fails if {str->ynext >= str->rows}.

    Pixel samples are converted from integer to {double} with the table {str->ftab}. */

/* ********************************************************************** */
/* OUTPUT STREAMS */

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
  );
  /* Creates a new stream structure {str} primed for writing to the
    open PBM/PGM/PPM image file {wr}.
    
    Initializes the image stream {str} setting the
    {maxval,rows,cols,chns,isMask,badval} from the given parameters.
    Chooses a suitable {format} and sets {raw,bits} depending on those
    parameters and {forceplain}, then writes the corresponding PNM
    image file header to {wr}. Allocates the integer sample buffer
    {str->smp} and the float sample buffer {str->buf}. The latter will
    have space for at most {bufrows} rows, and will be initially
    empty. The table {ftb} is not allocated. */

double *float_image_write_pnm_stream_get_row(FILE *rd, float_image_io_pnm_stream_t *str, int y);
  /* Makes sure that there is a sample vector in the buffer {str->buf} to store
   row {y} of the image, and returns the address of that vector.
   
   If row {y} is not yet in the buffer, allocates new rows from
   {str->buf->ylim} to {y} inclusive. This may require writing out
   some earlier rows, with {float_image_write_pnm_stream_dump_first_row}, and
   deallocating them. The newly allocated buffer rows are initialized
   with zeros. 
   
   Returns NULL if {y} is outside the range {0..rows-1}. Fails if
   requested to return a row that has already been written out, i.e.
   if {y < y->buf->yini}. */

void float_image_write_pnm_stream_dump_first_row(FILE *wr, float_image_io_pnm_stream_t *str);
  /* Quantizes the first pixel row in the stream, namely row
    {str->yfrst}, writes it to {wr}, and removes it from the stream. Fails if
    {str->yfrst >= str->rows}.  To flush the streamed lines,
    up to row {y}, call this procedue until {y < str->yfrst}.
   
   Pixel samples are converted from {double} to integers with
   {pnm_quantize(fval,str->maxval,str->isMask,str->badval)} from {jspnm.h}. Note that
   the {str->badval} integer value, if in {0..str->maxval}, is used to
   encode {NAN} float values, exclusively. Use {badval==PNM_NO_BADVAL} to avoid this. */

/* ********************************************************************** */
/* FREEING STREAM STORAGE */

void fpnm_stream_free(float_image_io_pnm_stream_t *str);
  /* Frees all storage associated with {str}, icluding {*str} itself. */

#endif
