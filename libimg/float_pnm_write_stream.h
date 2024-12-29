#ifndef float_pnm_write_stream_H
#define float_pnm_write_stream_H

/* Writing a PBM/PGM/PPB image file by rows from a float-format buffer. */
/* Last edited on 2024-12-26 12:42:18 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <float_pnm_stream.h>

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

double *float_pnm_write_stream_get_row(FILE *rd, float_pnm_stream_t *str, int32_t y);
  /* Makes sure that there is a sample vector in the float buffer {str->buf} to store
   row {y} of the image, and returns the address of that vector.
   
   If row {y} is not yet in the buffer, allocates new rows from
   {str->buf->ylim} to {y} inclusive. This may require writing out
   some earlier rows, with {float_pnm_write_stream_dump_first_row}, and
   deallocating them. The newly allocated buffer rows are initialized
   with zeros. 
   
   Returns NULL if {y} is outside the range {0..rows-1}. Fails if
   requested to return a row that has already been written out, i.e.
   if {y < y->buf->yini}. */

void float_pnm_write_stream_dump_first_row(FILE *wr, float_pnm_stream_t *str);
  /* Quantizes the first pixel row in the stream, namely row
    {str->yfrst}, writes it to {wr}, and removes it from the stream. Fails if
    {str->yfrst >= str->rows}.  To flush the streamed lines,
    up to row {y}, call this procedue until {y < str->yfrst}.
   
   Pixel samples are converted from {double} to integers with
   {pnm_quantize(fval,str->maxval,str->isMask,str->badval)} from {jspnm.h}. Note that
   the {str->badval} integer value, if in {0..str->maxval}, is used to
   encode {NAN} float values, exclusively. Use {badval==PNM_NO_BADVAL} to avoid this. */

#endif
