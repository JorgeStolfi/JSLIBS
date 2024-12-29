#ifndef float_pnm_read_stream_H
#define float_pnm_read_stream_H

/* Row-by-row reading of PBM/PGM/PPB image file with conversion to float samples. */
/* Last edited on 2024-12-26 12:37:22 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <float_pnm_stream.h>
#include <bool.h>

float_pnm_stream_t *float_pnm_read_stream_new(FILE *rd, bool_t isMask, uint32_t badval, uint32_t bufrows);
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
    
double *float_pnm_read_stream_get_row(FILE *rd, float_pnm_stream_t *str, int32_t y);
  /* Makes sure that row {y} of pixels from file {rd} is loaded in the
    buffer {str->buf}, and returns its address there.
   
   If row {y} it is not already stored in the buffer, loads into {str}
   every row from row {str->ynext} up to row {y}, with
   {float_pnm_read_stream_load_next_row} (using {str->smp} as a temproary
   buffer). May cause some earlier rows to be deallocated.
   
   Returns NULL if {y} is not in {0..str->rows-1}. Fails if {y} is a
   row that has been already deallocated from the buffer, i.e. if
   {y<str->buf->yini}. */
   
void float_pnm_read_stream_load_next_row(FILE *rd, float_pnm_stream_t *str);
  /* Reads the next image row (namely row {str->ynext}) from file {rd}
    into {str->smp}, converts its samples to floating point, stores them
    into {*(str->rowptr[k])} for the appropriate index {k}. May cause
    row {str->yfrst} to disappear from the stream. 
    Fails if {str->ynext >= str->rows}.

    Pixel samples are converted from integer to {double} with the table {str->ftab}. */

#endif
