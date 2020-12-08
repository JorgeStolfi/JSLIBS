/* Bit stream tools for toy cryptography. */
/* Last edited on 2013-09-04 18:21:34 by stolfi */

#ifndef cryptoy_bstream_H
#define cryptoy_bstream_H

#include <stdio.h>

#define cryptoy_bstream_H_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

#define cryptoy_bstream_BUF_SIZE 1024
  /* Size of read/write buffers. */
  
#define cryptoy_bstream_EOF (-1)
  /* Return result to signal the end of a bit stream. */

typedef struct cryptoy_bstream_rec_t *cryptoy_bstream_t;
  /* A handle to a bit stream descriptor (input or output). */

cryptoy_bstream_t cryptoy_bstream_from_file(FILE *f);
  /* Allocates a bit stream descriptor and initializes it to read bits
    from or write bits to file {f}. Assumes {f} has been opened in the
    proper direction. */
    
int cryptoy_bstream_copy(cryptoy_bstream_t rd, cryptoy_bstream_t wr);    
  /* The argument {rd} must be an input bitstream and {wr} must be
    an output bitstream. Reads {rd} to end-of-file
    and writes its contents to {wr}. Returns {cryptoy_ERROR} if an i/o error occurred
    on either stream, otherwise returns the number of bits copied. */

int cryptoy_bstream_read_bit(cryptoy_bstream_t rd);
  /* Reads the next bit from {rd}. Returns {cryptoy_bstream_EOF} if the
    stream is exhausted. In case of read error returns {cryptoy_bstream_ERROR}. */
    
int cryptoy_bstream_write_bit(cryptoy_bstream_t wr, int val);
  /* Writes {val} (either 0  or 1) to {wr} as the next bit. 
    In case of write error returns {cryptoy_bstream_ERROR}. */
    
void cryptoy_bstream_transfer_reminder(cryptoy_bstream_t bs0, cryptoy_bstream_t bs1);
  /* The streams {bs0} and {bs1} must be output bstreams, and the
    total number of bits written to them must be a multiple of 8. If
    the last byte of {bs1} has between 1 and 7 bits, then un-writes those
    remaining bits and writes them to {bs0}. Otherwise does nothing.
    Does not raise any i/o errors. */

int cryptoy_bstream_flush(cryptoy_bstream_t wr);
  /* The stream {wr} must be an output stream, and the number of bytes
    already written to it must be a multiple of 8. Forces all those
    bits to be written to the underlying UNIX file, if any. */
    
void cryptoy_bstream_free(cryptoy_bstream_t bs);
  /* Releases the storage used by descriptor {bs} including all
     internal working storage.  The attached UNIX file, if 
     any, is not closed.  If {bs} is an output descriptor,
     the client must call {cryptoy_bstream_flush} first, 
     otherwise a number of trailing bits may be lost. */

#define cryptoy_bstream_H_rights \
  cryptoy_bstream_H_copyright ".\n" \
  "This file is provided 'as is'; the author and his employer are" \
  " not liable for any damages that may result from its use.  This" \
  " file can be used and modified for any purpose as long as" \
  " the copyright and this 'rights' note are preserved in" \
  " all copies and derived versions.";

#endif
