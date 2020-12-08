/* See {cryptoy_bstream.h} */
/* Last edited on 2013-09-04 18:28:27 by stolfi */

#define cryptoy_bstream_C_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <cryptoy.h>
#include <cryptoy_bstream.h>

typedef struct cryptoy_bstream_rec_t 
  { FILE *f;         /* Attached file. */
    unsigned char buf[cryptoy_bstream_BUF_SIZE]; 
    int nbuf;        /* Number of useful bytes in buf. */
    int nchars;      /* Number of full bytes consumed/created */
    unsigned char C; /* Next icomplete byte. */
    int nbits;       /* Number of bits consumed/created in {C}. */
  } cryptoy_bstream_rec_t;
  /* In any case one has {0 <= nchars <= nbuf <= cryptoy_BUF_SIZE}.
    The bits in each byte are read/written from highest to lowest.
  
    In an input bit stream, the file {f} must be open for 
    reading. The last batch of characters that was read
    from the file {f} is {buf[0..nbuf-1]}. The part of that batch that
    remains to be consumed by the client consists of the HIGHEST
    {nbits} bits (0 to 7) of {C} followed by characters
    {buf[nchars..nbuf-1]}. The remaining {8-nbits} of {C} are zero. If
    {nbits = 0} then {C} is to be loaded from {buf[nchars]}, otherwise
    {nchars} is positive and {C} was taken from {buf[nchars-1]}.
    
    In an output stream, the file {f} must be open for 
    writing. Then {nbuf} is some positive number; the batch of
    data that is to be written to file {f} is bytes {buf[0..nchars-1]}
    followed by the LOWEST {nbits} bits (0 to 7) of {C}. If {nbits =
    0} then {nchars <= nbuf}, otherwise then {nchars < nbuf} and {C},
    when completed, will be stored into {buf[nchars]}. */
  
cryptoy_bstream_t cryptoy_bstream_from_file(FILE *f)
  { cryptoy_bstream_t bs = malloc(sizeof(cryptoy_bstream_t));
    assert(bs != NULL);
    bs->f = f;
    bs->nbuf = 0;  /* Ignored in output streams. */
    bs->nchars = 0;
    bs->nbits = 0;
    bs->C = 0;
    return bs;
  }
  
int cryptoy_bstream_read_bit(cryptoy_bstream_t rd)
  { 
    if (ferror(rd->f)) { return cryptoy_ERROR; }
    /* Make sure that {rd->nbits > 0}: */
    if (rd->nbits == 0)
      { if (rd->nchars >= rd->nbuf)
          { /* Read a new batch: */
            if (feof(rd->f)) { return cryptoy_bstream_EOF; }
            rd->nbuf = fread(rd->buf, 1, cryptoy_bstream_BUF_SIZE, rd->f);
            if (ferror(rd->f)) { return cryptoy_ERROR; }
            if (rd->nbuf == 0) { assert(feof(rd->f)); return cryptoy_bstream_EOF; }
            rd->nchars = 0;
          }
        assert(rd->nchars < rd->nbuf);
        rd->C = rd->buf[rd->nchars]; rd->nchars++;
        rd->nbits = 8;
      }
    /* Return the next bit of {rd->C}: */
    int res = (rd->C >> 7) & 1; rd->C <<= 1; rd->nbits--;
    return res;
  }
  
int cryptoy_bstream_write_bit(cryptoy_bstream_t wr, int val)
  {
    if (ferror(wr->f)) { return cryptoy_ERROR; }
    wr->nbuf = cryptoy_bstream_BUF_SIZE;
    /* Make sure that {wr->nbits < 8}: */
    if (wr->nbits == 8)
      { if (wr->nchars >= wr->nbuf)
          { /* write the current batch: */
            int nwr = fwrite(wr->buf, 1, wr->nbuf, wr->f);
            if (ferror(wr->f)) { return cryptoy_ERROR; }
            assert(! feof(wr->f));
            assert(nwr == wr->nbuf);
            wr->nchars = 0;
          }
        assert(wr->nchars < wr->nbuf);
        wr->buf[wr->nchars] = wr->C; wr->nchars++;
        wr->nbits = 0; wr->C = 0;
      }
    /* Append the bit {val} to {wr->C}: */
    wr->C <<= 1; wr->C |= (val & 1); wr->nbits++;
    return 0;
  }
  
void cryptoy_bstream_transfer_reminder(cryptoy_bstream_t wr0, cryptoy_bstream_t wr1)
  { 
    assert(((wr0->nbits + wr1->nbits) % 8) == 0);
    if ((wr1->nbits > 0) && (wr1->nbits < 8))
      { assert(wr0->nbits + wr1->nbits == 8);
        wr0->C <<= wr1->nbits; wr0->C |= wr1->C; wr0->nbits = 8; 
        wr1->nbits = 0; wr1->C = 0;
      }
  }

int cryptoy_bstream_copy(cryptoy_bstream_t rd, cryptoy_bstream_t wr)
  {
    int res = 0;
    while(1)
      { int bx = cryptoy_bstream_read_bit(rd);
        if (bx == cryptoy_ERROR) { return cryptoy_ERROR; }
        if (bx == cryptoy_bstream_EOF) { return res; }
        int resw = cryptoy_bstream_write_bit(wr, bx);
        if (resw == cryptoy_ERROR) { return cryptoy_ERROR; }
        assert(resw == 0);
        res++;
      }
  }

int cryptoy_bstream_flush(cryptoy_bstream_t wr)
  { 
    if (wr->nbits == 8)
      { assert(wr->nchars < wr->nbuf);
        wr->buf[wr->nchars] = wr->C; wr->nchars++; wr->nbits = 0;
      }
    assert(wr->nbits == 0);
    if (wr->nchars > 0)
      { /* write the current batch: */
        int nwr = fwrite(wr->buf, 1, wr->nchars, wr->f);
        if (ferror(wr->f)) { return cryptoy_ERROR; }
        assert(nwr == wr->nchars);
        assert(! feof(wr->f));
        wr->nchars = 0;
      }
    int res = fflush(wr->f);
    if (res != 0) { return cryptoy_ERROR; }
    return 0;
  }

void cryptoy_bstream_free(cryptoy_bstream_t bs)
  { 
    free(bs);
  }

#define cryptoy_bstream_C_rights \
  cryptoy_bstream_C_copyright ".\n" \
  "This file is provided 'as is'; the author and his employer are" \
  " not liable for any damages that may result from its use.  This" \
  " file can be used and modified for any purpose as long as" \
  " the copyright and this 'rights' note are preserved in" \
  " all copies and derived versions.";
