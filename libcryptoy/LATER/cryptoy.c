/* See {cryptoy.h} */
/* Last edited on 2013-09-04 18:29:23 by stolfi */

#define cryptoy_C_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <cryptoy.h>
#include <cryptoy_bstream.h>

int cryptoy_xor_files(FILE *f0, FILE *f1, FILE *f)
  { int N = cryptoy_bstream_BUF_SIZE; /* Buffer size */
    unsigned char abuf[N]; /* Buffer for input file {f0} and output file {f}. */
    uint64_t nw = 0;       /* Total bytes written to file {f}. */
    unsigned char bbuf[N]; /* Buffer for input file {f1}. */

    while(1)
      { /* Read a buffer-load from each file: */
        size_t ma = fread(abuf, 1, N, f0);
        if ((ma < N) && (! feof(f0))) { cryptoy_emsg("i/o error reading file 1"); return cryptoy_ERROR; }
        size_t mb = fread(bbuf, 1, N, f1);
        if ((mb < N) && (! feof(f1))) { cryptoy_emsg("i/o error reading file 2"); return cryptoy_ERROR; }
        if (ma > mb) { cryptoy_emsg("file 1 is longer than file 2"); return cryptoy_ERROR; }
        /* Combine buffers: */
        int mc = ma;
        int k;
        for (k = 0; k < mc; k++) { abuf[k] ^= bbuf[k]; }
        size_t mw = fwrite(abuf, 1, mc, f);
        nw += mw;
        if (mw < mc) { cryptoy_emsg("i/o error writing the output file"); return cryptoy_ERROR; } 
        if (mc < N) { break; }
      }
    int res = fflush(f);
    if (res != 0) { cryptoy_emsg("i/o error flushing the output file"); return cryptoy_ERROR; } 
    return nw;
  }

int cryptoy_bitsplit_file(FILE *f, FILE *k, FILE *f0, FILE *f1)
  {
    cryptoy_bstream_t sf = cryptoy_bstream_from_file(f);
    cryptoy_bstream_t sk = cryptoy_bstream_from_file(k);
    cryptoy_bstream_t sf0 = cryptoy_bstream_from_file(f0);
    cryptoy_bstream_t sf1 = cryptoy_bstream_from_file(f1);
    int res = 0;
    
    while(1)
      {
        /* Obtain a bit {fx} from {f}: */
        int fx = cryptoy_bstream_read_bit(sf);
        if (fx == cryptoy_ERROR) { cryptoy_emsg("i/o error reading data file"); res = cryptoy_ERROR; break; }
        if (fx == cryptoy_bstream_EOF) { break; }
        /* Obtain a bit {kx} from {k}: */
        int kx = cryptoy_bstream_read_bit(sk);
        if (kx == cryptoy_ERROR) { cryptoy_emsg("i/o error reading key file"); res = cryptoy_ERROR; break; }
        if (kx == cryptoy_bstream_EOF) { cryptoy_emsg("key file shorter than data file"); res = cryptoy_ERROR; break; }
        /* Write {fx} to the proper output file: */
        if (kx == 0)
          { int res0 = cryptoy_bstream_write_bit(sf0, fx);
            if (res0 < 0) { cryptoy_emsg("i/o error writing file 0"); res = cryptoy_ERROR; break; }
          }
        else
          { int res1 = cryptoy_bstream_write_bit(sf1, fx);
            if (res1 < 0) { cryptoy_emsg("i/o error writing file 1"); res = cryptoy_ERROR; break; }
          }
        res++;
      }
      
    if (res != cryptoy_ERROR)
      { /* File {f} is exhausted, complete {f0} with remainder bits of {f1} if any: */
        cryptoy_bstream_transfer_reminder(sf0, sf1);
        int res0 = cryptoy_bstream_flush(sf0);
        if (res0 < 0) { cryptoy_emsg("i/o error writing file 0"); }
        int res1 = cryptoy_bstream_flush(sf1);
        if (res1 < 0) { cryptoy_emsg("i/o error writing file 1"); }
        if ((res0 < 0) || (res1 < 0)) { res = cryptoy_ERROR; }
      }
    /* Free working storage: */
    cryptoy_bstream_free(sf1);
    cryptoy_bstream_free(sf0);
    cryptoy_bstream_free(sk);
    cryptoy_bstream_free(sf);
    
    return res;
  }
  
int cryptoy_bitmerge_files(FILE *f0, FILE *f1, FILE *k, FILE *f)
  {
    cryptoy_bstream_t sf0 = cryptoy_bstream_from_file(f0);
    cryptoy_bstream_t sf1 = cryptoy_bstream_from_file(f1);
    cryptoy_bstream_t sk = cryptoy_bstream_from_file(k);
    cryptoy_bstream_t sf = cryptoy_bstream_from_file(f);
    int res = 0;
    int fx0 = cryptoy_bstream_read_bit(sf0); /* Next bit from {f0}, or error/eof. */
    int fx1 = cryptoy_bstream_read_bit(sf1); /* Next bit from {f1}, or error/eof. */
    while(1)
      {
        /* Check for data exhaustion: */
        if ((fx0 < 0) || (fx1 < 0)) { break; }
        /* Obtain a bit {kx} that selects the source of the next bit: */
        int kx = cryptoy_bstream_read_bit(sk);
        if (kx == cryptoy_ERROR) { cryptoy_emsg("i/o error reading key file"); res = cryptoy_ERROR; break; }
        if (kx == cryptoy_bstream_EOF) { cryptoy_emsg("key file shorter than data files"); res = cryptoy_ERROR; break; }
        /* Get the next output but {fx} from the proper input file: */
        int fx;
        if (kx == 0)
          { fx = fx0;
            fx0 = cryptoy_bstream_read_bit(sf0);
          }
        else
          { fx = fx1;
            fx1 = cryptoy_bstream_read_bit(sf1);
          }
        /* Write to output stream: */
        int resf = cryptoy_bstream_write_bit(sf, fx);
        if (resf < 0) { cryptoy_emsg("i/o error writing the output file"); res = cryptoy_ERROR; break; }
        res++;
      }
     
    if (fx0 == cryptoy_ERROR) { cryptoy_emsg("i/o error reading file 0"); res = cryptoy_ERROR; }
    if (fx1 == cryptoy_ERROR) { cryptoy_emsg("i/o error reading file 1"); res = cryptoy_ERROR; }
    
    auto int copy_tail(int gx, cryptoy_bstream_t sg);
      /* Copies the bit {gx} and the rest of stream {sg} to {sf}. 
        If {gx == cryptoy_bstream_EOF} assumed {sg} is empty
        and does nothing. Returns the number of bits copied
        (possibly 0), or {cryptoy_ERROR} if an i/o error occurred. */
    
    int ntail = 0; /* Number of tail bits copied from {f0}, minus those copied from {f1}. */
    if (res != cryptoy_ERROR)
      { /* File {f1} is exhausted, copy the tail of {f0} the other file: */
        int resf = copy_tail(fx0, sf0);
        if (resf < 0) 
          { res = cryptoy_ERROR; }
        else
          { ntail = resf; res += resf; }
        if (res == cryptoy_ERROR) { cryptoy_emsg("i/o error while copying tail of file 0"); }
      }
    if (res != cryptoy_ERROR)
      { /* File {f0} is exhausted, copy the tail of {f1} the other file: */
        int resf = copy_tail(fx1, sf1);
        if (resf < 0) 
          { res = cryptoy_ERROR; }
        else
          { ntail -= resf; res += resf; }
        if (res == cryptoy_ERROR) { cryptoy_emsg("i/o error while copying tail of file 1"); }
      }

    if (res != cryptoy_ERROR)
      { if ((ntail < 0) || (ntail > 7)) { cryptoy_wmsg("input file lengths inconsistent with key file"); }
        res = cryptoy_bstream_flush(sf);
        if (res < 0) { cryptoy_emsg("i/o error while flushing output file"); res = cryptoy_ERROR; }
      }

    /* Free working storage: */
    cryptoy_bstream_free(sf);
    cryptoy_bstream_free(sk);
    cryptoy_bstream_free(sf1);
    cryptoy_bstream_free(sf0);
    return res;
       
    /* Internal implementations: */
    
    int copy_tail(int gx, cryptoy_bstream_t sg)
      { if (gx == cryptoy_bstream_EOF) 
          { return 0; }
        else
          { int resx = cryptoy_bstream_write_bit(sf, gx);
            if (resx < 0)
              { return cryptoy_ERROR; }
            else 
              { int resg = cryptoy_bstream_copy(sg, sf);
                if (resg < 0) 
                  { return cryptoy_ERROR; }
                else
                  { return 1 + resg; }
              }
          }
      }
  }


#define cryptoy_C_rights \
  cryptoy_C_copyright ".\n" \
  "This file is provided 'as is'; the author and his employer are" \
  " not liable for any damages that may result from its use.  This" \
  " file can be used and modified for any purpose as long as" \
  " the copyright and this 'rights' note are preserved in" \
  " all copies and derived versions.";
