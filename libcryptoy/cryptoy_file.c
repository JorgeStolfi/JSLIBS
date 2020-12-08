/* See {cryptoy.h} */
/* Last edited on 2013-10-31 01:45:40 by stolfilocal */

#define cryptoy_C_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <cryptoy_mem.h>
#include <cryptoy_file.h>

int64_t cryptoy_file_xor(FILE *a, FILE *b, FILE *r)
  { int N = 1024; /* Buffer size */
    unsigned char abuf[N]; /* Buffer for input file {a} and output file {r}. */
    uint64_t nw = 0;       /* Total bytes written to file {r}. */
    unsigned char bbuf[N]; /* Buffer for input file {b}. */

    while(1)
      { /* Read a buffer-load from each file: */
        size_t ma = fread(abuf, 1, N, a);
        if ((ma < N) && (! feof(a))) { fprintf(stderr, "i/o error reading file {a}\n"); return -1; }
        size_t mb = fread(bbuf, 1, N, b);
        if ((mb < N) && (! feof(b))) { fprintf(stderr, "i/o error reading file {b}\n"); return -1; }
        if (ma > mb) { fprintf(stderr, "file {a} is longer than file {b}\n"); return -1; }
        /* Combine buffers: */
        size_t mab = ma;
        cryptoy_mem_xor(mab, abuf, bbuf, abuf);
        size_t mw = fwrite(abuf, 1, mab, r);
        nw += mw;
        if (mw < mab) { fprintf(stderr, "i/o error writing the output file\n"); return -1; } 
        if (mab < N) { break; }
      }
    fflush(r);
    return nw;
  }

int64_t cryptoy_file_mix(FILE *a, FILE *b, FILE *x, FILE*r, FILE *s)
  { int N = 1024; /* Buffer size */
    unsigned char abuf[N]; /* Buffer for input file {a} and output file {r}. */
    unsigned char bbuf[N]; /* Buffer for input file {b} and output file {s}. */
    unsigned char xbuf[N]; /* Buffer for key file {x}. */
    uint64_t nw = 0;       /* Total bytes written to each output file. */

    while(1)
      { /* Read a buffer-load from each file: */
        size_t ma = fread(abuf, 1, N, a);
        if ((ma < N) && (! feof(a))) { fprintf(stderr, "i/o error reading file {a}\n"); return -1; }
        size_t mb = fread(bbuf, 1, N, b);
        if ((mb < N) && (! feof(b))) { fprintf(stderr, "i/o error reading file {b}\n"); return -1; }
        size_t mx = fread(xbuf, 1, N, x);
        if ((mx < N) && (! feof(x))) { fprintf(stderr, "i/o error reading file {x}\n"); return -1; }
        if (ma != mb) { fprintf(stderr, "files {a} and {b} have different lengths\n"); return -1; }
        if (ma > mx) { fprintf(stderr, "files {a} and {b} are longer than file {x}\n"); return -1; }
        /* Combine buffers: */
        size_t mab = ma;
        cryptoy_mem_mix(mab, abuf, bbuf, xbuf, abuf, bbuf);
        size_t mr = fwrite(abuf, 1, mab, r);
        if (mr < mab) { fprintf(stderr, "i/o error writing the output file {r}\n"); return -1; } 
        size_t ms = fwrite(bbuf, 1, mab, s);
        if (ms < mab) { fprintf(stderr, "i/o error writing the output file {s}\n"); return -1; } 
        nw += mr;
        if (mab < N) { /* Must be EOF: */ break; }
      }
    fflush(r);
    return nw;
  }

#define cryptoy_C_rights \
  "This file can be used and modified for any purpose as long as the source" \
  " is made available with any compiled code and the copyright and this" \
  " 'rights' note are preserved in all copies and derived versions."

#define cryptoy_C_warranty \
  "This file is provided 'as is'. The author and his employer are" \
  " not liable for any damages that may result from its use."
