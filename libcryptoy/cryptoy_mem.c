/* See {cryptoy_mem.h} */
/* Last edited on 2013-10-31 01:43:25 by stolfilocal */

#define cryptoy_C_copyright \
  "Copyright © 2013 Jorge Stolfi, State University of Campinas"

#define _GNU_SOURCE
#include <stdio.h>
#include <errno.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <cryptoy_mem.h>

void cryptoy_mem_xor(size_t n, byte_t a[], byte_t b[], byte_t r[])
  { int k;
    for (k = 0; k < n; k++) { r[k] = a[k] ^ b[k]; }
  }

void cryptoy_mem_mix(size_t n, byte_t a[], byte_t b[], byte_t x[], byte_t s[], byte_t r[])
  { int k;
    for (k = 0; k < n; k++) 
      { byte_t ak = a[k];
        byte_t bk = b[k];
        byte_t xk = x[k];
        r[k] = (byte_t)((ak & (~ xk)) | (bk & xk));
        s[k] = (byte_t)((ak & xk) | (bk & (~ xk)));
      }
  }

#define cryptoy_mem_C_rights \
  "This file can be used and modified for any purpose as long as the source" \
  " is made available with any compiled code and the copyright and this" \
  " 'rights' note are preserved in all copies and derived versions."

#define cryptoy_mem_C_warranty \
  "This file is provided 'as is'. The author and his employer are" \
  " not liable for any damages that may result from its use."
