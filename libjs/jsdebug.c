/* See jsstring.h */
/* Last edited on 2024-06-28 02:09:21 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <jsdebug.h>

void jsdebug_addr_span(char *name, void *ini, void *fin, int32_t n)
  { uint64_t aini = (uint64_t)ini;
    uint64_t afin = (uint64_t)fin;
    uint64_t tlen = afin - aini;
    uint64_t elen = tlen/n;
    fprintf(stderr, "  %-6s from %20lu to %20lu (%lu = %d*%lu bytes)\n", name, aini, afin,tlen, n, elen);
    if ((tlen % n) != 0) 
      { fprintf(stderr, "  ** oops - the length is not multiple of {n}\n");  }
  }

