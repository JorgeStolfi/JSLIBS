/* See jsstring.h */
/* Last edited on 2024-11-16 06:34:43 by stolfi */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <jsdebug.h>

void jsdebug_addr_span(char *name, void *ini, void *fin, uint32_t n)
  { demand(n > 0, "invalid {n}");
    demand(ini != fin, "invalid {ini/fin} arguments");
    uint64_t aini = (uint64_t)ini;
    uint64_t afin = (uint64_t)fin;
    uint64_t tlen = (uint64_t)(aini < afin ? afin - aini : aini - afin);
    uint64_t elen = tlen/n;
    fprintf(stderr, "  %-6s from %20lu to %20lu (%lu = %d*%lu bytes)\n", name, aini, afin,tlen, n, elen);
    if ((tlen % n) != 0) 
      { fprintf(stderr, "  ** oops - the length is not multiple of {n}\n");  }
  }

