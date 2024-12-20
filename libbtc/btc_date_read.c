/* See {btc_date_read.h} */
/* Last edited on 2024-12-05 10:23:28 by stolfi */

#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <fget.h>

#include <btc_date_lookup.h>

#include <btc_date_read.h>

int btc_date_read(FILE* rd, char* name, int nd, char* dt[], bool_t debug)
  {
    char* dts = fget_string(rd); /* Date of start of relevant period */
    if (debug) { fprintf(stderr, " %s", dts); }
    int id = btc_date_lookup(nd, dt, dts);
    free(dts);
    return id;
  }
  
