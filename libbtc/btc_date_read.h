#ifndef btc_date_read_H
#define btc_date_read_H

/* Reading ISO format dates. */
/* Last edited on 2024-12-05 10:23:30 by stolfi */

#include <stdio.h>

#include <bool.h>

int btc_date_read(FILE* rd, char* name, int nd, char* dt[], bool_t debug);
  /* Parses an ISO date from {rd}, skipping leading blanks if needed. Then 
    converts the date to an integer day index by looking it up in {dt[0..nd-1]}. */

#endif
