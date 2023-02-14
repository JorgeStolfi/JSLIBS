/* See {btc_price_series_read.h} */
/* Last edited on 2023-02-12 07:53:13 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <jsfile.h>

#include <btc_price_series_read.h>

void btc_price_series_read(FILE* rd, int* ndP, char*** dtP, double **apP)
  {
    bool_t debug = TRUE;
    
    if (debug) { fprintf(stderr, "reading price data from {stdin}\n"); }
    
    int nd = nget_int32(rd, "days"); fget_eol(rd); /* Number of days in series. */
    char** dt = notnull(malloc(nd*sizeof(char*)), "no mem");   /* Dates indexed {0..nd-1}. */
    double* ap = notnull(malloc(nd*sizeof(double)), "no mem"); /* Daily mean prices, indexed {0..nd-1}. */
    
    if (debug) { fprintf(stderr, "days = %d\n", nd); }
    
    int id;
    for (id = 0; id < nd; id++)
      { dt[id] = fget_string(rd);
        ap[id] = fget_double(rd);
        fget_eol(rd);
      }
    (*ndP) = nd;
    (*dtP) = dt;
    (*apP) = ap;
  }

