#ifndef btc_bubble_parms_write_H
#define btc_bubble_parms_write_H

/* Writing BTC price bubble parameters. */
/* Last edited on 2015-04-29 21:55:58 by stolfilocal */

#include <btc_bubble_t.h>

void btc_bubble_parms_write(char* outPrefix, char* tag, int nd, char* dt[], int nb, btc_bubble_t bp[]);
  /* Writes the parameters {bp[0..nb-1]}to file
    "{outPrefix}{tag}.parms", in the format described in {btc_bubble_parms_read_INFO}.
    (see {btc_bubble_parms_read.h}). The dates {dt[0..nd-1]} are used to convert
    the day indices in {bp} to ISO dates. */

#endif
