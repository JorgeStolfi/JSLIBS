#ifndef btc_bubble_parms_validate_H
#define btc_bubble_parms_validate_H

/* Validating BTC price bubble parameters. */
/* Last edited on 2015-04-20 01:03:07 by stolfilocal */

#include <bool.h>

#include <btc_bubble_t.h>

bool_t btc_bubble_parms_validate(char* fName, int nlin, int nd, char* dt[], btc_bubble_t* bpj);
  /* Checks the validity of the parameters in bubble {*bpj}. In case of
    problems, writes warning message(s) to {stderr} and returns FALSE.
    Otherwise returns TRUE. Uses the dates {dt[0..nd-1]} for error
    messages. */

#endif
