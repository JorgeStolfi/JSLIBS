#ifndef btc_date_lookup_H
#define btc_date_lookup_H

/* Finding a date in a sorted list of dates. */
/* Last edited on 2015-04-20 00:36:24 by stolfilocal */

int btc_date_lookup(int nd, char *dt[], char* dtx);
  /* Returns an integer {id} in {0..nd-1} such that {dt[id]} is equal (as string) to {dtx}. Assumes
    the dates {dt[0..nd-1]} are in increasing lexicographic order (same as chronological order
    for ISO format). */

#endif
