#ifndef btc_is_in_int_range_H
#define btc_is_in_int_range_H

/* Checking whether an integer BTC price bubble parameter is within its range. */
/* Last edited on 2015-04-20 01:09:13 by stolfilocal */

#include <bool.h>

bool_t btc_is_in_int_range(int vlo, int v, int vhi, int ib, char* vname, bool_t die);
  /* Checks whether {v} is in the interval {vlo..vhi}.
    If it is, returns {THUE}; if not, returns {FALSE} when {die} is {FALSE},
    fails if {die} is TRUE. In any case, checks that {vlo <= vhi}, fails if not.
    The bubble index {ib} and the field name {*vname} are printed if error. */

#endif
