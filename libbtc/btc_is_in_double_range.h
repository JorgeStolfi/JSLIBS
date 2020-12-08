#ifndef btc_is_in_double_range_H
#define btc_is_in_double_range_H

/* Checking whether a continuous BTC price bubble parameter is within its range. */
/* Last edited on 2015-04-20 01:07:18 by stolfilocal */
    
#include <bool.h>

bool_t btc_is_in_double_range(double vlo, double v, double vhi, int ib, char* vname, bool_t die);
  /* Checks whether {v} is in the interval {[vlo _ vhi]}.
    If it is, returns {THUE}; if not, returns {FALSE} when {die} is {FALSE},
    fails if {die} is TRUE. In any case, checks that {vlo <= vhi}, fails if not.
    The bubble index {ib} and the field name {*vname} are printed if error. */

#endif
