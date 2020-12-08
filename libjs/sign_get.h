#ifndef sign_get_H
#define sign_get_H

/* Extract sign of a value as a {sign_t} data type. */
/* Last edited on 2007-10-28 19:32:37 by stolfi */

#include <sign.h>
#include <stdint.h>

sign_t sign_int(int x);
sign_t sign_long_int(long int x);
sign_t sign_int32(int32_t x);
sign_t sign_int64(int64_t x);
sign_t sign_float(float x);
sign_t sign_double(double x);
  /* These procedures return {-1} if {x} is negative,
    {+1} if {x} is positive, and {0} if {x} is zero.
    Beware of {float} and {double} minus zero. */

#endif
