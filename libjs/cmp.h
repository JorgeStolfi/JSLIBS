#ifndef cmp_H
#define cmp_H

/* Compare values returning a {sign_t} result. */
/* Last edited on 2007-10-28 21:31:36 by stolfi */

#include <sign.h>
#include <stdint.h>

sign_t cmp_int(const int* x, const int *y);
sign_t cmp_long_int(const long int *x, const long int *y);
sign_t cmp_unsigned_int(const unsigned int* x, const unsigned int *y);
sign_t cmp_unsigned_long_int(const unsigned long int *x, const unsigned long int *y);
sign_t cmp_int8(const int8_t *x, const int8_t *y);
sign_t cmp_int16(const int16_t *x, const int16_t *y);
sign_t cmp_int32(const int32_t *x, const int32_t *y);
sign_t cmp_int64(const int64_t *x, const int64_t *y);
sign_t cmp_uint8(const uint8_t *x, const uint8_t *y);
sign_t cmp_uint16(const uint16_t *x, const uint16_t *y);
sign_t cmp_uint32(const uint32_t *x, const uint32_t *y);
sign_t cmp_uint64(const uint64_t *x, const uint64_t *y);
sign_t cmp_float(const float *x, const float *y);
sign_t cmp_double(const double *x, const double *y);
  /* These procedures return {-1} if {(*x) < (*y)}, {+1} if {(*x) > (*y)}, 
    and {0} if {(*x) == (*y)}. Beware of {float} and {double} minus zero. */

#endif
