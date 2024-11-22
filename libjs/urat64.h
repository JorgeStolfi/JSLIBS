#ifndef urat64_H
#define urat64_H

/* Unsigned rational fractions with 64-bit numerators and denominators. */
/* Last edited on 2024-11-15 19:16:48 by stolfi */

#include <stdint.h>

typedef struct urat64_t { uint64_t num, den; } urat64_t;
  /* A rational number {num/den}. The two numbers need not be
    relatively prime. If {den} is 0, the value is {+INF} if {num > 0},
    {NAN} if {num == 0}. */

#define urat64_ZERO  ((urat64_t){ 0, 1 })
#define urat64_ONE   ((urat64_t){ 1, 1 })
#define urat64_INF   ((urat64_t){ 1, 0 })
#define urat64_NAN   ((urat64_t){ 0, 0 })
  /* Some useful {urat64_t} constants. */

void urat64_reduce(urat64_t *x);
  /* Reduces {*x} to its lowest terms (with {num} and {den} relatively prime).  */

void urat64_sqr(urat64_t *x, urat64_t *z);
  /* Sets {*z} to {(*x)^2}.  The result is reduced iff {*x} is.  */

void urat64_add(urat64_t *x, urat64_t *y, urat64_t *z);
  /* Sets {z} to {x+y}.  The result is reduced iff {x,y} are.  */

int32_t urat64_compare(urat64_t *x, urat64_t *y);
  /* Returns {-1,00,+1} depending on whether {*x} is less than,
    equal to, or greater than {*y}, respectively.  Fails
    if either argument is {NAN}.  */

#endif
