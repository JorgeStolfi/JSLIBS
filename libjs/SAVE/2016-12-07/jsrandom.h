#ifndef jsrandom_H
#define jsrandom_H

/* Alternative random generator functions */
/* Last edited on 2007-01-03 11:42:40 by stolfi */

#include <stdint.h>

#if (defined(SunOS4) || defined(SunOS5))
extern long random(void);
extern int srandom(unsigned seed);
#endif

int32_t abrandom(int32_t a, int32_t b);
  /* Returns a random integer with uniform distr in {[a..b]}. */

float frandom(void);
double drandom(void);
  /* Returns a random number with uniform distr in [0.0 __ 1.0). */

float fgaussrand(void);
double dgaussrand(void);
  /* Returns a random number with Gaussian distr, zero mean, 
    unit deviation. */

uint64_t uint64_random(void);
  /* Returns a random 64-bit unsigned integer (all bits random). */

int64_t int64_random(void);
  /* Returns a random 64-bit signed integer (all bits random).
    Note that the result may be negative. */

void randchars(char s[], int ns, char a[], int na);
  /* Sets {s[0..ns-1]} to a string of letters picked from {a[0..na-1]},
    with uniform and independent random probabilities.
    NOTE: does not add the final 0. */

#endif
