#ifndef jsrandom_H
#define jsrandom_H

/* Alternative random generator functions */
/* Last edited on 2024-11-15 19:14:32 by stolfi */

#include <stdint.h>
#include <stdlib.h> 

/* 
  These procedure use {random} from {stdlib.h} that is supposed to
  return a random {long int} result in {0..RAND_MAX}, where {RAND_MAX}
  is hoped to be at least {2^22-1}. To get a repeatable or truly random
  sequence of values, call {srandom} from that library with a suitable
  {unsigned int} seed value.
*/

int32_t int32_random(void);
int64_t int64_random(void);
  /* Returns a random 32-bit or 64-bit signed integer (all bits random). 
    Note that the result may be negative. */

uint32_t uint32_random(void);
uint64_t uint64_random(void);
  /* Returns a random 32-bit or 64-bit unsigned integer (all bits random). */
  
uint32_t uint32_mrandom(uint32_t max);  
uint64_t uint64_mrandom(uint64_t max); 
  /* Returns a random integer with uniform distr in {[0..max]}. */

int32_t int32_abrandom(int32_t a, int32_t b);
int64_t int64_abrandom(int64_t a, int64_t b);
uint32_t uint32_abrandom(uint32_t a, uint32_t b);
uint64_t uint64_abrandom(uint64_t a, uint64_t b);
  /* Returns a random integer with uniform distr in {[a..b]}. */
  
int32_t abrandom(int32_t a, int32_t b);
  /* Same as {int32_abrandom} for compatibility (deprecated). */

float frandom(void);
double drandom(void);
  /* Returns a random number with uniform distr in [0.0 __ 1.0). */

float fabrandom(float a, float b);
double dabrandom(double a, double b);
  /* Returns a random number with uniform distr in {[a __ b)}. */

float fgaussrand(void);
double dgaussrand(void);
  /* Returns a random number with normal distribution (Gaussian distr
    with zero mean and unit deviation). */

float floggaussand(float avg, float dev);
double dloggaussrand(double avg, double dev);
  /* Returns a random number with log-Gaussian distribution
    of mean {avg} and deviation {dev}.
    
    If {dev} is zero, the result will be {avg}. Otherwise,
    the natural log of the result will have a Gaussian distribution
    with suitable mean and deviation. */

void randchars(char s[], uint32_t ns, char a[], uint32_t na);
  /* Sets {s[0..ns-1]} to a string of letters picked from {a[0..na-1]},
    with uniform and independent random probabilities.
    NOTE: does not add the final 0. */
    
uint64_t *uint64_choose(uint64_t n, size_t k, uint64_t *perm);
  /* Chooses {k} distinct integers in the range {0..n-1}.  Fails
    if {k > n}.
    
    If {perm} is not null, it must be the address of an array with at
    least {k} elements. If {perm} is null, the procedure allocate an
    array of {k} elements from the heap. Either way, procedure will
    store the chosen numbers in elements {0..k-1} if that array, sorted
    in increasing order, and return its address.
    
    The current implementation takes time proportional to {k^2} in the
    worst and average case. */


#endif
