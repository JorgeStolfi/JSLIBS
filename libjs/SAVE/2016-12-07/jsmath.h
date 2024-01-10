#ifndef jsmath_H
#define jsmath_H

#include <stdint.h>

/* J. Stolfi's miscellaneous general-purpose definitions */
/* Last edited on 2016-04-18 19:38:20 by stolfilocal */

#ifndef INF
#define INF INFINITY
#endif

int64_t ipow(int64_t x, unsigned int y);
  /* Returns {x^y}. Does not check for overflow. */

int64_t imod(int64_t x, int64_t y);
  /* Mathematical remainder of {x ÷ y}, always in {0..y-1}.
    Fails with divide-by-zero if {y} is zero or negative. */

int64_t ifloor(int64_t x, int64_t y);
  /* Rounds {x} down to the nearest multiple of {y}.
    Fails with divide-by-zero if {y} is zero or negative.
    Beware of overflow. */

int64_t iceil(int64_t x, int64_t y);
  /* Rounds {x} up to the nearest multiple of {y}.
    Fails with divide-by-zero if {y} is zero or negative.
    Beware of overflow. */

int64_t iround(double x, double eps, int64_t d, int64_t r, int64_t imax);
  /* Returns the fraction {x/eps} rounded to an integer. 
  
    If {r < 0}, or {d <= 1}, rounds to the nearest integer. 
    In case of ties, rounds {x/eps} to the nearest even number.
    
    If {r >= 0} and {d >= 2}, then {r} must be in {0..d-1}.
    Rounds to the nearest integer that congruent to {r} modulo {d}.
    In case of ties, rounds to the nearest integer of the form 
    {2*d*k + r}. 
    
    Fails if the ratio {x/eps} is too large for the rouding to be
    accurate. If {imax} is not negative, also fails if the result 
    is greater than {imax} in absolute value. */

uint64_t gcd(uint64_t a, uint64_t b);
  /* Greatest common divisor of {a} and {b}.
    If both arguments are 0, returns 0. */
  
uint64_t lcm(uint64_t a, uint64_t b);
  /* Least common multiple of {a} and {b}.
    If either argument is 0, returns 0. */

uint64_t comb(int64_t n, int64_t k);
  /* Number of subset of {k} elements in a set with {n}
    elements. Not efficient - number of operations
    is at least {imin(k,n-k)}. Will fail if {n} or {k} is too large,
    even if the result fits in {uint64_t}. Returns 0
    if {k < 0} or {k > n}. */ 

int64_t imin(int64_t x, int64_t y);
int64_t imax(int64_t x, int64_t y);
  /* Returns maximum/minimum of {x,y}. */

uint64_t umin(uint64_t x, uint64_t y);
uint64_t umax(uint64_t x, uint64_t y);
  /* Returns maximum/minimum of {x,y}. */
  
void uint64_mul(uint64_t x, uint64_t y, uint64_t *Z1, uint64_t *Z0);
  /* Multiplies two 64-bit unsigned integers {x,y} to give a
    128-bit unsiged int.  The result is stored into {Z1,Z0}
    and should be interpreted as the integer {Z1*2^64 + Z0}. */

void int64_mul(int64_t x, int64_t y, int64_t *Z1, uint64_t *Z0);
  /* Multiplies two 64-bit signed integers {x,y} to give a
    128-bit siged int.  The result is stored into {Z1,Z0}
    and should be interpreted as the integer {Z1*2^64 + Z0}. */

unsigned int digits(uint64_t x);
  /* Number of decimal digits needed to write {x}. 
    Namely, {1} if {x==0}, else {floor(log10(x)) + 1} */

double falpow(double x, int n);
  /* Returns {x} to the power {n} falling, that is, that is, 
     {x*(x-1)*...*(x-n+1)}.  In particular, {falpow(n,n)} is {n!}. */

double rel_diff(double x, double y);
  /* Returns the relative difference between {x} and {y}, namely
    {(x-y)/sqrt(x^2 + y^2)/2}. This number is zero when {x == y}, is
    0.5 when only one of them is zero, and has extremal values -1
    and +1 when {x} and {y} are equal and opposite).  It is {NaN}
    when both {x} and {y} are zero. */

double abs_rel_diff(double x, double y, double abs_tol, double rel_tol);
  /* Returns the difference {d = x - y} divided by 
    {D = hypot(abs_tol, rel_tol*hypot(x,y)/sqrt(2))}.
    
    If {fabs(x),fabs(y)} are both very small, the denominator {D} is
    approximately {abs_tol}. In that case, if the result {d/D} is 1.0 in
    absolute value, then {fabs(d)} is approximately {abs_tol}.
    
    If either of {fabs(x),fabs(y)} is large compared to {abs_tol}, {D}
    is between {0.7*M} and {M}, where {M = max{fabs(x),fabs(y)}};
    thus, if the result is less that 1.0 in absolute value, then
    {fabs(d) < rel_tol*fabs(x)} and {fabs(d) < rel_tol*fabs(y)}; if
    the result is greater than 1 in absolute value, then {fabs(d) >
    0.7*rel_tol*fabs(x)} and {fabs(d) > 0.7*rel_tol*fabs(y)}. */

void expand_range(double *zP, double zlo, double zhi, double *dzP);
  /* Applies to {*zP} a rational quadratic map that expands the interval
    {[zlo_zhi]} to the whole {R^2}.  In particular, if 
    {z==zlo} the result is {-INF}, if {z==zhi} the result is {+INF}.
    The result is {NAN} if {*zP} is outside that interval,
    or if {zlo >= zhi}, or if {*zP} is {NAN}.
    
    If {dzP} is not null, also stores into {*dzP} the derivative
    of the map at {*zP}. */

void contract_range(double *zP, double zlo, double zhi, double *dzP);
  /* Applies to {*zP} a map that contracts the whole {R^2} to {[zlo_zhi]}.
    The map is the inverse of the {expand_range}.  In particular, if 
    {z==-INF} the result is {zlo}, if {z==+INF} the result is {zhi}.
    The result is {NAN} if {zlo >= zhi} or if {*zP} is {NAN}.
    
    If {dzP} is not null, also stores into {*dzP} the derivative
    of the map at {*zP}. */


#endif
