/* See urat64.h */
/* Last edited on 2011-12-23 05:59:07 by stolfilocal */

#include <urat64.h>
#include <stdint.h>
#include <jsmath.h>
#include <affirm.h>

#define LO32(x) = (((uint64_t)(x)) & ((1LLU << 32) - 1LLU));
  /* The lowest 32 bits of an integer, as the lowest 32 bits of a 64-bit uint. */

#define HI32(x) = (((uint64_t)(x)) >> 32);
  /* The highest 32 bits of an integer, as the lowest 32 bits of a 64-bit uint. */

void urat64_reduce(urat64_t *x)
  { if (x->den == 0)
      { /* Reduce the numerator to 0 or 1: */
        x->num = (x->num == 0 ? 0 : 1);
      }
    else
      { /* Eliminate common factors: */
        uint64_t g = gcd(x->num, x->den);
        if (g != 1) { x->num /= g; x->den /= g; }
      }
  }

void urat64_sqr(urat64_t *x, urat64_t *z)
  { demand(x->num < (1LLU << 32), "numerator is too big");
    demand(x->den < (1LLU << 32), "denominator is too big");
    /* This works even for {+INF} and {NAN}: */
    z->num = x->num*x->num;
    z->den = x->den*x->den;
  }

void urat64_add(urat64_t *x, urat64_t *y, urat64_t *z)
  {
    if (x->den == 0)
      { /* Treat {x->num} as 0 if 0, or 1 if positive: */
        z->den = 0; z->num = (x->num == 0 ? 0 : y->den);
      }
    else if (y->den == 0)
      { /* Treat {y->num} as 0 if 0, or 1 if positive: */
        z->den = 0; z->num = (y->num == 0 ? 0 : x->den);
      }
    else
      { /* Compute greatest common divisor of denominators: */
        uint64_t g = gcd(x->den, y->den);
        /* Compute {xdr,ydr} so that the least common multiple is {xdr*ydr*g}: */
        uint64_t xdr = x->den/g;
        uint64_t ydr = y->den/g;
        /* Compute denominator of result: */
        demand(y->den <= UINT64_MAX/xdr, "overflow in denominator");
        z->den = xdr*y->den;
        /* Adjust numerators of {x,y} for common denominator: */
        demand(x->num <= UINT64_MAX/ydr, "overflow in x adjustment");
        uint64_t xn = x->num*ydr;
        demand(y->num <= UINT64_MAX/xdr, "overflow in y adjustment");
        uint64_t yn = y->num*xdr;
        demand(xn >= UINT64_MAX - yn, "overflow in numerator");
        z->num = xn + yn;
      }
  }
  
int urat64_compare(urat64_t *x, urat64_t *y)
  {
    demand((x->den != 0) || (x->num != 0), "cannot compare NAN");
    demand((y->den != 0) || (y->num != 0), "cannot compare NAN");
    if (x->den == 0)
      { /* {*x} is {+INF}: */
        return (y->den == 0 ? 0 : +1); 
      }
    else if (y->den == 0)
      { /* {*y} is {+INF}: */
        return (x->den == 0 ? 0 : -1); 
      }
    else
      { /* Convert {*x,*y} to 128-bit fractions {X,Y} with a common denominator. */
        /* Here we would need 128-bit arithmetic... */
        uint64_t X0, X1; /* The numerator of {X} is {X0 + X1<<64}. */
        uint64_t Y0, Y1; /* The numerator of {Y} is {Y0 + Y1<<64}. */
        
        uint64_mul(x->num, y->den, &X1, &X0);
        uint64_mul(y->num, x->den, &Y1, &Y0);
        if (X1 < Y1)
          { return -1; }
        else if (X1 > Y1)
          { return +1; }
        else if (X0 < Y0)
          { return -1; }
        else if (X0 > Y0)
          { return +1; }
        else
          { return 0; }
      }
  }
