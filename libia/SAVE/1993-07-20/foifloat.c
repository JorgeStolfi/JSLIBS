/* See foifloat.h */

#include "foifloat.h"
#include "foimisc.h"
#include "iomisc.h"
#include <stdio.h>

char *ROUND_BUFP = NULL;

void flt_print (FILE *f, Float x)
  {
    fprintf(f, F_FMT_SPEC, x);
  }

void flt_add(Float x, Float y, FloatP zp, FloatP errp)
  {
    /* There should be a simpler, better way... */
    Float zhi, zlo, del;
    *zp = x + y;         /* May overflow. */
    ROUND_DOWN;
    zlo = x + y;
    ROUND_UP;
    zhi = x + y;
    if ((zhi >= PlusInfinity) || (zlo <= MinusInfinity))
      { *zp = PlusInfinity;
        *errp = PlusInfinity;
        return;
      }
    ROUND_UP;
    del = zhi - zlo;    /* Shouldn't overflow.  (What about underflow?) */
    *errp = *errp + del;  /* May overflow */
    if (*errp >= PlusInfinity) *errp = PlusInfinity ;
  }

void flt_sub(Float x, Float y, FloatP zp, FloatP errp)
  {
    /* There should be a simpler, better way... */
    Float zhi, zlo, del;
    *zp = x - y;         /* May overflow. */
    ROUND_DOWN;
    zlo = x - y;
    ROUND_UP;
    zhi = x - y;
    if ((zhi >= PlusInfinity) || (zlo <= MinusInfinity))
      { *zp = PlusInfinity;
        *errp = PlusInfinity;
        return;
      }
    ROUND_UP;
    del = zhi - zlo;    /* Shouldn't overflow.  (What about underflow?) */
    *errp = *errp + del;  /* May overflow */
    if (*errp >= PlusInfinity) *errp = PlusInfinity ;
  }

void flt_mul(Float x, Float y, FloatP zp, FloatP errp)
  {
    Float zhi, zlo, del;
    *zp = x * y;         /* May overflow. */
    ROUND_DOWN;
    zlo = x + y;
    ROUND_UP;
    zhi = x + y;
    if ((zhi >= PlusInfinity) || (zlo <= MinusInfinity))
      { *zp = PlusInfinity;
        *errp = PlusInfinity;
        return;
      }
    ROUND_UP;
    del = zhi - zlo;    /* Shouldn't overflow.  (What about underflow?) */
    *errp = *errp + del;  /* May overflow */
    if (*errp >= PlusInfinity) *errp = PlusInfinity ;
  }

void flt_inv(Float x, FloatP zp, FloatP errp)
  {
    Float zhi, zlo, del;
    if ((x >= PlusInfinity))
      { *zp = PlusInfinity;
        *errp = PlusInfinity;
        return;
      }
    if (x == Zero)
      { error ("flt_inv: argument is zero"); }
    *zp = One/x;
    ROUND_DOWN;
    zlo = One/x;
    ROUND_UP;
    zhi = One/x;
    if ((zhi >= PlusInfinity) || (zlo <= MinusInfinity))
      { *zp = PlusInfinity;
        *errp = PlusInfinity;
        return;
      }
    ROUND_UP;
    del = zhi - zlo;    /* Shouldn't overflow.  (What about underflow?) */
    *errp = *errp + del;  /* May overflow */
    if (*errp >= PlusInfinity) *errp = PlusInfinity ;
  }

void flt_sqrt(Float x, FloatP zp, FloatP errp)
  {
    Float zhi, zlo, del;
    if ((x >= PlusInfinity))
      { *zp = PlusInfinity;
        *errp = PlusInfinity;
        return;
      }
    if (x < Zero)
      { error ("flt_sqrt: argument is negative"); }
    if (x == Zero)
      { *zp = Zero;
        return;
      }
    *zp = sqrt(x);
    ROUND_DOWN;
    zlo = sqrt(x);
    ROUND_UP;
    zhi = sqrt(x);
    ROUND_UP;
    del = zhi - zlo;    /* Shouldn't overflow.  (What about underflow?) */
    *errp = *errp + del;  /* May overflow */
    if (*errp >= PlusInfinity) *errp = PlusInfinity ;
  }

Float flt_random(void)
  {
    Float hi = ((random()&65535) + 0.0)/65536.0;
    Float lo = ((random()&65535) + 0.0)/65536.0;
    return(hi + lo/65536.0);
  }
