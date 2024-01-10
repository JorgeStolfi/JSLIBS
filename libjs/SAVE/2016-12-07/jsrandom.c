/* See jsrandom.h */
/* Last edited on 2013-10-25 01:30:35 by stolfilocal */

#define _GNU_SOURCE
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#include <affirm.h>

#include <jsrandom.h>

int32_t abrandom(int32_t a, int32_t b)
  {
    demand(a <= b, "bad interval");
    if (a == b) { return a; }
    uint32_t W = (uint32_t)(b - a);
    uint32_t RMAX = (uint32_t)RAND_MAX;
    demand (W <= RMAX, "interval too big");
    W++;
    /* Partition {0..RMAX} into {W} blocks of size {D} plus reminder: */
    uint32_t D = (RMAX+1)/W;
    uint32_t SLIM = D*W;
    /* Generate a random integer in {0..D*W-1}: */
    uint32_t s;
    do { s = (uint32_t)random(); } while (s >= SLIM); 
    return a + (int32_t)(s / D);
  }

/* The implementations of {frandom} and {drandom} below
  are not entirely correct, because normalization may 
  produce trailing zero bits with probability slightly 
  larger than 1/2. */
  
float frandom(void)
  {
    double rnd = (double)(random() & 8388607);
    return (float)(rnd/8388608.0);
  }

double drandom(void)
  {
    double d = 0.0;
    double rnd1 = (double)(random() & 536870911);
    d = (d + rnd1)/536870912.0;
    double rnd2 = (double)(random() & 8388607);
    d = (d + rnd2)/8388608.0;
    return (d);
  }

float fgaussrand(void)
  {
    return (float)(dgaussrand());
  }

static double js_min_gauss_R = 0.0;

double dgaussrand(void)
  { /* Polar method [Knuth II:3.4.1, Algorithm P] */
    if (js_min_gauss_R == 0.0)
      { js_min_gauss_R = sqrt(log(DBL_MAX))/DBL_MAX; }
    while(1)
      { double u = 2.0*drandom() - 1.0;
        double v = 2.0*drandom() - 1.0;
        double r = hypot(u, v);
        if ((r < 1.0) && (r > js_min_gauss_R))
          { /* In unit circle, and far enough from overflow: */
            return (M_SQRT2 * (u + v) * sqrt(-log(r))/r);
          }
      }
  }

#define MASK_22 ((1 << 22)-1)
#define MASK_20 ((1 << 22)-1)
  /* Mask that selects the lowest 22 bits, i.e. {2^22-1}. */

uint64_t uint64_random(void)
  { affirm(RAND_MAX >= MASK_22, "range of {random} is insufficient");
    affirm((RAND_MAX & MASK_22) == MASK_22, "range of {random} is not power of 2");
    uint64_t a = (random() & MASK_20);
    uint64_t b = (random() & MASK_22);
    uint64_t c = (random() & MASK_22);
    return (a << 44) | (b << 22) | c; 
  }

int64_t int64_random(void)
  { uint64_t x = uint64_random();
    return (int64_t)x;
  }

void randchars(char s[], int ns, char a[], int na)
  { int i;
    for (i = 0; i < ns; i++)
      { int k = abrandom(0, na-1);
        s[i] = a[k];
      } 
  }

