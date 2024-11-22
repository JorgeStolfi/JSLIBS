/* See jsrandom.h */
/* Last edited on 2024-11-22 02:21:18 by stolfi */

#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>

#include <jsrandom.h>

/* SIMPLE INTEGER RANDOM FUNCTIONS */

#define MASK_22 ((1 << 22)-1)
  /* Mask that selects the lowest 22 bits, i.e. {2^22-1}. */

#define MASK_20 ((1 << 20)-1)
  /* Mask that selects the lowest 20 bits, i.e. {2^20-1}. */

#define MASK_12 ((1 << 12)-1)
  /* Mask that selects the lowest 12 bits, i.e. {2^12-1}. */
  
#define U32RMAX ((uint32_t)RAND_MAX)
  /* {RAND_MAX} as a 32-bit unsigned integer. */

uint32_t uint32_random(void)
  { affirm(U32RMAX >= MASK_22, "range of {random} is insufficient");
    affirm(((U32RMAX + 1) & U32RMAX) == 0, "range of {random} is not power of 2");
    uint32_t a = (random() & MASK_12);
    uint32_t b = (random() & MASK_22);
    return (a << 22) | b; 
  }

int32_t int32_random(void)
  { uint32_t x = uint32_random();
    return (int32_t)x;
  }

uint64_t uint64_random(void)
  { affirm(U32RMAX >= MASK_22, "range of {random} is insufficient");
    affirm(((U32RMAX + 1) & U32RMAX) == 0, "range of {random} is not power of 2");
    uint64_t a = (random() & MASK_20);
    uint64_t b = (random() & MASK_22);
    uint64_t c = (random() & MASK_22);
    return (a << 44) | (b << 22) | c; 
  }

int64_t int64_random(void)
  { uint64_t x = uint64_random();
    return (int64_t)x;
  }

/* RANDOM INTEGERS WITH GIVEN MAXIMUM: */

uint32_t uint32_mrandom(uint32_t max)
  {
    /* Some simple cases: */
    if (max == 0) { return 0; }
    if (max == (uint32_t)RAND_MAX) { return (uint32_t)random(); }
    /* If {max} is {2^k-1}, just take the lower bits: */
    if ((max & (max+1)) == 0) { return (uint32_random() & max); }
    
    /* We need to generate random in the proper range and take mod: */
    assert(max < UINT32_MAX);
    /* Split the range {0..UINT32_MAX} into {D+1} pieces of size {max+1} plus some remainder: */
    uint32_t W = max + 1;
    assert(W >= 2);
    uint32_t D = ((uint32_t)(0 - W))/W; /* Result should be OK even if {W=2}. */
    /* We need a random number in {0..SMAX} where {SMAX = (D+1)*W-1}: */
    uint32_t SMAX = D*W + max; /* Should not overflow. */
    if ((max > SMAX) || (UINT32_MAX - SMAX > SMAX))
      { fprintf(stderr, "SMAX = %u  max = %u\n", SMAX, max); }
    assert(max <= SMAX);
    assert(UINT32_MAX - SMAX <= SMAX);
    uint32_t s;
    do { s = uint32_random(); } while (s > SMAX); 
    return s % W;
  }

uint64_t uint64_mrandom(uint64_t max)
  {
    /* Some simple cases: */
    if (max == 0) { return 0; }
    if (max == (uint64_t)RAND_MAX) { return (uint64_t)random(); }
    if (max == (uint64_t)UINT32_MAX) { return (uint64_t)uint32_random(); }
    if (max == UINT64_MAX) { return uint64_random(); }
    /* If {max} is small enough, use the 32-bit version: */
    if (max <= (uint64_t)UINT32_MAX) { return (uint64_t)uint32_mrandom((uint32_t)max); }
    /* If {max} is {2^k-1}, just take the lower bits: */
    if ((max & (max+1)) == 0) { return (uint64_random() & max); }
    
    /* We need to generate random in the proper range and take mod: */
    assert(max < UINT64_MAX);
    /* Split the range {0..UINT64_MAX} into {D+1} pieces of size {max+1} plus some remainder: */
    uint64_t W = max + 1;
    assert(W >= 2);
    uint64_t D = ((uint64_t)(0 - W))/W; /* Result should be OK even if {W=2}. */
    /* We need a random number in {0..SMAX} where {SMAX = (D+1)*W-1}: */
    uint64_t SMAX = D*W + max; /* Should not overflow. */
    if ((max > SMAX) || (UINT64_MAX - SMAX > SMAX))
      { fprintf(stderr, "SMAX = %lu  max = %lu\n", SMAX, max); }
    assert(max <= SMAX);
    assert(UINT64_MAX - SMAX <= SMAX);
    uint64_t s;
    do { s = uint64_random(); } while (s > SMAX); 
    return s % W;
  }

/* RANDOM INTEGERS IN GIVEN RANGE: */

uint32_t uint32_abrandom(uint32_t a, uint32_t b)
  {
    demand(a <= b, "bad interval");
    if (a == b) { return a; }
    uint32_t rmax = (uint32_t)(b - a);
    uint32_t r = uint32_mrandom(rmax);
    return a + r;
  }

int32_t int32_abrandom(int32_t a, int32_t b)
  {
    demand(a <= b, "invalid interval");
    if (a == b) { return a; }
    /* Convert {a,b} to {uint32_t} by adding the zero shift: */
    uint32_t zshift = ((uint32_t)1) << 31; 
    uint32_t ua = ((uint32_t)a) + zshift;
    uint32_t ub = ((uint32_t)b) + zshift;
    /* Compute a random {uint32_t} in the range {ua..ub}: */
    uint32_t ur = uint32_abrandom(ua, ub);
    /* Shift back: */
    return (int32_t)(ur - zshift);
  }
 
int32_t abrandom(int32_t a, int32_t b)
  { return int32_abrandom(a, b); }

uint64_t uint64_abrandom(uint64_t a, uint64_t b)
  {
    demand(a <= b, "bad interval");
    if (a == b) { return a; }
    uint64_t rmax = (uint64_t)(b - a);
    uint64_t r = uint64_mrandom(rmax);
    return a + r;
  }

int64_t int64_abrandom(int64_t a, int64_t b)
  {
    demand(a <= b, "invalid interval");
    if (a == b) { return a; }
    /* Convert {a,b} to {uint64_t} by adding the zero shift: */
    uint64_t zshift = ((uint64_t)1) << 63; 
    uint64_t ua = ((uint64_t)a) + zshift;
    uint64_t ub = ((uint64_t)b) + zshift;
    /* Compute a random {uint64_t} in the range {ua..ub}: */
    uint64_t ur = uint64_abrandom(ua, ub);
    /* Shift back: */
    return (int64_t)(ur - zshift);
  }
  
/* RANDOM FLOATS */  

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

float fabrandom(float a, float b)
  { return (float)dabrandom((double)a, (double)b); } 
  
double dabrandom(double a, double b)
  { double d = drandom();
    /* This should never overflow: */
    return (1-d)*a + d*b;
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

float floggaussand(float avg, float dev)
  { 
    return (float)(dloggaussrand(avg, dev));
  }
  
double dloggaussrand(double avg, double dev)
  {
    if (dev == 0.0)
      { return avg; }
    else
      { demand(avg != 0.0, "invalid {avg}");
        double sig = (avg < 0.0 ? -1.0 : +1.0); /* Sign of result. */
        avg = fabs(avg); /* Sign will be set later. */
        dev = fabs(dev); /* Ignore sign of {dev}. */
        /* Convert {avg,dev} to the mean and deviation {avg_z,dev_z} of {z = ln(|x|)}: */
        double rdv = dev/avg;   /* Deviation relative to the mean. */
        double avg_z = log(avg/hypot(1.0, rdv));
        double dev_z = sqrt(log(1.0 + rdv*rdv));
        /* Generate the Gaussian variable {z}: */
        double z = avg_z + dev_z*dgaussrand();
        /* Map to the result {x}: */
        return sig*exp(z);
      }
  }

/* RANDOM CHARS: */

void randchars(char s[], uint32_t ns, char a[], uint32_t na)
  { for (int32_t i = 0; i < ns; i++)
      { uint32_t k = uint32_abrandom(0, na-1);
        s[i] = a[k];
      } 
  }

/* RANDOM SUBSETS */

uint64_t *uint64_choose(uint64_t n, size_t k, uint64_t *perm)
  {
    demand(k <= n, "invalid arguments");
    if ((k > 0) && (perm == NULL)) 
      { perm = notnull(malloc(k*sizeof(uint64_t)), "no mem"); }
    for (int64_t i = 0; i < k; i++)
      { int64_t pi = int64_abrandom(i, (int64_t)n-1);
        int64_t j = i;
        while ((j > 0) && (perm[j-1] >= pi))
          { perm[j] = perm[j-1];
            assert(pi > 0);
            pi = pi - 1;
            j = j - 1;
          }
        perm[j] = pi;
      }
    return perm;
  }
    
