/* See jsmath.h */
/* Last edited on 2016-04-18 19:38:13 by stolfilocal */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>

#include <jsmath.h>
#include <affirm.h>

int64_t ipow(int64_t x, unsigned int y)
  { if (y == 0) { return 1; }
    if (x == 0) { return 0; }
    if (x == 1) { return 1; }
    if (y == 1) { return x; }
    /* Find largest power of 2 {b} not exceeding {y}: */ 
    unsigned int h = y/2;
    unsigned int b = 1;
    while (b <= h) { b += b; }
    /* Compute result by binary method: */ 
    int64_t p = x;
    y -= b; b /= 2;
    while (b > 0)
      { p = p*p; 
        if (y >= b) { p *= x; y -= b; } 
        b /= 2;
      }
    return p;
  }

int64_t imod(int64_t x, int64_t y)
  { /* If {y} is negative, force a divide-by-zero: */
    if (y < 0) { y = 0; }
    x -= (x/y)*y;
    if (x < 0) { x += y; }
    return x;
  }

int64_t ifloor(int64_t x, int64_t y)
  { /* If {y} is negative, force a divide-by-zero: */
    if (y < 0) { y = 0; }
    int64_t z = (x/y)*y;
    if (z > x) { z -= y; }
    return z;
  }

int64_t iceil(int64_t x, int64_t y)
  { /* If {y} is negative, force a divide-by-zero: */
    if (y < 0) { y = 0; }
    int64_t z = -((-x)/y)*y;
    if (z < x) { z += y; }
    return z;
  }
  
#define double_EXACT_MAX ((double)((((int64_t)1) << 53) - ((int64_t)1)))
 /* The max integer {n} such that {1..n} can be safely stored into a 
   IEEE {double} without any rounding or truncation.  The {-1UL} 
   is just for extra safety, to account for the round to even/odd
   option. */

int64_t iround(double x, double eps, int64_t d, int64_t r, int64_t imax)
  { 
    demand(eps > 0.0, "{eps} must be positive");
    demand(d <= 2*((int64_t)double_EXACT_MAX), "invalid divisor {d}");
    
    double div; /* 1 if {d <= 1}, {d} otherwise. */
    double rem; /* {r} converted to double. */
    if ((r >= 0) && (d >= 2))
      { demand((0 <= r) && (r < d), "invalid remainder {r}");
        div = (double)d; rem = (double)r;
      }
    else
      { div = 1.0; rem = 0.0; }
    
    /* Compute {qd = (x/eps - rem)/div}, the amount that should be rounded to integer: */
    double qd = (((double)x)/((double)eps) - rem)/div;
    demand(isinf(qd) == 0, "overflow in rounding");
    demand(fabs(qd) <= double_EXACT_MAX, "result is too large for accurate rounding");
    
    /* Round {qd} to an integer {qn}, still represented as double: */
    double qlo = floor(qd + 0.5); /* Rounds to nearest; half-integers up. */
    double qhi = ceil(qd - 0.5); /* Rounds to nearest; half-integers down. */
    double qn;
    if (qlo == qhi)
      { qn = qlo; }
    else 
      { /* Tie in rounding: round {qn} to even. */
        if (fmod(qlo,2.0) != 0)
          { qn = qlo; }
        else
          { qn = qhi; }
      }
    assert(qn == floor(qn));
    
    /* Compute the result, still as double: */
    double qr = div*qn + rem;
    
    /* Convert to integer: */
    int64_t qk = (int64_t)qr;
    assert(qr == (double)qk);

    if (imax >= 0)
      { /* Check user-defined range. The negation of {qk} can't oversflow: */
        demand((qk < 0 ? -qk : qk) <= imax, "result is too large");
      }

    return qk;
  }

uint64_t gcd(uint64_t a, uint64_t b)
  { while (b != 0)  { uint64_t r = a % b; a = b; b = r; }
    return a;
  }

uint64_t lcm(uint64_t a, uint64_t b)
  { uint64_t g = gcd(a, b);
    if (g == 0) { return 0; }
    /* The following math overflows only if overflow is unavoidable: */
    return ((uint64_t)a/g) * (uint64_t)b;
  }

uint64_t comb(int64_t n, int64_t k)
  { 
    uint64_t p = 1;
    if (k < 0) { return 0; }
    if (k > n) { return 0; }
    if (k > n-k) { k = n-k; }
    uint64_t a = n, b = 0;
    while(b < k)
      { b++;
        uint64_t q = p*a;
         /* Overflow check. */
        if (q/a != p) { assert(0); }
        p = q/b;
        a--;
      }
    return p;
  }

int64_t imin(int64_t x, int64_t y)
  { return (x < y ? x : y); }

int64_t imax(int64_t x, int64_t y) 
  { return (x > y ? x : y); }
  
uint64_t umin(uint64_t x, uint64_t y)
  { return (x < y ? x : y); }

uint64_t umax(uint64_t x, uint64_t y) 
  { return (x > y ? x : y); }
  
#define LO32(x) (((uint64_t)(x)) & ((1LLU << 32) - 1LLU))
  /* The lowest 32 bits of an integer, as the lowest 32 bits of a 64-bit uint. */

#define HI32(x) (((uint64_t)(x)) >> 32)
  /* The highest 32 bits of an integer, as the lowest 32 bits of a 64-bit uint. */

void uint64_mul(uint64_t x, uint64_t y, uint64_t *Z1, uint64_t *Z0)
  {
    uint64_t x0 = LO32(x), x1 = HI32(x);
    uint64_t y0 = LO32(y), y1 = HI32(y);
    
    uint64_t A00 = x0*y0;
    uint64_t A01 = x0*y1;
    uint64_t A10 = x1*y0;
    uint64_t A11 = x1*y1;
    
    (*Z0) = A00 + ((A01 + A10) << 32);
    uint64_t P = HI32(A00) + LO32(A10) + LO32(A01);
    (*Z1) = (P >> 32) + HI32(A10) + HI32(A01) + A11; /* Should not overflow */
  }

void int64_mul(int64_t x, int64_t y, int64_t *Z1, uint64_t *Z0)
  {
    uint64_mul((uint64_t)x, (uint64_t)y, (uint64_t *)Z1, Z0);
    if (x < 0) { (*Z1) -= y; }
    if (y < 0) { (*Z1) -= x; }
  }

unsigned int digits(uint64_t x)
  { unsigned int d = 1;
    while (x > 9) { x /= 10; d++; }
    return d;
  }

double falpow(double x, int n)
  { 
    /* !!! Improve speed and generalize to fractional {n} using {lngamma} !!! */
    double p = 1;
    while(n > 0) { p = p*x; x = x-1; n--; }
    return p;
  }

double rel_diff(double x, double y)
  { double mx = fabs(x), my = fabs(y);
    double m = (mx > my ? mx : my);
    if (m == 0.0) 
      { return 0.0; }
    else
      { double rx = x/m, ry = y/m;
        return 0.5*(rx - ry)/sqrt(rx*rx + ry*ry);
      }
  }
  
double abs_rel_diff(double x, double y, double abs_tol, double rel_tol)
  { double x2 = x*x;
    double y2 = y*y;
    double D2 = abs_tol*abs_tol + rel_tol*rel_tol*(x2+y2)/2 + 1.0e-300;
    return (x - y)/sqrt(D2);
  }

void expand_range(double *zP, double zlo, double zhi, double *dzP)
  { 
    double z = (*zP);
    if (isnan(z) || (z < zlo) || (z > zhi) || (zlo >= zhi))
      { (*zP) = NAN; if (dzP != NULL) { (*dzP) = NAN; } }
    else if (z == zlo)
      { (*zP) = -INF; if (dzP != NULL) { (*dzP) = NAN; } }
    else if (z == zhi)
      { (*zP) = +INF; if (dzP != NULL) { (*dzP) = NAN; } }
    else 
      { double h = zhi - zlo;
        double t = 2*(z - zlo)/h - 1;
        assert(fabs(t) <= 1.0);
        if (fabs(t) < 1.0e-10)
          { (*zP) = t; if (dzP != NULL) { (*dzP) = 2.0/h; } }
        else
          { double w = 1 - t*t;
            assert(w >= 0);
            (*zP) = t/w;
            if (dzP != NULL) { (*dzP) = 2*(1 + t*t)/(w*w*h); }
          }
      }
  }
                   
void contract_range(double *zP, double zlo, double zhi, double *dzP)
  { 
    double z = (*zP);
    if (isnan(z) || (zlo >= zhi))
      { (*zP) = NAN; if (dzP != NULL) { (*dzP) = NAN; } }
    else if (z == -INF)
      { (*zP) = zlo; if (dzP != NULL) { (*dzP) = NAN; } }
    else if (z == +INF)
      { (*zP) = zhi; if (dzP != NULL) { (*dzP) = NAN; } }
    else
      { 
        double h = zhi - zlo;
        double t, dt = 0;
        if (fabs(z) < 1.0e-10)
          { t = z; dt = 1; }
        else if (fabs(z) < 1.0)
          { double d = copysign(sqrt(1/(z*z) + 4), z);
            t = 0.5*(d - 1/z);
            if (dzP != NULL) { dt = 0.5*(1 - 1/(d*z))/(z*z); }
          }
        else
          { double d = sqrt(1 + 4*z*z);
            t = 0.5*(d - 1)/z;
            if (dzP != NULL) { dt = 0.5*(4/d - (d-1)/(z*z)); }
          }
        (*zP) = zlo + h*(t + 1)/2;
        if (dzP != NULL) { (*dzP) = h*dt/2; }
      }
  }
