#define PROG_NAME "test_jsmath"
#define PROG_DESC "test of {jsmath.h}, {ball.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-24 03:57:14 by stolfi */ 
/* Created on 2007-01-02 by J. Stolfi, UNICAMP */

#define test_jsmath_COPYRIGHT \
  "Copyright � 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <math.h>
#include <values.h>
#include <assert.h>
#include <float.h>

#include <jsmath.h>
#include <ball.h>
#include <jsrandom.h>
#include <jswsize.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <bool.h>

#define bug(FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

double zrandom(void); /* A random double with random magnitude. */

void test_iroundup(uint32_t nt);

void test_gcd(uint32_t nt);

void test_lcm(uint32_t nt);

void test_imod(uint32_t nt);

/* !!! should test ifloor, iceil, iroundfrac */

void test_ipow(uint32_t nt);

void test_upow(uint32_t nt);

void test_comb(uint32_t nt);

void test_64bit_mul(uint32_t nt);
  /* Tests {uint64_mul} and {int64_mul}. */
  
void test_minbits(uint32_t nt);
  
void test_digits(uint32_t nt);
  
void test_frac_digits(uint32_t nt);

void test_expand_contract(uint32_t nt);

void test_erf_inv(uint32_t nt);

void test_ball_vol(uint32_t nt);

typedef void range_map_t(double *zP, double zlo, double zhi, double *dzP);

double compute_map_derivative(double z, double zlo, double zhi, range_map_t *map);
  /* Computes the numeric derivative of {map} at {z}. The argumetns {zlo,zhi} 
    are passed on to {map}. */

void check_expand_contract_derivative(double z, double dz, double zlo, double zhi, range_map_t *map, char *mapname);
  /* Checks the analytic derivative {dz} against the numeric derivative of {map} at {z},
    compuetd by {compute_map_derivative(z,zlo,zhi,map)}. */

// #define N_double_nice 10
// 
// static double double_nice[N_double_nice] = 
//   { 
//     4.45014771701440277e-308, /* 8/DBL_MAX. */                   
//     7.45834073120020674e-155, /* 1/sqrt(DBL_MAX). */
//     5.00000000000000000e-001,
//     9.99999999999999889e-001,
//     1.00000000000000000e+000, 
//     1.00000000000000022e+000,
//     2.00000000000000000e+000,
//     1.34078079299425971e+154, /* sqrt(DBL_MAX). */
//     2.24711641857789464e+307, /* DBL_MAX/8. */                   
//     0
//   };

#define N_uint32_nice 13

static uint32_t uint32_nice[N_uint32_nice] = 
  { 1, 
    2, 
    3, 
    4, 
    12, 
    1023, 
    1024, 
    1025,
    2147483647UL,  /* 2^31-1 */
    2147483648UL, /* 2^31 */
    2147483649UL, /* 2^31+1 */
    4294967295UL, /* 2^32-1 */
    0
  };

#define N_uint64_nice 14

static uint64_t uint64_nice[N_uint64_nice] = 
  { 1, 
    2, 
    3, 
    4, 
    12, 
    1023, 
    1024, 
    1025, 
    4294967295ULL,      /* 2^32-1 */
    4294967296ULL,      /* 2^32 */
    4294967297ULL,      /* 2^32+1 */
    UINT64_MAX - 1ULL,  /* 2^64-2 */
    UINT64_MAX,         /* 2^64-1 */
    0
  };

#define N_int64_nice 28

static int64_t int64_nice[N_int64_nice] = 
  { +1,                       -1, 
    +2,                       -2, 
    +3,                       -3, 
    +4,                       -4, 
    +12,                      -12,
    +1023,                    -1023,                  
    +1024,                    -1024,                  
    +1025,                    -1025,                  
    +4294967295LL,            -4294967295LL,             /* �2^32-1 */
    +4294967296LL,            -4294967296LL,             /* �2^32 */
    +4294967297LL,            -4294967297LL,             /* �2^32+1 */
    INT64_MAX - 1LL,          1LL - INT64_MAX,           /* �2^63-2 */
    INT64_MAX,                -INT64_MAX,                /* �2^63-1 */
    INT64_MIN,                                           /* -2^63 */
    0
  };

int32_t main (int32_t argn, char **argv)
  { test_digits(1000);
    test_frac_digits(1000);
    test_minbits(1000);
    test_iroundup(100);
    test_gcd(200);
    test_lcm(200);
    test_imod(200);
    test_ipow(200);
    test_upow(200);
    test_comb(200);
    test_64bit_mul(200);
    test_expand_contract(200);
    test_erf_inv(2000);
    test_ball_vol(20);

    return 0;
  }

void test_iroundup(uint32_t nt)
  { fprintf(stderr, "Checking {iroundup,addrsync}...\n");
    int32_t i, j;
    for (i = 0; i < N_uint64_nice + nt; i++)
      { uint64_t a = (i < N_uint64_nice ? uint64_nice[i] : uint64_random());
        for (j = 0; j < N_uint64_nice + nt; j++)
          { uint64_t d = (j < N_uint64_nice ? uint64_nice[j] : uint64_abrandom(1,1000));
            if (d > 0)
              { uint64_t r = a % d;
                if ((r == 0) || (a <= UINT64_MAX - (d - r)))
                  { uint64_t b = iroundup(a, d);
                    demand(b >= a, "not monotonic");
                    demand(b % d == 0, "not multiple");
                    char *c = addrsync((char *)a, d);
                    demand((uint64_t)c == b, "addrsync inconsistent");
                  }
              }
          }
      }
  }


void test_gcd(uint32_t nt)
  { fprintf(stderr, "Checking {gcd}...\n");
    int32_t i,j;
    for (i = 0; i < N_uint64_nice + nt; i++)
      { uint64_t a = (i < N_uint64_nice ? uint64_nice[i] : uint64_random());
        for (j = 0; j < N_uint64_nice + nt; j++)
          { uint64_t b = (j < N_uint64_nice ? uint64_nice[j] : uint64_random());
            uint64_t c = gcd(a, b);
            if (c == 0)
              { affirm((a == 0) && (b == 0), "zero gcd"); }
            else
              { affirm(a % c == 0, "gcd does not divide a");
                affirm(b % c == 0, "gcd does not divide b");
                affirm(gcd(a/c, b/c) == 1, "gcd(a/c,b/c) is not 1");
                /* Should check maximality... */
              }
          }
      }
  }

void test_lcm(uint32_t nt)
  { fprintf(stderr, "Checking {lcm}...\n");
    int32_t i,j;
    for (i = 0; i < N_uint64_nice + nt; i++)
      { uint64_t a = (i < N_uint64_nice ? uint64_nice[i] : uint64_random());
        for (j = 0; j < N_uint64_nice + nt; j++)
          { uint64_t b = (j < N_uint64_nice ? uint64_nice[j] : uint64_random());
            uint64_t g = gcd(a, b);
            if ((a == 0) || (b == 0))
              { affirm(lcm(a,b) == 0, "lcm of zero is not zero"); }
            else if (a/g <= UINT64_MAX/b)
              { uint64_t c = lcm(a, b);
                affirm(c != 0, "zero lcm");
                affirm(c % a == 0, "lcm does not divide a");
                affirm(c % b == 0, "lcm does not divide b");
                /* Should check minimality... */
              }
          }
      }
  }

void test_imod(uint32_t nt)
  { fprintf(stderr, "Checking {imod}...\n");
    int32_t i,j;
    for (i = 0; i < N_int64_nice + nt; i++)
      { int64_t a = (i < N_int64_nice ? int64_nice[i] : int64_random());
        for (j = 0; j < N_int64_nice + nt; j++)
          { int64_t b = (j < N_int64_nice ? int64_nice[j] : int64_random());
            if (b > 0)
              { int64_t c = imod(a, b);
                affirm(c >= 0, "imod is negative");
                affirm(c < b, "imod is too big");
                /* Check {(a - c)%b}, but beware of overflow: */
                if (a < 0)
                  { affirm((a + (b - c)) % b == 0, "imod is not remainder"); }
                else
                  { affirm((a - c) % b == 0, "imod is not remainder"); }
                  
                /* Should check minimality... */
              }
          }
      }
  }

void test_ipow(uint32_t nt)
  { fprintf(stderr, "Checking {ipow}...\n");
    int32_t iy,ix;
    for (iy = 0; iy < N_uint32_nice + nt; iy++)
      { uint32_t y = (iy < N_uint32_nice ? uint32_nice[iy] : (uint32_t)int64_random());
        for (ix = 0; ix < N_int64_nice + nt; ix++)
          { int64_t x = (ix < N_int64_nice ? int64_nice[ix] : int64_random());
            if ((y > 1) && ((x > +1) || (x < -1)))
              { /* Reduce {y} if needed so that the power will not overflow: */
                uint32_t ymax;
                if ((x > 0) || (y % 2 == 0))
                  { int64_t pmax = INT64_MAX; /* 2^63-1 */
                    ymax = 0; do { pmax /= x; ymax++; } while ((pmax/x) != 0);
                  }
                else
                  { int64_t pmax = INT64_MIN; /* -2^63 */
                    ymax = 0; do { pmax /= x; ymax++; } while ((pmax/x) != 0);
                  }
                if (y > ymax) { y = ymax; }  
              }
            /* Compute the power {z}: */
            int64_t z = ipow(x, y);
            /* Check whether {x} divides {z} {y} times exactly, leaving 1: */
            if (y == 0)
              { affirm(z == 1, "x^0 is not 1"); }
            else if (x == 0)
              { affirm(z == 0, "0^y is not 0"); }
            else if (x == +1)
              { affirm(z == +1, "1^y is not 1"); }
            else if (x == -1)
              { affirm(z == (y % 2 == 0 ? +1 : -1), "(-1)^y is not (-1)^(y%2)"); }
            else if (y == 1)
              { affirm(z == x, "x^1 is not x"); }
            else
              { int64_t p = z;
                while ((p > 1) || (p < -1))
                  { int64_t q = p/x;
                    affirm(p - q*x == 0, "x^y is not a power of x"); 
                    p = q;
                  }
                if (p != 1) { fprintf(stderr, "x = %ld y = %u z = %ld p = %ld\n", x, y, z, p); }
                affirm(p == 1, "power is not divisible y times by x");
              }
          }
      }
  }
  
void test_upow(uint32_t nt)
  { fprintf(stderr, "Checking {upow}...\n");
    int32_t iy,ix;
    for (iy = 0; iy < N_uint32_nice + nt; iy++)
      { uint32_t y = (iy < N_uint32_nice ? uint32_nice[iy] : (uint32_t)int64_random());
        for (ix = 0; ix < N_int64_nice + nt; ix++)
          { uint64_t x = (ix < N_uint64_nice ? uint64_nice[ix] : uint64_random());
            if ((y > 1) && (x > 1))
              { /* Reduce {y} if needed so that the power will not overflow: */
                uint32_t ymax;
                uint64_t pmax = UINT64_MAX; /* 2^64-1 */
                ymax = 0; do { pmax /= x; ymax++; } while ((pmax/x) != 0);
                if (y > ymax) { y = ymax; }  
              }
            /* Compute the power {z}: */
            uint64_t z = upow(x, y);
            /* Check whether {x} divides {z} {y} times exactly, leaving 1: */
            if (y == 0)
              { affirm(z == 1, "x^0 is not 1"); }
            else if (x == 0)
              { affirm(z == 0, "0^y is not 0"); }
            else if (x == 1)
              { affirm(z == 1, "1^y is not 1"); }
            else if (y == 1)
              { affirm(z == x, "x^1 is not x"); }
            else
              { uint64_t p = z;
                while (p > 1)
                  { uint64_t q = p/x;
                    affirm(p == q*x, "x^y is not a power of x"); 
                    p = q;
                  }
                if (p != 1) { fprintf(stderr, "x = %lu y = %u z = %lu p = %lu\n", x, y, z, p); }
                affirm(p == 1, "power is not divisible y times by x");
              }
          }
      }
  }
  
void test_comb(uint32_t nt)
  { fprintf(stderr, "Checking {comb}...\n");
    int32_t t;
    for (t = 0; t < nt; t++)
      { int32_t maxpos = 30;
        int32_t maxneg = 2;
        int32_t n = (17*t + t*t) % (maxpos + maxneg + 1) - maxneg;
        int32_t k = (43*t + t*t/8) % (n + 2*maxneg) - maxneg;
        /* Compute {comb(n,k)} as {double} to get an idea of the magnitude: */
        double fr;
        bool_t safe; /* Believed to be safe to compute. */
        if ((k < 0) || (k > n))
          { fr = 0; safe = TRUE; }
        else if ((k == 0) || (k == n))
          { fr = 1; safe = TRUE; }
        else
          { int32_t m = (k <= n-k ? k : n-k);
            fr = 1;
            int32_t i;
            for (i = 1; i <= m; i++) { fr = (fr*((double)n-m+i))/((double)i); }
            /* Guard against overflow: */
            double maxfr = exp2(63)/((double)m); /* Max internal temp value sems to be {comb(n,m)*m}. */
            safe = (fr < 0.9*maxfr);
          }
        if (safe)
          { /* Seems small enough o avoid overflow: */
            uint64_t r = comb((uint64_t)n, (uint64_t)k);
            if((double)r != fr)
              { fprintf(stderr, ("n = %d  k = %d  r = %" uint64_u_fmt "  fr = %28.3f\n"), n, k, r, fr);
                affirm(FALSE, "{comb} failed (3)");
              }
          }
      }
  }
 
void test_frac_digits(uint32_t nt)
  { fprintf(stderr, "Checking {frac_digits}...\n");
    double x = 0; /* Test arg. */
    for (uint32_t k = 0;  k < nt; k++)
      { uint32_t m_cmp = frac_digits(x);
        if (x == 0)
          { demand(m_cmp == UINT32_MAX, "{frac_digits} failed for {x=0}"); }
        else if (fabs(x) <= DBL_MIN)
          { demand(m_cmp == 308, "{frac_digits} failed for {x<=DBL_MIN}"); }
        else
          { demand(m_cmp <= 308, "invalid {frac_digits} result"); 
            double p = pow(0.1, (double)m_cmp);
            if (p <= 0) { p = 0.1*DBL_MIN; /* Denenormalized but positive. */ }
            demand(fabs(x) >= p, "{frac_digits} failed");
          }
        /* Advance to next {x}: */
        x = 5*x;
        if (! isfinite(x)) { x = sqrt(k)*DBL_MIN; }
      }
   }        
 
void test_digits(uint32_t nt)
  { bool_t debug = FALSE;
    fprintf(stderr, "Checking {digits}...\n");
    uint64_t x = 0; /* Test arg. */
    for (uint32_t k = 0;  k < nt; k++)
      { uint32_t m_cmp = digits(x);
        if (debug) { fprintf(stderr, "  digits(%ld) = %d", x, m_cmp); } 
        uint32_t m_exp; /* Expected result. */
        if (x == 0)
          { m_exp = 1; }
        else
          { char *sx = jsprintf("%ld", x);  
            m_exp = (uint32_t)strlen(sx);
            free(sx);
          }
        if (debug) { fprintf(stderr, " expected %d\n", m_exp); } 
        demand(m_cmp == m_exp, "{digits} failed");
        /* Advance to next {x}: */
        if (k < nt/4)
          { x = x + 1; }
        else
          { x = (7*x/5) | 1; }
      }
   }        
 
void test_minbits(uint32_t nt)
  { fprintf(stderr, "Checking {minbits}...\n");
    uint64_t x = 0; /* Test arg. */
    uint32_t d = 0; /* Expected {minbits(x)}. */
    uint64_t p = 0; /* {2^d-1}. */
    for (uint32_t k = 0;  k < nt; k++)
      { /* Check {minbits(x) == d}: */
        uint32_t m = minbits(x);
        if (m != d)  { fprintf(stderr, "x = %lu  m = %u  d = %u\n", x, m, d); }
        demand(m == d, "minbits failed");
        /* Advance to next {x}: */
        if (k < nt/4)
          { x = x + 1; }
        else
          { x = (7*x/5) | 1; }
        /* Adjust {d,p}: */
        while (x > p) { d++; p = (p << 1) | 1; }
        while (x <= (p >> 1)) { d--; p = (p >> 1); }
        assert((x >= (p/2)+1) && (x <= p));
      }
   }        
  
void test_64bit_mul(uint32_t nt)
  { fprintf(stderr, "Checking {uint64_mul,int64_mul}...\n");
    bool_t verbose = FALSE;
    int32_t i,j;
    for (i = 0; i < N_uint64_nice + nt; i++)
      { uint64_t x = (i < N_uint64_nice ? uint64_nice[i] : uint64_random());
        int64_t sx = (int64_t)x;
        for (j = 0; j < N_uint64_nice + nt; j++)
          { uint64_t y = (j < N_uint64_nice ? uint64_nice[j] : uint64_random());
            int64_t sy = (int64_t)y;
            
            /* Test {uint64_mul}: */
            uint64_t Z1, Z0;
            uint64_mul(x, y, &Z1, &Z0);
            /* Slow multiplication: */
            int32_t i;
            uint64_t P1 = 0, P0 = 0;
            uint64_t M = 1;
            for (i = 0; i < 64; i++)
              { if ((M & x) != 0)
                  { /* Add {y} shifted by {i} onto {P}: */
                    uint64_t Q = P0 + (y << i);
                    if (Q < P0) { /* carry happened: */ P1 += 1; }
                    P0 = Q;
                    if (i > 0) { P1 += (y >> (64 - i)); }
                  }
                M = (M << 1);
              }
                    
            if ((P0 != Z0) || (P1 != Z1) || verbose)
              { fprintf(stderr, ("x = %" uint64_u_fmt "  y = %" uint64_u_fmt), x, y);
                fprintf(stderr, ("  Z = %" uint64_u_fmt " *2^64 + %" uint64_u_fmt), Z1, Z0);
                fprintf(stderr, ("  P = %" uint64_u_fmt " *2^64 + %" uint64_u_fmt "\n"), P1, P0);
                if ((P0 != Z0) || (P1 != Z1)) { affirm(FALSE, "bug"); }
              }
            
            /* Test {int64_mul}: */
            int64_t W1; uint64_t W0;
            int64_mul(sx, sy, &W1, &W0);
            
            int64_t Q1 = (int64_t)P1; uint64_t Q0 = P0;
            if (sx < 0) { Q1 = Q1 - sy; }
            if (sy < 0) { Q1 = Q1 - sx; }
            if ((Q0 != W0) || (Q1 != W1) || verbose)
              { fprintf(stderr, ("sx = %" int64_d_fmt "  sy = %" int64_d_fmt), x, y);
                fprintf(stderr, ("  W = %" int64_d_fmt " *2^64 + %" uint64_u_fmt), W1, W0);
                fprintf(stderr, ("  Q = %" int64_d_fmt " *2^64 + %" uint64_u_fmt "\n"), Q1, Q0);
                if ((Q0 != W0) || (Q1 != W1)) { affirm(FALSE, "bug"); }
              }
          }
      }
  }

void test_expand_contract(uint32_t nt)
  { 
    fprintf(stderr, "Checking {expand_range,contract_range}...\n");
    int32_t i;
    for (i = 0; i < nt; i++)
      { 
        double zmid = 10*drandom() - 5;
        double zrad = 3*drandom() + 0.01;
        double zlo = zmid - zrad;
        double zhi = zmid + zrad;
        
        if (i == 0)
          { 
            /* Plot maps: */
            FILE *wr = open_write("/tmp/foo.dat", TRUE);
            int32_t ns = 101;
            int32_t k;
            for (k = 0; k < ns; k++)
              { double zk = zlo + (zhi - zlo)*(k + 0.5)/((double)ns);
                double zke = zk;
                double dzke_num = compute_map_derivative(zke, zlo, zhi, &expand_range);
                double dzke_alg;
                expand_range(&zke, zlo,zhi, &dzke_alg);
                double zkc = zke;
                double dzkc_num = compute_map_derivative(zkc, zlo, zhi, &contract_range);
                double dzkc_alg;
                contract_range(&zkc, zlo,zhi, &dzkc_alg);
                fprintf
                  ( wr, "%16.8f  %16.8f %16.8f %16.8f   %16.8f %16.8f %16.8f\n",
                    zk,  zke, dzke_alg, dzke_num,  zkc, dzkc_alg, dzkc_num
                  );
              }
            fclose(wr);
          }
        
        double za = (drandom() < 0.800 ? zlo + (zhi-zlo)*drandom() : (drandom() < 0.500 ? zlo : zhi));
        
        /* Check the expansion map: */
        double zb = za, dzb;
        expand_range(&zb, zlo,zhi, &dzb);
        if (((za == zlo) && (zb != -INF)) || ((za == zhi) && (zb != +INF)))
          { fprintf(stderr, "za = %24.16e  zlo = %24.16e  zhi = %24.16e  zb = %24.16e\n", za, zlo, zhi, zb);
            affirm(FALSE, "bug");
          }
        if ((za > zlo) && (za < zhi))
          { /* Check the derivative: */
            check_expand_contract_derivative(za, dzb, zlo, zhi, &expand_range, "expand_range");
          }
        
        /* Check the contraction map: */
        double zc = zb, dzc;
        contract_range(&zc, zlo,zhi, &dzc);
        if (((zb == -INF) && (zc != zlo)) || ((zb == +INF) && (zc != zhi)))
          { fprintf(stderr, "za = %24.16e  zb = %24.16e  zlo = %24.16e  zhi = %24.16e  zc = %24.16e\n", za, zb, zlo, zhi, zc);
            affirm(FALSE, "bug");
          }
        if ((za > zlo) && (za < zhi))
          { /* Check the derivative: */
            check_expand_contract_derivative(zb, dzc, zlo, zhi, &contract_range, "contract_range");
          }
          
        if (fabs(zc - za) > 1.0e-8)
          { fprintf(stderr, "za = %24.16e  zb = %24.16e  zlo = %24.16e  zhi = %24.16e  zc = %24.16e\n", za, zb, zlo, zhi, zc);
            affirm(FALSE, "bug");
          }
      }
  }

void test_erf_inv(uint32_t nt)
  {
   fprintf(stderr, "Checking {erf_inv}...\n");
   double emax = 0.0; /* Max error {p - erf(erf_inv(p))}. */
   double p1emax = 0.0; /* {p1} for which the erroris maximum. */
   double z1emax = 0.0; /* {z1} for which the erroris maximum. */
   double p2emax = 0.0; /* {p2} for which the erroris maximum. */
   for (uint32_t k = 0;  k <= nt; k++)
     {
       double z0 = 2*((double)k)/((double)nt) - 1.0;
       double p1 = erf(10*z0);
       double z1 = erf_inv(p1);
       double p2 = erf(z1);
       double err = p1 - p2;
       if (fabs(emax) < fabs(err)) 
         { emax = err; p1emax = p1; z1emax = z1; p2emax = p2; }
       if (err > 1.0e-12)
         { fprintf(stderr, "  p1 = %25.20f  z1 = %25.20f p2 = %25.20f", p1, z1, p2);
           fprintf(stderr, "  err = %25.16e\n", err); 
         }
     }
    fprintf(stderr, "  max error:\n");
    fprintf(stderr, "  p1 = %25.20f  z1 = %25.20f p2 = %25.20f", p1emax, z1emax, p2emax);
    fprintf(stderr, "  err = %25.16e\n", emax); 
  }

void test_ball_vol(uint32_t nt)
  {
    for (uint32_t d = 0;  d <= nt; d++)
      { bool_t verbose = (d <= 5);
        /* TEST: double ball_vol(int32_t d); */
        /* TEST: double ball_cap_vol_frac_pos(int32_t d, double u); */
        /* TEST: double ball_zone_vol_frac_ang(int32_t d, double w); */

        if (verbose) { fprintf(stderr, "Checking {ball_vol} d = %d...\n", d); }
        { 
          double vball = ball_vol(d);
          if (verbose) { fprintf(stderr, "  measure of %d-ball with unit radius = %12.7f\n", d, vball); } 
          /* Checking: */
          if (d <= 5)
            { double vv;
              if (d == 0)
                { vv = 1; }
              else if (d == 1)
                { vv = 2; }
              else if (d == 2)
                { vv = M_PI; }
              else if (d == 3)
                { vv = M_PI*4/3; }
              else if (d == 4)
                { vv = M_PI*M_PI/2; }
              else if (d == 5)
                { vv = M_PI*M_PI*8/15; }
              else 
                { assert(FALSE); vv = 0.0; }
              affirm(fabs(vv - vball) <= 1.0e-6, "{ball_vol} error");
            }
        }

        if (verbose) { fprintf(stderr, "Checking {ball_zone_vol_frac_ang} d = %d...\n", d); }
        { int32_t NW = 10; /* Number of latitude steps. */
          if (verbose) { fprintf(stderr, "  volume fraction between lat = 0 and lat = z:\n"); }
          for (int32_t i = -NW; i <= NW; i++)
            { double w = M_PI/2*((double)i)/((double)NW);
              double f = ball_zone_vol_frac_ang(d, w);
              if (verbose) { fprintf(stderr, "    %+8.5f  %+8.5f\n", w, f); } 
              /* Checking: */
              if (d <= 4)
                { double ff;
                  if (i <= -NW)
                    { ff = -0.5; }
                  else if (i >= +NW)
                    { ff = +0.5; }
                  else if (d == 0)
                    { ff = 0.0; }
                  else if (d == 1)
                    { ff = sin(w)/2; }
                  else  if (d == 2)
                    { ff = (w + sin(w)*cos(w))/M_PI; }
                  else if (d == 3)
                    { ff = sin(w)*(3 - sin(w)*sin(w))/4; }
                  else if (d == 4)
                    { ff = (w + 2*sin(2*w)/3 + sin(4*w)/12)/M_PI; }
                  else 
                    { assert(FALSE); ff = 0.0; }
                  affirm(fabs(ff - f) <= 1.0e-6, "{ball_zone_vol_frac_ang} error");
                }
              if (i == -NW)
                { affirm (fabs(f + 0.5) < 1.0e-8, "{ball_zone_vol_frac_ang} error for {w = -PI/2}"); }
              else if (i == +NW)
                { affirm (fabs(f - 0.5) < 1.0e-8, "{ball_zone_vol_frac_ang} error for {w = +PI/2}"); }
              else if (i == 0)
                { affirm (f == 0, "{ball_zone_vol_frac_ang} error for {w = 0}"); }
            }
          if (verbose) { fprintf(stderr, "\n"); }
        }

        if (verbose) { fprintf(stderr, "Checking {ball_cap_vol_frac_pos} d = %d\n", d); }
        { int32_t NX = 10; /* Number of position steps in each hemisphere. */
          int32_t imin = -NX-1;
          int32_t imax = +NX+1;
          if (verbose) { fprintf(stderr, "  volume fraction between x = -1 and x = u:\n"); }
          for (int32_t i = imin; i <= imax; i++)
            { double u = ((double)i)/((double)NX);
              double f = ball_cap_vol_frac_pos(d, u);
              if (verbose) { fprintf(stderr, "  %8.5f  %8.5f\n", u, f); } 
              double ff;
              if (u < -1)
                { ff = 0.0; }
              else if (u == -1)
                { ff = (d == 0 ? 0.25 : 0.0); }
              else if (u > +1)
                { ff = 1.0; }
              else if (u == +1)
                { ff = (d == 0 ? 0.75 : 1.0); }
              else
                { /* Checking consistency with {ball_zone_vol_frac_ang}: */
                  double w = asin(u);
                  double vw = ball_zone_vol_frac_ang(d, w);
                  ff = 0.5 + vw;
                }
              affirm(fabs(ff - f) <= 1.0e-6, "{ball_cap_vol_frac_pos} error");
            }
          if (verbose) { fprintf(stderr, "\n"); }
        }
    }
  }

double compute_map_derivative(double z, double zlo, double zhi, range_map_t *map)
  { 
    double step = 1.0e-6; /* Numeric derivative step. */
    double zp = z+step, zm = z-step;
    map(&zp, zlo,zhi, NULL);
    map(&zm, zlo,zhi, NULL);
    if (isfinite(zm) && isfinite(zp))
      { return (zp - zm)/(2*step); }
    else
      { return NAN; }
  }

void check_expand_contract_derivative(double z, double dz, double zlo, double zhi, range_map_t *map, char *mapname)
  { 
    double rtol = 1.0e-4;   /* Error tolerance relative to derivative values. */
    double atol = 1.0e-10;  /* Error tolerance relative to {z} values. */
    double dnum = compute_map_derivative(z, zlo, zhi, map);
    if (isfinite(dnum))
      { double err = fabs(dnum-dz);
        double err_rmax = rtol*(fabs(dnum)+fabs(dz));
        double err_amax = atol*fmax(fabs(zlo),fabs(zhi));
        double err_max = fmax(err_rmax, err_amax);
        if ((! isfinite(dz)) || (err > err_max))
          { fprintf(stderr, " ** Derivative of %s is inconsistent (err = %+12.6e)\n", mapname, err); 
            fprintf(stderr, "  z = %24.16e  zlo = %24.16e  zhi = %24.16e\n", z, zlo, zhi);
            fprintf(stderr, "  dz = %24.16e numeric = %24.16e", dz, dnum);
            fprintf(stderr, "\n");
            affirm(FALSE, "Aborted");
          }
      }
  }

#define LOG_DBL_MAX (709.78271289338397)

double zrandom(void)
  { 
    return drandom()*exp(LOG_DBL_MAX*(2*drandom() - 1));
  }
