#define PROG_NAME "test_jsmath"
#define PROG_DESC "test of {jsmath.h}"
#define PROG_VERS "1.0"

/* Last edited on 2012-12-20 20:16:19 by stolfilocal */ 
/* Created on 2007-01-02 by J. Stolfi, UNICAMP */

#define test_jsmath_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <values.h>

#include <jsmath.h>
#include <jsrandom.h>
#include <jswsize.h>
#include <jsfile.h>
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

void test_gcd(int nt);

void test_lcm(int nt);

void test_imod(int nt);
/* !!! should test ifloor, iceil */

void test_ipow(int nt);

void test_comb(int nt);

void test_64bit_mul(int nt);
  /* Tests {uint64_mul} and {int64_mul}. */
  
void test_expand_contract(int nt);

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
    4294967295ULL, /* 2^32-1 */
    4294967296ULL, /* 2^32 */
    4294967297ULL, /* 2^32+1 */
    18446744073709551614ULL, /* 2^64-2 */
    18446744073709551615ULL, /* 2^64-1 */
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
    +4294967295LL,            -4294967295LL,             /* ±2^32-1 */
    +4294967296LL,            -4294967296LL,             /* ±2^32 */
    +4294967297LL,            -4294967297LL,             /* ±2^32+1 */
    +9223372036854775806LL,   -9223372036854775806LL,    /* ±2^63-2 */
    +9223372036854775807LL,   -9223372036854775807LL,    /* ±2^63-1 */
    ((int64_t)-9223372036854775807LL)-1,                 /* -2^63 */
    0
  };

int main (int argn, char **argv)
  { test_gcd(200);
    test_lcm(200);
    test_imod(200);
    test_ipow(200);
    test_comb(200);
    test_64bit_mul(200);
    test_expand_contract(200);

    return 0;
  }

void test_gcd(int nt)
  { fprintf(stderr, "Checking {gcd}...\n");
    int i,j;
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

void test_lcm(int nt)
  { fprintf(stderr, "Checking {lcm}...\n");
    int i,j;
    for (i = 0; i < N_uint64_nice + nt; i++)
      { uint64_t a = (i < N_uint64_nice ? uint64_nice[i] : uint64_random());
        for (j = 0; j < N_uint64_nice + nt; j++)
          { uint64_t b = (j < N_uint64_nice ? uint64_nice[j] : uint64_random());
            uint64_t g = gcd(a, b);
            if ((a == 0) || (b == 0))
              { affirm(lcm(a,b) == 0, "lcm of zero is not zero"); }
            else if (a/g <= UINT64_MAX/b)
              { int64_t c = lcm(a, b);
                affirm(c != 0, "zero lcm");
                affirm(c % a == 0, "lcm does not divide a");
                affirm(c % b == 0, "lcm does not divide b");
                /* Should check minimality... */
              }
          }
      }
  }

void test_imod(int nt)
  { fprintf(stderr, "Checking {imod}...\n");
    int i,j;
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

void test_ipow(int nt)
  { fprintf(stderr, "Checking {ipow}...\n");
    int iy,ix;
    for (iy = 0; iy < N_uint32_nice + nt; iy++)
      { uint32_t y = (iy < N_uint32_nice ? uint32_nice[iy] : (uint32_t)int64_random());
        for (ix = 0; ix < N_int64_nice + nt; ix++)
          { int64_t x = (ix < N_int64_nice ? int64_nice[ix] : int64_random());
            if ((y > 1) && ((x > +1) || (x < -1)))
              { /* Reduce {y} if needed so that the power will not overflow: */
                uint32_t ymax;
                if ((x > 0) || (y % 2 == 0))
                  { int64_t pmax = +9223372036854775807LL; /* 2^63-1 */
                    ymax = 0; do { pmax /= x; ymax++; } while ((pmax/x) != 0);

                  }
                else
                  { int64_t pmax = ((int64_t)-9223372036854775807LL)-1; /* -2^63 */
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
                if (p != 1) { fprintf(stderr, "x = %lld y = %u z = %lld p = %lld\n", x, y, z, p); }
                affirm(p == 1, "power is not divisible y times by x");
              }
          }
      }
  }
  
void test_comb(int nt)
  { fprintf(stderr, "Checking {comb}...\n");
    int t;
    for (t = 0; t < nt; t++)
      { int maxpos = 30;
        int maxneg = 2;
        int n = (17*t + t*t) % (maxpos + maxneg + 1) - maxneg;
        int k = (43*t + t*t/8) % (n + 2*maxneg) - maxneg;
        /* Compute {comb(n,k)} as {double} to get an idea of the magnitude: */
        double fr;
        bool_t safe; /* Believed to be safe to compute. */
        if ((k < 0) || (k > n))
          { fr = 0; safe = TRUE; }
        else if ((k == 0) || (k == n))
          { fr = 1; safe = TRUE; }
        else
          { int m = (k <= n-k ? k : n-k);
            fr = 1;
            int i;
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
  
void test_64bit_mul(int nt)
  { fprintf(stderr, "Checking {uint64_mul,int64_mul}...\n");
    bool_t verbose = FALSE;
    int i,j;
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
            int i;
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

void test_expand_contract(int nt)
  { 
    fprintf(stderr, "Checking {expand_range,contract_range}...\n");
    int i;
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
            int ns = 101;
            int k;
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
    double tol = 1.0e-4;  /* Relative error tolerance. */
    double dnum = compute_map_derivative(z, zlo, zhi, map);
    if (isfinite(dnum))
      { double err = fabs(dnum-dz);
        if (!isfinite(dz) || (err > tol*(fabs(dnum)+fabs(dz))))
          { fprintf(stderr, " ** Derivative of %s is inconsistent (err = %+12.6e)\n", mapname, err); 
            fprintf(stderr, "  z = %24.16e  zlo = %24.16e  zhi = %24.16e\n", z, zlo, zhi);
            fprintf(stderr, "  dz = %14.8f numeric = %14.8f", dz, dnum);
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
