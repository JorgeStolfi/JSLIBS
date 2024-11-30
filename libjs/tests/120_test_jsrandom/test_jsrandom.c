#define PROG_NAME "test_jsrandom"
#define PROG_DESC "test of {jsmath.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-26 22:54:37 by stolfi */ 
/* Created on 2011-09-20 by J. Stolfi, UNICAMP */

#define test_jsrandom_COPYRIGHT \
  "Copyright © 2011  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsmath.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <bool.h>

#include <jsrandom.h>

#define bug(FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

double zrandom(void); /* A random double with random magnitude. */

#define MAX_HIST_BINS 1024
  /* Max number of bins in histograms. */
  
void test_int32_random(uint32_t nb, uint32_t ntpb, char *xtest);
void test_uint32_random(uint32_t nb, uint32_t ntpb, char *xtest);
void test_int64_random(uint32_t nb, uint32_t ntpb, char *xtest);
void test_uint64_random(uint32_t nb, uint32_t ntpb, char *xtest);
  /* Tests {function = {int,uint}{32,64}_random}. Writes to disk
    histograms of the high and low {nb} bits, with {ntpd}
    samples per bin on average.  File names are 
    "out/{function}_{"lo","hi"}_{xtest}.his". */

void test_stdint_random(uint32_t size, bool_t sgn, uint32_t nb, uint32_t ntpb, char *xtest);
  /* Tests {function = {xtype}{size}_random} where {size} is 32 or 64
    and {xtype} is {int} if {sgn} is true, {uint} if {sgn} is false. 
    Writes to disk histograms of the high and low {nb} bits, with {ntpd}
    samples per bin on average.  File names are "out/{function}_{"lo","hi"}_{xtest}.his". */

void test_int32_abrandom(int32_t a, int32_t b, uint32_t q, uint32_t ntpb, char *xtest);
void test_uint32_abrandom(uint32_t a, uint32_t b, uint32_t q, uint32_t ntpb, char *xtest);
void test_int64_abrandom(int64_t a, int64_t b, uint64_t q, uint32_t ntpb, char *xtest);
void test_uint64_abrandom(uint64_t a, uint64_t b, uint64_t q, uint32_t ntpb, char *xtest);
  /* Tests {function(a,b) = {int,uint}{32,64}_abrandom(a,b)}. Writes to disk
    a histogram of values with {q} bins and
    {ntpd} samples per bin on average.  The parameter {q} must divide {b-a}.
    File names are "out/{function}_qt_{xtest}.his". */

void test_stdint_abrandom(uint32_t size, bool_t sgn, uint64_t a, uint64_t b, uint64_t q, uint32_t ntpb, char *xtest);
  /* Tests {function(a,b) = {xtype}{size}_abrandom(a,b)} where {size} is 32 or 64
    and {xtype} is {int} if {sgn} is true, {uint} if {sgn} is false.
    Assumes that {a,b} are the original arguments cast to the 
    type {uint64_t} without any shifting. 
    
    Writes to disk histograms of the values, with {q} bins and
    {ntpd} samples per bin on average.  The parameter {q} must divide {b-a}.
    File name is "out/{function}_qt_{xtest}.his". */

void test_drandom(uint32_t ntpb, char *xtest);
void test_dgaussrand(uint32_t ntpb, char *xtest);
void test_dloggaussrand(double avg, double dev, uint32_t ntpb, char *xtest);
  /* Test {function = drandom,dgaussrand,dloggaussrand}. Writes to disk a histogram
    of values with {ntpd} samples per bin on average.
    File names are "out/{function}_qt_{xtest}.his". */

void test_uint64_choose(uint64_t n, size_t k, uint64_t nt, char *xtest);
  /* Tests {uint64_choose} with parameters {n} and {k}, {nt} times. */

void test_random_perm(uint32_t n, uint64_t nt, char *xtest);
  /* Tests {random_perm} with parameter {n}, {nt} times. */

FILE *open_test_file(char *xtype, uint32_t size, char *func, char *xend, char *xtest);
  /* Opens "out/{xtype}{size}_{func}_hist_{xend}_{xtest}.his" for writing.
    If {xtype} is {NULL}, omits the "{xtype}{size}_" part. */

#define P03 (1000)
#define P06 (1000000)
#define P07 (10000000)
#define P08 (100000000)
#define P09 (1000000000)
#define P12 (1000000000000L)
#define P15 (1000000000000000L)
#define P18 (1000000000000000000L)
  /* Powers of 10. */

#define T08 (1<<8)
#define T10 (1<<10)
#define T40 (1LL<<40)
#define T60 (1LL<<60)
  /* Powers of 2. */

int32_t main (int32_t argn, char **argv)
  { 
    test_int32_random(6, 10000, "A");
    test_uint32_random(6, 10000, "A");
    test_int64_random(6, 10000, "A");
    test_uint64_random(6, 10000, "A");
    
    test_int32_abrandom(-200, +99, 1, 10000, "A");
    test_uint32_abrandom(20000, 50000-1, 100, 10000, "A");
    test_int64_abrandom(-200*P07, +100*P07-1, P09, 10000, "A");
    test_uint64_abrandom(200*P12, 500*P12-1, P12, 10000, "A");

    test_int32_abrandom(0-30, +T08-1-30, T08/256, 10000, "B");
    test_uint32_abrandom(0+2000, T10-1+2000, T10/256, 10000, "B");
    test_int64_abrandom(0-P12, +T40-1-P12, T40/256, 10000, "B");
    test_uint64_abrandom(0+P18, T60-1+P18, T60/256, 10000, "B");

    test_drandom(10000, "A");
    test_dgaussrand(10000, "A");
    test_dloggaussrand(7.0, 2.0, 10000, "A");
    
    test_uint64_choose(10,5,160000,"C");
    test_uint64_choose(100,1,160000,"C");
    test_uint64_choose(100,10,160000,"C");
    test_uint64_choose(100,100,160000,"C");
    test_uint64_choose(1000,10,160000,"C");
    
    test_random_perm(10,500,"D");
    test_random_perm(10000,500,"D");

    return 0;
  }

void test_int32_random(uint32_t nb, uint32_t ntpb, char *xtest)
  {
    test_stdint_random(32, TRUE, nb, ntpb, xtest);
  }

void test_uint32_random(uint32_t nb, uint32_t ntpb, char *xtest)
  {
    test_stdint_random(32, FALSE, nb, ntpb, xtest);
  }

void test_int64_random(uint32_t nb, uint32_t ntpb, char *xtest)
  {
    test_stdint_random(64, TRUE, nb, ntpb, xtest);
  }

void test_uint64_random(uint32_t nb, uint32_t ntpb, char *xtest)
  {
    test_stdint_random(64, FALSE, nb, ntpb, xtest);
  }

void test_stdint_random(uint32_t size, bool_t sgn, uint32_t nb, uint32_t ntpb, char *xtest)
  { char *xtype = (sgn ? "int" : "uint");
    fprintf(stderr, "Checking {%s%d_random}...\n", xtype, size);
    srandom(19501129);
    
    /* Allocate and initialize histogram: */
    uint32_t nh = (uint32_t)ipow(2, nb);  /* Number of histogram bins, a power of 2. */
    uint64_t hi[nh], lo[nh]; /* Histograms of high and low bits. */
    for (uint32_t j = 0;  j < nh; j++) { lo[j] = hi[j] = 0; }

    /* Collect histogram data: */
    uint64_t nt = ntpb*nh;      /* Total number of samples. */
    for (int64_t i = 0; i < nt; i++)
      { uint32_t jlo, jhi; /* Lower and higher {nb} bits. */
        if (size == 32)
          { uint32_t x = (sgn ? (uint32_t)int32_random() : uint32_random());
            jlo = (uint32_t)(x & (uint32_t)(nh - 1));
            jhi = (uint32_t)(x >> (32 - nb));
          }
        else if (size == 64)
          { uint64_t x = (sgn ? (uint64_t)int64_random() : uint64_random());
            jlo = (uint32_t)(x & (uint32_t)(nh - 1));
            jhi = (uint32_t)(x >> (64 - nb));
          }
        else
          { assert(FALSE); }
          
        lo[jlo]++;
        hi[jhi]++;
      }

    /* Write histograms: */
    FILE *wrhi = open_test_file(xtype, size, "random", "hi", xtest); 
    for (uint32_t j = 0;  j < nh; j++) 
      { double r = ((double)hi[j])/((double)ntpb);
        fprintf(wrhi, "%03d %12ld %12.8f\n", j, hi[j], r);
      }
    fclose(wrhi);

    FILE *wrlo = open_test_file(xtype, size, "random", "lo", xtest); 
    for (uint32_t j = 0;  j < nh; j++) 
      { double r = ((double)lo[j])/((double)ntpb);
        fprintf(wrlo, "%03d %12ld %12.8f\n", j, lo[j], r);
      }
    fclose(wrlo);
  }

void test_int32_abrandom(int32_t a, int32_t b, uint32_t q, uint32_t ntpb, char *xtest)
  {
    test_stdint_abrandom(32, TRUE, (uint64_t)a, (uint64_t)b, (uint64_t)q, ntpb, xtest);
  }

void test_uint32_abrandom(uint32_t a, uint32_t b, uint32_t q, uint32_t ntpb, char *xtest)
  {
    test_stdint_abrandom(32, FALSE, (uint64_t)a, (uint64_t)b, (uint64_t)q, ntpb, xtest);
  }

void test_int64_abrandom(int64_t a, int64_t b, uint64_t q, uint32_t ntpb, char *xtest)
  {
    test_stdint_abrandom(64, TRUE, (uint64_t)a, (uint64_t)b, (uint64_t)q, ntpb, xtest);
  }

void test_uint64_abrandom(uint64_t a, uint64_t b, uint64_t q, uint32_t ntpb, char *xtest)
  {
    test_stdint_abrandom(64, FALSE, (uint64_t)a, (uint64_t)b, (uint64_t)q, ntpb, xtest);
  }

void test_stdint_abrandom(uint32_t size, bool_t sgn, uint64_t a, uint64_t b, uint64_t q, uint32_t ntpb, char *xtest)
  { char *xtype = (sgn ? "int" : "uint");
    fprintf(stderr, "Checking {%s%d_abrandom}...\n", xtype, size);
    srandom(19501129);
    
    /* Print range {a..b} recast to original type, compute num values {nv1+1}: */
    uint64_t nv1; /* Num values minus one. */
    if (sgn)
      { int64_t sa, sb;
        if (size == 32)
          { int32_t ta = (int32_t)a; sa = (int64_t)ta;
            int32_t tb = (int32_t)b; sb = (int64_t)tb;
          }
        else if (size == 64) 
          { sa = (int64_t)a;
            sb = (int64_t)b; 
          }
        else
          { assert(FALSE); }
        fprintf(stderr, "range = {%ld..%ld}\n", sa, sb);
        demand(sa <= sb, "invalid range {a..b}");
        nv1 = (uint64_t)(sb - sa);  /* Number of values. */
      }
    else
      { fprintf(stderr, "range = {%lu..%lu}\n", a, b);
        demand(a <= b, "invalid range {a..b}");
        nv1 = (uint64_t)(b - a);  /* Number of values. */
      }
    fprintf(stderr, "num values = %lu + 1\n", nv1);
        
    /* Allocate and initialize histogram: */
    demand((q >= 1) && ((q - 1) <= nv1), "invalid bin width {q}"); 
    uint64_t nvq = nv1 - (q - 1); /* Number of values minus {q}. */
    demand((nvq % q) == 0, "invalid bin width {q}"); 
    uint64_t nh = nvq/q + 1;  /* Number of histogram bins. */
    demand(nh <= MAX_HIST_BINS, "too many histogram bins");
    uint64_t hs[nh];  /* Histogram of values. */
    for (uint32_t j = 0;  j < nh; j++) { hs[j] = 0; }
    
    /* Collect histogram data: */
    uint64_t nt = (uint32_t)(ntpb*nh);  /* Total number of samples. */
    for (int64_t i = 0; i < nt; i++)
      { uint64_t d; /* Random number minus {a}. */
        if (size == 32)
          { if (sgn)
              { int32_t sa = (int32_t)a;
                int32_t sb = (int32_t)b;
                int32_t x = int32_abrandom(sa, sb);
                demand((x >= sa) && (x <= sb), "int32_abrandom: out of range");
                d = (uint64_t)(x - sa);
              }
            else
              { uint32_t sa = (uint32_t)a;
                uint32_t sb = (uint32_t)b;
                uint32_t x = uint32_abrandom(sa, sb);
                demand((x >= sa) && (x <= sb), "uint32_abrandom: out of range");
                d = (uint64_t)(x - sa);
              }
          }
        else if (size == 64) 
          { if (sgn)
              { int64_t sa = (int64_t)a;
                int64_t sb = (int64_t)b; 
                int64_t x = int64_abrandom(sa, sb);
                demand((x >= sa) && (x <= sb), "int64_abrandom: out of range");
                d = (uint64_t)(x - sa);
              }
            else
              { uint64_t sa = (uint64_t)a;
                uint64_t sb = (uint64_t)b; 
                uint64_t x = uint64_abrandom(sa, sb);
                demand((x >= sa) && (x <= sb), "uint64_abrandom: out of range");
                d = (uint64_t)(x - sa);
              }
          }
        else
          { assert(FALSE); }
        
        uint32_t j = (uint32_t)(d / q); /* Bin number. */
        assert(j < nh);
        hs[j]++;
      }

    /* Write histogram: */
    FILE *wrhs = open_test_file(xtype, size, "abrandom", "qt", xtest); 
    for (uint32_t j = 0;  j < nh; j++) 
      { double r = ((double)hs[j])/((double)ntpb);
        fprintf(wrhs, "%03d %12ld %12.8f\n", j, hs[j], r);
      }
    fclose(wrhs);
  }

void test_drandom(uint32_t ntpb, char *xtest)
  { fprintf(stderr, "Checking {drandom}...\n");
    srandom(19501129);
    uint32_t nh = 300;  /* Number of histogram bins. */
    uint64_t hs[nh];         /* Histogram of values. */
    for (uint32_t j = 0;  j < nh; j++) { hs[j] = 0; }
    
    uint64_t nt = (uint64_t)(ntpb*nh);  /* Total number of samples. */
    uint64_t nzero = 0; /* Number of times that zero appeared. */
    for (int64_t i = 0; i < nt; i++)
      { double x = drandom();
        assert(x >= 0.0);
        assert(x < 1.0);
        if (x == 0) { nzero++; }
        uint32_t j = (uint32_t)floor(x * (double)nh);
        assert(j < nh);
        hs[j]++;
      }
    if (nzero > 0)
      { fprintf(stderr, "the result of {drandom} was zero %lu times\n", nzero);
         demand(FALSE, "that is extremely unlikely");
      }

    FILE *wrhs = open_test_file(NULL, 0, "drandom", "qt", xtest); 
    for (uint32_t j = 0;  j < nh; j++) 
      { double xlo = ((double)j)/((double)nh);
        double xhi = ((double)j+1)/((double)nh);
        double r = ((double)hs[j])/((double)ntpb);
        fprintf(wrhs, "%03d %12lu %12.8f %12.8f %12.8f\n", j, hs[j], xlo, xhi, r);
      }
    fclose(wrhs);
  }

void test_dgaussrand(uint32_t ntpb, char *xtest)
  { fprintf(stderr, "Checking {dgaussrand}...\n");
    srandom(19501129);
    uint32_t nh = 300;  /* Number of histogram bins. */
    uint64_t hs[nh];    /* Histogram of values. */
    for (uint32_t j = 0;  j < nh; j++) { hs[j] = 0; }
    
    double xmin = -6.0; /* Low end of histogram range. */
    double xmax = +6.0; /* High end of histogram range. */
    uint64_t nt = (uint64_t)(ntpb*nh);  /* Total number of samples. */
    uint64_t nzero = 0; /* Number of times that zero appeared. */
    uint64_t nhuge = 0; /* Number of times that huge values appeared. */
    double avg_exp = 0.0; /* Expected average. */
    double dev_exp = 1.0; /* Expected deviation. */
    double sum_x = 0;  /* Sum of values. */
    double sum_dx2 = 0; /* Sum of squared deviations. */
    for (int64_t i = 0; i < nt; i++)
      { double x = dgaussrand();
        sum_x += x;
        double dx = x - avg_exp;
        sum_dx2 += dx*dx;
        if (x == 0) { nzero++; }
        if ((x < xmin) || (x > xmax)) { nhuge++; }
        uint32_t j = (uint32_t)floor(((x - xmin)/(xmax - xmin))*(double)nh);
        if (j < 0) { j = 0; }
        if (j >= nh) { j = nh-1; }
        hs[j]++;
      }
    if ((nzero > 0) || (nhuge > 0))
      { fprintf(stderr, "the result of {dgaussrand} was zero %lu times", nzero);
        fprintf(stderr, " and outside [ %.3f _ %.3f ] %lu times\n", xmin, xmax, nhuge);
        demand(FALSE, "that is extremely unlikely!");
      }

    /* Check average and deviation: */
    double avg_obs = sum_x/((double)nt); /* Observed average. */
    double dev_obs = sqrt(sum_dx2/((double)nt)); /* Observed deviation. */
    fprintf(stderr, "{dgaussrand}: avg = %.3f dev = %.3f", avg_obs, dev_obs);
    fprintf(stderr, "  should be avg = %.3f dev = %.3f\n", avg_exp, dev_exp);
    
    FILE *wrhs = open_test_file(NULL, 0, "dgaussrand", "qt", xtest); 
    for (uint32_t j = 0;  j < nh; j++) 
      { double xlo = xmin + ((double)j)/((double)nh)*(xmax - xmin);
        double xhi = xmin + ((double)j+1)/((double)nh)*(xmax - xmin);
        double r = ((double)hs[j])/((double)nt)/(xhi-xlo);
        fprintf(wrhs, "%03d %12lu %12.8f %12.8f %12.8f\n", j, hs[j], xlo, xhi, r);
      }
    fclose(wrhs);
  }

void test_dloggaussrand(double avg, double dev, uint32_t ntpb, char *xtest)
  { fprintf(stderr, "Checking {dloggaussrand(%.8f,%8f)}...\n", avg, dev);
    srandom(19501129);
    uint32_t nh = 300;  /* Number of histogram bins. */
    uint64_t hs[nh];    /* Histogram of values. */
    for (uint32_t j = 0;  j < nh; j++) { hs[j] = 0; }
    
    /* Convert {avg,dev} of {x} to the avg and dev of {z=log(x)}: */
    double rdv = dev/avg;   /* Deviation relative to the mean. */
    double avg_z = log(avg/hypot(1.0, rdv));
    double dev_z = sqrt(log(1.0 + rdv*rdv));

    double xmin = exp(avg_z - 6.0*dev_z); /* Low end of histogram range. */
    double xmax = exp(avg_z + 6.0*dev_z); /* High end of histogram range. */
    fprintf(stderr, "{dloggaussrand}: range = [ %.6f _ %.6f ]\n", xmin, xmax);
    uint64_t nt = (uint64_t)(ntpb*nh);  /* Total number of samples. */
    uint64_t nhuge = 0; /* Number of times that huge values appeared. */
    uint64_t ntiny = 0; /* Number of times that tiny values appeared. */
    double sum_x = 0;  /* Sum of values. */
    double sum_dx2 = 0; /* Sum of squared deviations. */
    for (int64_t i = 0; i < nt; i++)
      { double x = dloggaussrand(avg, dev);
        assert(x > 0.0);
        sum_x += x;
        double dx = x - avg;
        sum_dx2 += dx*dx;
        if (x > xmax) { nhuge++; }
        if (x < xmin) { nhuge++; }
        uint32_t j = (uint32_t)floor(((x - xmin)/(xmax - xmin))*(double)nh);
        if (j < 0) { j = 0; }
        if (j >= nh) { j = nh-1; }
        hs[j]++;
      }
    if ((ntiny > 0) || (nhuge > 0))
      { fprintf(stderr, "the result of {dloggaussrand} was");
        fprintf(stderr, " smaller than %.3f %lu times", xmin, ntiny);
        fprintf(stderr, " and bigger than %.3f %lu times\n", xmax, nhuge);
        demand(FALSE, "that is extremely unlikely");
      }

    /* Check average and deviation: */
    double avg_obs = sum_x/((double)nt); /* Observed average. */
    double dev_obs = sqrt(sum_dx2/((double)nt)); /* Observed deviation. */
    fprintf(stderr, "{dloggaussrand}: avg = %.3f dev = %.3f", avg_obs, dev_obs);
    fprintf(stderr, "  should be avg = %.3f dev = %.3f\n", avg, dev);

    FILE *wrhs = open_test_file(NULL, 0, "dloggaussrand", "qt", xtest); 
    for (uint32_t j = 0;  j < nh; j++) 
      { double xlo = xmin + ((double)j)/((double)nh)*(xmax - xmin);
        double xhi = xmin + ((double)j+1)/((double)nh)*(xmax - xmin);
        double r = ((double)hs[j])/((double)nt)/(xhi-xlo);
        fprintf(wrhs, "%03d %12lu %12.8f %12.8f %12.8f\n", j, hs[j], xlo, xhi, r);
      }
    fclose(wrhs);
  }

void test_uint64_choose(uint64_t n, size_t k, uint64_t nt, char *xtest)
  { 
    fprintf(stderr, "Checking {uint64_choose}...\n");
    srandom(19501129);
    /* Histogram to count choices: */
    uint32_t nh = ((n < 1000000) && (n <= nt/10) ? (uint32_t)n : 0);
    uint32_t *hist = (nh == 0 ? NULL : notnull(malloc(nh*sizeof(uint32_t)), "no mem"));
    for (uint32_t k = 0;  k < nh; k++) { hist[k] = 0; }
    
    uint64_t *perm = NULL;
    for (int64_t i = 0; i < nt; i++)
      { perm = uint64_choose(n, k, perm);
        /* Check order and range: */
        for (uint32_t i = 0;  i < k; i++)
          { assert(perm[i] < n);
            if (i > 0) { assert(perm[i] > perm[i-1]); }
            if ((i < nh) && (hist != NULL)) { hist[perm[i]]++; }
          }
      }
    if ((nh > 0) && (hist != NULL))
      { /* Check average and deviation of histogram entries: */
        double sum_hk = 0.0;
        for (uint32_t k = 0;  k < nh; k++) 
          { double hk = (double)hist[k]; sum_hk += hk; }
        double avg_h = sum_hk/nh;
        double sum_dk2 = 0.0;
        for (uint32_t k = 0;  k < nh; k++) 
          { double dk = (double)hist[k] - avg_h; sum_dk2 += dk*dk; }
        double dev_h = sqrt(sum_dk2/nh);
        fprintf(stderr, "avg = %.3f dev = %.3f", avg_h, dev_h);
        /* Compute expected average and deviation of histogram entries: */
        double p1 = ((double)k)/((double)n); /* Avg of bin value after 1 try. */
        double v1 = p1*(1-p1);               /* Var of bin value after 1 try. */
        double avg_h_exp = ((double)nt)*p1;
        double dev_h_exp = sqrt(((double)nt)*v1);
        fprintf(stderr, "  should be avg = %.3f dev = %.3f\n", avg_h_exp, dev_h_exp);
      }
        
    if (perm != NULL) { free(perm); }
    if (hist != NULL) { free(hist); }
  }

void test_random_perm(uint32_t n, uint64_t nt, char *xtest)
  { 
    fprintf(stderr, "Checking {random_perm}...\n");
    srandom(19501129);
    
    /* Histogram {hist[i]} is the sum of {perm[i]} for each {i}: */
    uint32_t nh = ((n < 1000000) && (n <= nt/10) ? n : 0);
    uint64_t *hist = (nh == 0 ? NULL : talloc(nh, uint64_t));
    for (uint32_t k = 0;  k < nh; k++) { hist[k] = 0; }
    
    uint32_t *perm = NULL;
    for (uint64_t i = 0; i < nt; i++)
      { perm = random_perm(n, perm);
        /* Check order and range: */
        for (uint32_t i = 0;  i < n; i++)
          { assert(perm[i] < n);
            if ((i < nh) && (hist != NULL)) { hist[i] += perm[i]; }
          }
      }
    if ((nh > 0) && (hist != NULL))
      { /* Check average and deviation of histogram entries: */
        double sum_hk = 0.0;
        for (uint32_t k = 0;  k < nh; k++) 
          { double hk = (double)hist[k]; sum_hk += hk; }
        double avg_h = sum_hk/nh;
        double sum_dk2 = 0.0;
        for (uint32_t k = 0;  k < nh; k++) 
          { double dk = (double)hist[k] - avg_h; sum_dk2 += dk*dk; }
        double dev_h = sqrt(sum_dk2/nh);
        fprintf(stderr, "avg = %.3f dev = %.3f", avg_h, dev_h);
        /* Compute expected average and deviation of histogram entries: */
        double a1 = ((double)n-1)/2;  /* Avg of bin value after 1 try. */
        double v1 = ((double)n-1)*((double)n+1)/12; /* Var of bin value after 1 try. */
        /* !!! NOT RIGHT !!! */
        double avg_h_exp = ((double)nt)*a1;
        double dev_h_exp = sqrt(((double)nt)*v1);
        fprintf(stderr, "  should be avg = %.3f dev = %.3f\n", avg_h_exp, dev_h_exp);
      }
        
    if (perm != NULL) { free(perm); }
    if (hist != NULL) { free(hist); }
  }

FILE *open_test_file(char *xtype, uint32_t size, char *func, char *xend, char *xtest)
  { char *fname = NULL;
    if (xtype != NULL)
      { fname = jsprintf("out/%s%d_%s_%s_%s.his", xtype, size, func, xend, xtest); }
    else
      { fname = jsprintf("out/%s_%s_%s.his", func, xend, xtest); }
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    return wr;
  }
