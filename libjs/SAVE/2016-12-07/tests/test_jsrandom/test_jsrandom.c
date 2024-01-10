#define PROG_NAME "test_jsrandom"
#define PROG_DESC "test of {jsmath.h}"
#define PROG_VERS "1.0"

/* Last edited on 2013-12-17 00:57:59 by stolfilocal */ 
/* Created on 2011-09-20 by J. Stolfi, UNICAMP */

#define test_jsrandom_COPYRIGHT \
  "Copyright © 2011  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsmath.h>
#include <jsrandom.h>
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

void test_uint64_random(int nb, int ntpb, char *tag);
  /* Tests {uint64_random}. Writes to disk
    histograms of the high and low {nb} bits, with {ntpd}
    samples per bin on average. */

void test_abrandom(int32_t a, int32_t b, int ntpb, char *tag);
  /* Tests {abrandom} with interval {a..b}. Writes to disk
    a histogram of values with {ntpd} samples per bin on average. */

void test_drandom(int ntpb, char *tag);
  /* Tests {drandom}. Writes to disk a histogram of values
    with {ntpd} samples per bin on average. */

FILE *open_test_file(char *prefix, char *tag);

int main (int argn, char **argv)
  { 
    test_uint64_random(6, 10000, "A");
    
    test_abrandom(-200, +100, 10000, "A");
    test_abrandom(   0,    1, 10000, "B");
    test_drandom(10000, "A");
    // test_dgaussrandom(10000);

    return 0;
  }

void test_uint64_random(int nb, int ntpb, char *tag)
  { fprintf(stderr, "Checking {uint64_random}...\n");
    int nh = (int)ipow(2, nb);  /* Number of histogram bins, a power of 2. */
    int hi[nh], lo[nh]; /* Histograms of high and low bits. */
    int i,j;

    for (j = 0; j < nh; j++) { lo[j] = hi[j] = 0; }

    int nt = ntpb*nh;      /* Total number of samples. */
    for (i = 0; i < nt; i++)
      { uint64_t x = uint64_random();
        int jlo = (int)(x & (uint64_t)(nh - 1));
        lo[jlo]++;
        int jhi = (int)(x >> (64 - nb));
        hi[jhi]++;
      }

    FILE *wrhi = open_test_file("uint64_random-hi", tag); 
    for (j = 0; j < nh; j++) 
      { double r = ((double)hi[j])/((double)ntpb);
        fprintf(wrhi, "%03d %8d %12.8f\n", j, hi[j], r);
      }
    fclose(wrhi);
    FILE *wrlo = open_test_file("uint64_random-lo", tag); 
    for (j = 0; j < nh; j++) 
      { double r = ((double)lo[j])/((double)ntpb);
        fprintf(wrlo, "%03d %8d %12.8f\n", j, lo[j], r);
      }
    fclose(wrlo);
  }

void test_abrandom(int32_t a, int32_t b, int ntpb, char *tag)
  { fprintf(stderr, "Checking {abrandom}...\n");
    uint64_t nh = (uint64_t)(b - a) + 1;  /* Number of histogram bins. */
    int hs[nh];  /* Histogram of values. */
    int i,j;

    for (j = 0; j < nh; j++) { hs[j] = 0; }
    
    int nt = (int)(ntpb*nh);  /* Total number of samples. */
    for (i = 0; i < nt; i++)
      { int32_t x = abrandom(a, b);
        assert(x >= a);
        assert(x <= b);
        int j = x - a; 
        assert(j < nh);
        hs[j]++;
      }

    FILE *wrhs = open_test_file("abrandom", tag); 
    for (j = 0; j < nh; j++) 
      { double r = ((double)hs[j])/((double)ntpb);
        fprintf(wrhs, "%03d %8d %12.8f\n", j, hs[j], r);
      }
    fclose(wrhs);
  }

void test_drandom(int ntpb, char *tag)
  { fprintf(stderr, "Checking {drandom}...\n");
    srandom(19501129);
    uint64_t nh = 300;  /* Number of histogram bins. */
    int hs[nh];         /* Histogram of values. */
    int i,j;

    for (j = 0; j < nh; j++) { hs[j] = 0; }
    
    int nt = (int)(ntpb*nh);  /* Total number of samples. */
    int nzero = 0; /* Number of times that zero appeared. */
    for (i = 0; i < nt; i++)
      { double x = drandom();
        assert(x >= 0.0);
        assert(x < 1.0);
        if (x == 0) { nzero++; }
        int j = (int)floor(x * (double)nh);
        assert(j < nh);
        hs[j]++;
      }
    fprintf(stderr, "the result of {drandom} was zero %d times\n", nzero);
    demand(nzero == 0, "that is extremely unlikely");

    FILE *wrhs = open_test_file("drandom", tag); 
    for (j = 0; j < nh; j++) 
      { double r = ((double)hs[j])/((double)ntpb);
        fprintf(wrhs, "%03d %8d %12.8f\n", j, hs[j], r);
      }
    fclose(wrhs);
  }

FILE *open_test_file(char *prefix, char *tag)
  { char *fname = NULL;
    asprintf(&fname, "out/%s-%s.his", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    return wr;
  }
