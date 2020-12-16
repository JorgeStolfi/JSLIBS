#define PROG_NAME "test_fget"
#define PROG_DESC "test of {fget.h}"
#define PROG_VERS "1.0"

/* Last edited on 2020-12-15 00:21:05 by jstolfi */ 
/* Created on 2007-01-02 by J. Stolfi, UNICAMP */

#define test_jsmath_COPYRIGHT \
  "Copyright © 2013  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <jsfile.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <affirm.h>
#include <bool.h>
#include <fget.h>

double zrandom(void); /* A random double with random magnitude. */

#define bug(FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

void test_fget_double(void);
void test_fget_chars(void);
void test_fget_int(void);

int main (int argn, char **argv)
  { srandom(19501129);
    test_fget_chars();
    test_fget_int();
    test_fget_double();
    return 0;
  }
  
void test_fget_chars(void)
  { /* Write a test file: */
    char* fname = "out/testc.txt";
    FILE* wr = open_write(fname, TRUE);
    fprintf(wr, "1X YZ\n");
    fprintf(wr, "1XYZ \n");
    fprintf(wr, "2  X YZ\n");
    fprintf(wr, "2  XY # Z\n");
    fprintf(wr, "3X  YZ\n");
    fprintf(wr, "3X  YZ # Comment\n");
    fprintf(wr, "3# Comment\n");
    fclose(wr);

    /* Read back with {fget_test_char,fget_skip_spaces,fget_skip_and_test_char} and check: */
    FILE* rd = open_read(fname, TRUE);
    while (! fget_test_char(rd,EOF))
      { int c = fget_char(rd);
        fprintf(stderr, "[%c]", c);
        if (c == '1')
          { /* Just read all characters with {fgetc}. */
            while (! fget_test_char(rd, '\n'))
              { c = fgetc(rd);
                fprintf(stderr, "[%c]", c); 
              }
            fprintf(stderr, "[EOL]\n"); 
          }
        else if (c == '2')
          { /* Skip leading spaces then read all with {fget_char}. */
            fget_skip_spaces(rd);
            while (! fget_skip_and_test_char(rd, '\n')) 
              { c = fget_char(rd);
                fprintf(stderr, "[%c]", c);
              }
            fprintf(stderr, "[EOL]\n");
          }
        else if (c == '3')
          { /* Skip leading spaces then read until comment. */
            while (TRUE)
              { if (fget_skip_and_test_char(rd, '#'))  { fget_skip_to_eol(rd); break; }
                if (fget_skip_and_test_char(rd, '\n')) { break; }
                c = fget_char(rd);
                fprintf(stderr, "[%c]", c);
              }
            fprintf(stderr, "[EOL]\n");
          }
      }
    fclose(rd);
  }

void test_fget_double(void)
  { fprintf(stderr, "Checking {fget_double}...\n");
    int nt = 20;
    double x[nt]; /* Numbers to write and read back. */
    int i, k;
    i = 0;
    x[i] = 0.0;         i++;
    x[i] = +17.0;       i++;
    x[i] = -15.0;       i++;
    x[i] = 0.1;         i++;
    x[i] = 0.7e-40;     i++;
    x[i] = DBL_MAX;     i++;
    x[i] = -INF;        i++;
    x[i] = NAN;         i++;
    x[i] = 1.0/DBL_MAX; i++;
    while (i < nt) { x[i] = zrandom(); i++; }
    
    /* Choose the output formats: */
    int nfmt = 5; /* Number of copies of each value. */
    char* fmt[nfmt];
    k = 0;
    fmt[k] = "%.16g"; k++;
    fmt[k] = "%26.16e"; k++;
    fmt[k] = "%+38.16e"; k++;
    fmt[k] = "%0802.400f"; k++;
    fmt[k] = "%.400f"; k++;
    assert(k == nfmt);
    
    /* Write {x[0..nt]} to disk in various formats: */
    char* fname = "out/test.txt";
    FILE* wr = open_write(fname, TRUE);
    for (i = 0; i < nt; i++)
      { fprintf(stderr, "x[%2d] = %+26.16e\n", i, x[i]);
        for (k = 0; k < nfmt; k++)
          { fprintf(wr, " ");
            fprintf(wr, fmt[k], x[i]);
          }
        fprintf(wr, "\n");
      }
    fclose(wr);
    
    /* Read back with {fget_double} and check: */
    FILE* rd = open_read(fname, TRUE);
    for (i = 0; i < nt; i++)
      { for (k = 0; k < nfmt; k++)
          { double yi = fget_double(rd);
            bool_t bug = (isnan(yi) != isnan(x[i]));
            if ((! isnan(yi)) && (! isnan(x[i]))) { bug |= (yi != x[i]); }
            if (bug)
              { fprintf(stderr, "** fget_double error? written = %26.16e fmt = \"%s\" read = %26.16e\n", x[i], fmt[k], yi); }
          }
        fget_eol(rd);
      }
    fclose(rd);
  }

void test_fget_int(void)
  { fprintf(stderr, "Checking {fget_int64,fget_uint64}...\n");

    /* !!! Should check other bases !!! */

    int nt = 20;
    uint64_t xu[nt]; /* Unsigned numbers to write and read back. */
    int64_t  xs[nt]; /* Signed numbers to write and read back. */
    int i, k;
    i = 0;
    xu[i] = 0;           xs[i] = 0;            i++;
    xu[i] = 17;          xs[i] = +17;          i++;
    xu[i] = 15;          xs[i] = -15;          i++;
    xu[i] = 12345;       xs[i] = -12345;       i++;
    xu[i] = 10000000000; xs[i] = -10000000000; i++;
    xu[i] = UINT64_MAX;  xs[i] = INT64_MAX;    i++;
    xu[i] = 0;           xs[i] = INT64_MIN;    i++;
    while (i < nt) { xu[i] = uint64_random();  xs[i] = int64_random(); i++; }
    
    /* Choose the output formats: */
    int nfmt = 5; /* Number of copies of each value. */
    char *fmtu[nfmt]; char *fmts[nfmt];
    k = 0;
    fmtu[k] = "%lu";      fmts[k] = "%ld";      k++;
    fmtu[k] = "%030lu";   fmts[k] = "%030ld";   k++;
    fmtu[k] = "%030lu";   fmts[k] = "%+030ld";  k++;
    fmtu[k] = "%200lu";   fmts[k] = "%+200ld";  k++;
    fmtu[k] = "%lu";      fmts[k] = "%+ld";     k++;
    assert(k == nfmt);
    
    /* Write {xu[0..nt], xs[0..nt]} to disk in various formats: */
    char* fname = "out/test.txt";
    FILE* wr = open_write(fname, TRUE);
    for (i = 0; i < nt; i++)
      { fprintf(stderr, "xu[%2d] = %lu  xs[%2d] = %ld\n", i, xu[i], i, xs[i]);
        for (k = 0; k < nfmt; k++)
          { fprintf(wr, fmtu[k], xu[i]);
            fprintf(wr, " ");
            fprintf(wr, fmts[k], xs[i]);
            fprintf(wr, "\n");
          }
        fprintf(wr, "\n");
      }
    fclose(wr);
    
    /* Read back with {fget_double} and check: */
    FILE* rd = open_read(fname, TRUE);
    for (i = 0; i < nt; i++)
      { for (k = 0; k < nfmt; k++)
          { uint64_t yui = fget_uint64(rd, 10);
            if (yui != xu[i])
              { fprintf(stderr, "** fget_uint64 error");
                fprintf(stderr, " written = %lu fmt = \"%s\" read = %lu\n", xu[i], fmtu[k], yui);
              }
            int64_t ysi = fget_int64(rd);
            if (ysi != xs[i])
              { fprintf(stderr, "** fget_int64 error");
                fprintf(stderr, " written = %ld fmt = \"%s\" read = %ld\n", xs[i], fmts[k], ysi);
              }
            fget_eol(rd);
          }
        fget_eol(rd);
      }
    fclose(rd);
  }

#define LOG_DBL_MAX (709.78271289338397)

double zrandom(void)
  { 
    double e = 2*drandom() - 1;
    double m = drandom();
    double r = m*exp(LOG_DBL_MAX*e);
    return r;
  }
