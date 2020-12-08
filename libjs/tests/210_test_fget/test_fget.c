#define PROG_NAME "test_fget"
#define PROG_DESC "test of {fget.h}"
#define PROG_VERS "1.0"

/* Last edited on 2015-04-04 09:41:49 by stolfilocal */ 
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


int main (int argn, char **argv)
  { srandom(19501129);
    test_fget_chars();
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

#define LOG_DBL_MAX (709.78271289338397)

double zrandom(void)
  { 
    double e = 2*drandom() - 1;
    double m = drandom();
    double r = m*exp(LOG_DBL_MAX*e);
    return r;
  }
