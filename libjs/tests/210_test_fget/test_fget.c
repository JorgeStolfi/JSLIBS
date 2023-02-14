#define PROG_NAME "test_fget"
#define PROG_DESC "test of {fget.h}, {fget_data.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-02-12 07:52:39 by stolfi */ 
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
#include <fget_data.h>

double zrandom(void); 
  /* A random double with random magnitude. */
  
void check_doubles(double x_wr, char *fmt, double y_rd);
  /* Compares the value {y_rd} with the the result of converting 
    the value {x_wr} to string with the format {fmt}
    and then back to {double} with {strtod}.  Fails with error
    if they don't match. */

#define bug(FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

void test_fget_skip_spaces(bool_t verbose);
void test_fget_double(bool_t verbose);
void test_fget_chars(bool_t verbose);
void test_fget_int32_t(bool_t verbose);
void test_fget_data(bool_t verbose);

int32_t main (int32_t argn, char **argv)
  { srandom(19501129);
    int32_t nt = 3; /* Number of test passes. */
    test_fget_skip_spaces(TRUE);
    test_fget_chars(TRUE);
    for (int32_t it = 0; it < nt; it++)
      { bool_t verbose = (it == 0); 
        test_fget_int32_t(verbose);
        test_fget_double(verbose);
        test_fget_data(verbose);
      }
    return 0;
  }

void test_fget_skip_spaces(bool_t verbose)
  { if (verbose) { fprintf(stderr, "=== Checking {fget_skip_spaces} and related stuff... ===\n"); }

    assert(fget_is_space('\040'));
    assert(fget_is_space('\011'));
    assert(fget_is_space('\240'));

    assert(fget_is_formatting_char('\040'));
    assert(fget_is_formatting_char('\011'));
    assert(fget_is_formatting_char('\012'));
    assert(fget_is_formatting_char('\013'));
    assert(fget_is_formatting_char('\014'));
    assert(fget_is_formatting_char('\015'));
    assert(fget_is_formatting_char('\240'));
  
    char chc = '\377';
    char *str = "\377";
    char chs = (*str);
    int32_t ich = (int32_t)chc;
    
    fprintf(stderr, "chc = %+d\n", chc);
    fprintf(stderr, "chs = %+d\n", chs);
    fprintf(stderr, "ich = %+d\n", ich);
    fprintf(stderr, "EOF = %+d\n", EOF);
    assert(! fget_is_space('\377'));
    assert(! fget_is_formatting_char('\377'));

    /* Write three lines of stuff to disk: */
    char cmtc = '#'; /* Allowed comment chars. */
    char* fname = "out/test.txt";
    FILE* wr = open_write(fname, verbose);
    fprintf(wr, " \t \t \t \240 \n");
    fprintf(wr, " \t \t \t \240 FOO \n");
    fprintf(wr, " \t \t \t \240 %c FOO \n", cmtc);
    fclose(wr);
    
    /* Read back with {fget_skip_spaces} etc: */
    FILE* rd = open_read(fname, verbose);
    
    /* First line should be all spaces: */
    fprintf(stderr, "calling {fget_skip_spaces}...\n");
    fget_skip_spaces(rd);
    fprintf(stderr, "checking next char...");
    int32_t r = fgetc(rd);
    demand(r != EOF, " unexpected EOF");
    char ch = (char)r;
    fprintf(stderr, " got '%c' = \\%03o\n", ch, ch);
    demand(! fget_is_space(ch), "should have skipped it");
    demand(ch != '\240', "should have skipped NBSP");
    demand(ch == '\n', "should have skipped to eol");
    ungetc(r, rd);
    fprintf(stderr, "calling {fget_eol}...\n");
    fget_eol(rd); /* Skip the newline. */
    
    /* Second line should be bunch of blanks plus something: */
    fprintf(stderr, "calling {fget_skip_spaces}...\n");
    fget_skip_spaces(rd);
    fprintf(stderr, "checking next char...");
    r = fgetc(rd);
    demand(r != EOF, " unexpected EOF");
    ch = (char)r;
    fprintf(stderr, " got '%c' = \\%03o\n", ch, ch);
    demand(! fget_is_space(ch), "should have skipped it");
    demand(ch != '\240', "should have skipped NBSP");
    demand(ch == 'F', "should have skipped to 'F' ");
    ungetc(r, rd);
    fprintf(stderr, "calling {fget_skip_to_eol}...\n");
    fget_skip_to_eol(rd); /* Should consume the newline newline. */
    
    /* Third line should be bunch of blanks plus the comment char {cmtc}: */
    fprintf(stderr, "calling {fget_skip_spaces}...\n");
    fget_skip_spaces(rd);
    fprintf(stderr, "checking next char...");
    r = fgetc(rd);
    demand(r != EOF, " unexpected EOF");
    ch = (char)r;
    fprintf(stderr, " got '%c' = \\%03o\n", ch, ch);
    demand(! fget_is_space(ch), "should have skipped it");
    demand(ch != '\240', "should have skipped NBSP");
    demand(ch == cmtc, "should have skipped to {cmtc}");
    ungetc(r, rd);
    fprintf(stderr, "calling {fget_test_comment_or_eol}...\n");
    bool_t ok = fget_test_comment_or_eol(rd, cmtc); /* Should consume the newline newline. */
    demand(ok, "should have found a comment + newline");

    /* Next should be EOF: */
    r = fgetc(rd);
    demand(r == EOF, "expected EOF, not found");
    fclose(rd);
  }
  
void test_fget_chars(bool_t verbose)
  { if (verbose) { fprintf(stderr, "=== Checking {fget_char} etc... ===\n"); }

    /* Write a test file: */
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
    while (! fget_test_eof(rd))
      { int32_t c = fget_char(rd);
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

void test_fget_double(bool_t verbose)
  { if (verbose) { fprintf(stderr, "=== Checking {fget_double}... ===\n"); }
    int32_t nt = 20;
    double x[nt]; /* Numbers to write and read back. */
    int32_t i, k;
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
    int32_t nfmt = 5; /* Number of copies of each value. */
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
    FILE* wr = open_write(fname, verbose);
    for (i = 0; i < nt; i++)
      { if (verbose) { fprintf(stderr, "x[%2d] = %+26.16e\n", i, x[i]); }
        for (k = 0; k < nfmt; k++)
          { fprintf(wr, " ");
            fprintf(wr, fmt[k], x[i]);
          }
        fprintf(wr, "\n");
      }
    fclose(wr);
    
    /* Read back with {fget_double} and check: */
    FILE* rd = open_read(fname, verbose);
    for (i = 0; i < nt; i++)
      { for (k = 0; k < nfmt; k++)
          { double yi = fget_double(rd);
            check_doubles(x[i], fmt[k], yi);
          }
        fget_eol(rd);
      }
    fclose(rd);
  }

void test_fget_int32_t(bool_t verbose)
  { if (verbose) { fprintf(stderr, "=== Checking {fget_int64,fget_uint64}... ===\n"); }

    /* !!! Should check other bases !!! */

    int32_t nt = 20;
    uint64_t xu[nt]; /* Unsigned numbers to write and read back. */
    int64_t  xs[nt]; /* Signed numbers to write and read back. */
    int32_t i, k;
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
    int32_t nfmt = 5; /* Number of copies of each value. */
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
    FILE* wr = open_write(fname, verbose);
    for (i = 0; i < nt; i++)
      { if (verbose) { fprintf(stderr, "xu[%2d] = %lu  xs[%2d] = %ld\n", i, xu[i], i, xs[i]); }
        for (k = 0; k < nfmt; k++)
          { fprintf(wr, fmtu[k], xu[i]);
            fprintf(wr, " ");
            fprintf(wr, fmts[k], xs[i]);
            fprintf(wr, "\n");
          }
        fprintf(wr, "\n");
      }
    fclose(wr);
    
    /* Read back with {fget_uint32,fget_int32} and check: */
    FILE* rd = open_read(fname, verbose);
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

void test_fget_data(bool_t verbose)
  { if (verbose) { fprintf(stderr, "=== Checking {fget_data_fields}... ===\n"); }
  
    int32_t nt = 20;
    double   xd[nt]; /* Float numbers to write and read back. */
    char*    xa[nt]; /* Strings to write and read back. */
    int32_t i = 0;
    xd[i] = 0.0;         i++;
    xd[i] = +17.0;       i++;
    xd[i] = -15.0;       i++;
    xd[i] = 0.1;         i++;
    xd[i] = 0.7e-40;     i++;
    xd[i] = DBL_MAX;     i++;
    xd[i] = -INF;        i++;
    xd[i] = NAN;         i++;
    xd[i] = 1.0/DBL_MAX; i++;
    while (i < nt) { xd[i] = zrandom(); i++; }

    for (int32_t j = 0; j < nt; j++)
      { char *lab = NULL;
        asprintf(&lab, "[%03d]", j);
        xa[j] = lab;
      }
    
    /* Choose the output formats for doubles: */
    int32_t nfmt = 5; /* Number of copies of each value. */
    char* fmt[nfmt];
    int32_t k = 0;
    fmt[k] = "%0.16g"; k++;
    fmt[k] = "%26.16e"; k++;
    fmt[k] = "%+38.16e"; k++;
    fmt[k] = "%032.40f"; k++;
    fmt[k] = "%0.20f"; k++;
    assert(k == nfmt);
    
    /*Set up the data type vector: */
    int32_t nf = 3*nt; /* Number of fields per line: */
    int8_t type[nf]; /* Defines type of each field. */
    for (int32_t kf = 0; kf < nf; kf++) { type[kf] = 0; }
    fget_data_set_field_type(-1, 2, FALSE, nf, type); /* Should be ignored. */
    for (int32_t kt = 0; kt < nt; kt++)
      { fget_data_set_field_type(3*kt + 1, 1, FALSE, nf, type);
        fget_data_set_field_type(3*kt + 2, 2, FALSE, nf, type);
        fget_data_set_field_type(3*kt + 2, 2, TRUE, nf, type); /* Should be OK. */
      }
    
    /* Write {xd[0..nt], xa[0..nt]} to disk in various formats: */
    char cmtc = '#'; /* Allowed comment chars. */
    char* fname = "out/test.txt";
    FILE* wr = open_write(fname, verbose);
    fprintf(wr, " \t \t \t \240 \n");
    int32_t nd = nfmt; /* Number of data records. */
    for (int32_t kd = 0; kd < nd; kd++)
      { if (drandom() < 0.25) { fprintf(wr, " \240 # chances of rain are 25%%\n"); }
        if (drandom() < 0.25) { fprintf(wr, " \t \n"); }
        if (drandom() < 0.25) { fprintf(wr, "      "); }
        char *fmtk = fmt[kd % nfmt]; 
        for (int32_t kf = 0; kf < nf; kf++)
          { fprintf(wr, (drandom() < 0.25 ? "\t" : " "));
            if (type[kf] == 2)
              { fprintf(wr, fmtk, xd[kf/3]); } 
            else
              { fprintf(wr, "%s", xa[kf/3]); } 
          }
        if (drandom() < 0.25) 
          { fprintf(wr, " # bla bla bla."); }
        fprintf(wr, "\n");
        if (drandom() < 0.25) { fprintf(wr, " \t \n"); }
      }
    fclose(wr);
    
    /* Read back with {fget_data_fields} and check: */
    FILE* rd = open_read(fname, verbose);
    double num[nf];
    char *alf[nf];
    for (int32_t kd = 0; kd < nd; kd++)
      { bool_t ok = fget_data_fields(rd, cmtc, nf, type, alf, num);
        demand(ok, "unexpected EOF");
        char *fmtk = fmt[kd % nfmt]; 
        for (int32_t kf = 0; kf < nf; kf++)
          { /* Check data fields that should NOT be read: */
            if (type[kf] != 1)
              { demand(alf[kf] == NULL, "read as alpha a field that should have been skipped"); } 
            if (type[kf] != 2)
              { demand(isnan(num[kf]), "read as double a field that should have been skipped"); }
            
            /* Check that fields that SHOULD be read: */
            if (type[kf] == 1)
              { demand(strcmp(alf[kf], xa[kf/3]) == 0, "mismatched alpha field readback"); }
            else if (type[kf] == 2)
              { double xk = xd[kf/3];
                check_doubles(xk, fmtk, num[kf]);
              }
          }
      }
    bool_t nok = fget_data_fields(rd, cmtc, nf, type, alf, num);
    demand(! nok, "expected EOF, not seen");
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

void check_doubles(double x_wr, char *fmt, double y_rd)
  { 
    char *x_fm = NULL;
    asprintf(&x_fm, fmt, x_wr);
    char *rest = NULL;
    double x_ex = strtod(x_fm, &rest); /* Expected read back value. */
    assert((rest != NULL) && ((*rest) == '\000'));
    bool_t bug = (isnan(x_ex) != isnan(y_rd));
    if ((! isnan(x_ex)) && (! isnan(y_rd))) { bug |= (x_ex != y_rd); }
    if (bug)
      { fprintf(stderr, "** mismatched double field readback\n");
        fprintf(stderr, " written = %24.16e fmt = \"%s\" read = %24.16e expected = %24.16e\n", x_wr, fmt, y_rd, x_ex);
        demand(FALSE, "aborted");
      }
    free(x_fm);
  }
