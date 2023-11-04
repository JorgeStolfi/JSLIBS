#define PROG_NAME "test_fget"
#define PROG_DESC "test of {fget.h}, {fget_data.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-10-15 03:06:37 by stolfi */ 
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
#include <jsstring.h>
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
void test_fget_int(bool_t verbose);
void test_fget_data(bool_t verbose);

int32_t main (int32_t argn, char **argv)
  { srandom(19501129);
    int32_t nt = 3; /* Number of test passes. */
    test_fget_skip_spaces(TRUE);
    test_fget_chars(TRUE);
    for (int32_t it = 0; it < nt; it++)
      { bool_t verbose = (it == 0); 
        test_fget_int(verbose);
        test_fget_double(verbose);
        test_fget_data(verbose);
      }
    return 0;
  }

void test_fget_skip_spaces(bool_t verbose)
  { if (verbose) { fprintf(stderr, "=== Checking {fget_skip_spaces} and related stuff... ===\n"); }

    assert(fget_is_space('\040'));  assert(fget_is_formatting_char('\040'));
    assert(fget_is_space('\011'));  assert(fget_is_formatting_char('\011'));
    assert(fget_is_space('\240'));  assert(fget_is_formatting_char('\240'));
    
    assert(! fget_is_space('\012'));  assert(fget_is_formatting_char('\012'));
    assert(! fget_is_space('\013'));  assert(fget_is_formatting_char('\013'));
    assert(! fget_is_space('\014'));  assert(fget_is_formatting_char('\014'));
    assert(! fget_is_space('\015'));  assert(fget_is_formatting_char('\015'));
  
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

    auto void check_cmt(bool_t ok_cmp, bool_t ok_exp, char *text_cmp, char *text_exp);
    auto bool_t eq_cmt_text(char *text_cmp, char *text_exp);
    auto void print_cmt_res(char *which, bool_t ok, char *text);
    auto void print_cmt_text(char *text);

    /* Write six lines of stuff to disk: */
    int32_t nlines = 6;
    char *x[nlines];
    x[0] = " \t \t \t \240 ";                
    x[1] = " \t \t \t \240 FOO ";            
    x[2] = " \t \t \t \240 # FOO ";   
    x[3] = " # BAR ";                 
    x[4] = " ";                              
    x[5] = " \t @ FOO ";
    assert(nlines == 6);
    char* fname = "out/test_skip.txt";
    FILE* wr = open_write(fname, verbose);
    for (int32_t i = 0; i < nlines; i++) { fprintf(wr, "%s\n", x[i]); }                  
    fclose(wr);
    
    /* Read back with {fget_skip_spaces} etc: */
    char cmtc = '#'; /* Comment char. */
    FILE* rd = open_read(fname, verbose);
    
    fprintf(stderr, "\n-- line 0 = \"%s\":\n", escapify(x[0]));
    /* Should be all spaces: */
    fprintf(stderr, "calling {fget_skip_spaces}...\n");
    fget_skip_spaces(rd);
    fprintf(stderr, "checking next char...");
    int32_t r0 = fgetc(rd);
    demand(r0 != EOF, " unexpected EOF");
    char ch0 = (char)r0;
    fprintf(stderr, " got '%c' = \\%03o\n", ch0, ch0);
    demand(! fget_is_space(ch0), "should have skipped it");
    demand(ch0 != '\240', "should have skipped NBSP");
    demand(ch0 == '\n', "should have skipped to eol");
    ungetc(r0, rd);
    fprintf(stderr, "calling {fget_eol}...\n");
    fget_eol(rd); /* Skip the end-of-line. */
    
    fprintf(stderr, "\n-- line 1 = \"%s\":\n", escapify(x[1]));
    /* Should be bunch of blanks plus something: */
    fprintf(stderr, "calling {fget_skip_spaces}...\n");
    fget_skip_spaces(rd);
    fprintf(stderr, "checking next char...");
    int32_t r1 = fgetc(rd);
    demand(r1 != EOF, " unexpected EOF");
    char ch1 = (char)r1;
    fprintf(stderr, " got '%c' = \\%03o\n", ch1, ch1);
    demand(! fget_is_space(ch1), "should have skipped it");
    demand(ch1 != '\240', "should have skipped NBSP");
    demand(ch1 == 'F', "should have skipped to 'F' ");
    ungetc(r1, rd);
    fprintf(stderr, "calling {fget_skip_to_eol}...\n");
    fget_skip_to_eol(rd); /* Should consume the end-of-line. */
    
    fprintf(stderr, "\n-- line 2 = \"%s\":\n", escapify(x[2]));
    /* Should be bunch of blanks plus the comment "{cmtc} FOO ": */
    fprintf(stderr, "calling {fget_skip_spaces}...\n");
    fget_skip_spaces(rd);
    fprintf(stderr, "checking next char...");
    int32_t r2 = fgetc(rd);
    demand(r2 != EOF, " unexpected EOF");
    char ch2 = (char)r2;
    fprintf(stderr, " got '%c' = \\%03o\n", ch2, ch2);
    demand(! fget_is_space(ch2), "should have skipped it");
    demand(ch2 != '\240', "should have skipped NBSP");
    demand(ch2 == cmtc, "should have skipped to {cmtc}");
    ungetc(r2, rd);
    fprintf(stderr, "calling {fget_test_comment_or_eol}...\n");
    char *text2 = "UNCHANGED";
    bool_t ok2 = fget_test_comment_or_eol(rd, cmtc, &text2); /* Should consume the end-of-line too. */
    check_cmt(ok2, TRUE, text2, " FOO ");

    fprintf(stderr, "\n-- line 3 = \"%s\":\n", escapify(x[3]));
    /* Should be bunch of blanks plus the comment "{cmtc} BAR ": */
    char *text3 = "UNCHANGED";
    fget_comment_or_eol(rd, cmtc, &text3); /* Should consume the end-of-line too. */
    fprintf(stderr, "computed: text = "); print_cmt_text(text3); fprintf(stderr, "\n");
    demand(eq_cmt_text(text3, " BAR "), "did not return the comment text");

    fprintf(stderr, "\n-- line 4 = \"%s\":\n", escapify(x[4]));
    /* Should be just blanks: */
    char *text4 = "UNCHANGED";
    fget_comment_or_eol(rd, cmtc, &text4); /* Should consume the end-of-line too. */
    fprintf(stderr, "computed: text = "); print_cmt_text(text4); fprintf(stderr, "\n");
    demand(eq_cmt_text(text4, NULL), "returned a non-{NULL} comment text");

    fprintf(stderr, "\n-- line 5 = \"%s\":\n", escapify(x[5]));
    /* Should be blanks then a '@': */
    char *text5 = "UNCHANGED";
    bool_t ok5 = fget_test_comment_or_eol(rd, cmtc, &text5); /* Should consume the end-of-line too. */
    check_cmt(ok5, FALSE, text5, "UNCHANGED");
    int32_t r5 = fgetc(rd);
    demand(r5 != EOF, " unexpected EOF");
    char ch5 = (char)r5;
    fprintf(stderr, " got '%c' = \\%03o\n", ch5, ch5);
    demand(ch5 == '@', "was expecting '@'");
    ungetc(r5, rd);
    fprintf(stderr, "calling {fget_skip_to_eol}...\n");
    fget_skip_to_eol(rd); /* Should consume the end-of-line. */

    fprintf(stderr, "\n-- end of file:\n");
    /* Should be at EOF: */
    int32_t rEOF = fgetc(rd);
    demand(rEOF == EOF, "expected EOF, not found");
    fclose(rd);

    return;

    void check_cmt(bool_t ok_cmp, bool_t ok_exp, char *text_cmp, char *text_exp)
      { 
        print_cmt_res("computed", ok_cmp, text_cmp);
        print_cmt_res("expected", ok_exp, text_exp);
        demand(ok_cmp == ok_exp, "wrong test result");
        demand(eq_cmt_text(text_cmp, text_exp), "wrong {text} returned");
      }

    bool_t eq_cmt_text(char *text_cmp, char *text_exp)
      { if ((text_cmp == NULL) && (text_exp == NULL))
          { return TRUE; }
        else if ((text_cmp != NULL) && (text_exp != NULL))
          { return strcmp(text_cmp, text_exp) == 0; }
        else
          { return FALSE; }
      }

    void print_cmt_res(char *which, bool_t ok, char *text)
      { fprintf(stderr, "%s: result = %c text = ", which, "FT"[ok]);
        print_cmt_text(text);
        fprintf(stderr, "\n");
      }

    void print_cmt_text(char *text)
      { if (text == NULL) 
          { fprintf(stderr, "{NULL}"); }
        else if (strcmp(text, "UNCHANGED") == 0)
           { fprintf(stderr, "unchanged"); }
        else
          { fprintf(stderr, "\"%s\"", text); }
      }


  }
  
void test_fget_chars(bool_t verbose)
  { if (verbose) { fprintf(stderr, "=== Checking {fget_char} etc... ===\n"); }

    /* Write a test file: */
    char* fname = "out/test_chars.txt";
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
    char* fname = "out/test_double.txt";
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

void test_fget_int(bool_t verbose)
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
    char* fname = "out/test_int.txt";
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
    FILE *rd = open_read(fname, verbose);
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
  
    int32_t nt = 20; /* Number of data fields of each type. */
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
    char* fname = "out/test_data.txt";
    FILE* wr = open_write(fname, verbose);
    int32_t nd = nfmt; /* Number of data records. */
    int32_t kd_to_kl[nd+1]; /* The file lines which have data records. */
    int32_t nl = 0; /* Number of lines (including blanks): */
    for (int32_t kd = 0; kd < nd; kd++)
      { /* Maybe write some blank or comment lines: */
        if (drandom() < 0.25) { fprintf(wr, " \240 # chances of rain are 25%%\n"); nl++; }
        if (drandom() < 0.25) { fprintf(wr, " \t \n"); nl++; }
        /* Write a data line: */
        kd_to_kl[kd] = nl;
        if (drandom() < 0.25) { fprintf(wr, "      "); }
        char *fmtk = fmt[kd % nfmt]; 
        for (int32_t kf = 0; kf < nf; kf++)
          { /* Print a space or tab: */
            fprintf(wr, (drandom() < 0.25 ? "\t" : " "));
            if (type[kf] == 0)
              { /* Ignored: */ fprintf(wr, "<%s>", xa[kf/3]); } 
            else if (type[kf] == 1)
              { /* Alpha: */ fprintf(wr, "%s", xa[kf/3]); } 
            else if (type[kf] == 2)
              { /* Numeric: */ fprintf(wr, fmtk, xd[kf/3]); } 
            else
              { assert(FALSE); }
          }
        if (drandom() < 0.25) 
          { fprintf(wr, " # bla bla bla."); }
        fprintf(wr, "\n"); nl++;
        if (drandom() < 0.25) { fprintf(wr, " \t \n"); nl++; }
      }
    kd_to_kl[nd] = nl; /* Sentinel. */
    fclose(wr);
    fprintf(stderr, "wrote %d file lines with %d data lines\n", nl, nd);
    
    /* Read back each line as a string: */
    FILE* rd = open_read(fname, verbose);
    char *x[nl];
    for (int32_t kl = 0; kl < nl; kl++)
      { x[kl] = fget_line(rd);
        fprintf(stderr, "-- file line %d = \"%s\"\n", kl, escapify(x[kl]));
      }
    if (! fget_test_eof(rd))
      { fprintf(stderr, "** expected EOF, not seen\n");
        char *xtra = fget_line(rd);
        fprintf(stderr, "-- extra line = \"%s\"\n", escapify(xtra));
        assert(FALSE);
      }
    fclose(rd);

    /* Read back with {fget_data_fields} and check: */
    rd = open_read(fname, verbose);
    double num[nf];
    char *alf[nf];
    int32_t kl = 0;
    for (int32_t kd = 0; kd <= nd; kd++)
      { while ((kl < nl) && (kl < kd_to_kl[kd]))
          { fprintf(stderr, "-- file line %d = \"%s\"\n", kl, escapify(x[kl])); kl++; }
        if (kd >= nd) {break; }
        assert(kl < nl);
        fprintf(stderr, "\n!! file line %d (data record %d) = \"%s\"\n", kl, kd, escapify(x[kl])); kl++;
        bool_t ok = fget_data_fields(rd, cmtc, nf, type, alf, num);
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
