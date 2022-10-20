#define PROG_NAME "nmsim_test_003_write_read"
#define PROG_DESC "basic tests of {limnmism} write, read, and compare procs"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:34:44 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-03-19"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <fget.h>
#include <affirm.h>
#include <jsmath.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

void do_tests(double prec);
  /* Calls {test_write_read_compare} with {prec} twice 
    to write stuff to disk "out/file0_{prec}.txt", read it back, 
    compare and write it again to "out/file1_{prec}.txt". */

void test_write_read_compare(FILE *rd, FILE *wr, double prec, int32_t phase);
  /* If {phase} is 0, writes a bunch of values to {wr}, integer
    and double, with various options. The {rd} parameter is 
    ignored.
    
    If {phase} is 1, reads the same values from {rd},
    and compares them to the original values.  Then
    writes to {wr} the values that were read from {rd}. */
  
void test_int64(FILE *rd, FILE *wr, int32_t phase, int64_t v);
  /* Same as {test_write_read_compare}, but for the single value {v}. */
  
void test_double(FILE *rd, FILE *wr, int32_t phase, double v, double prec);
  /* Same as {test_write_read_compare}, but for the single value {v}
    with all eight {sgn,fudge_0,fudge_1} options. */

void test_double_case
  ( FILE *rd, FILE *wr, 
    int32_t phase,
    double v,
    double prec,
    bool_t sgn,
    bool_t fudge_0,
    bool_t fudge_1
  );
  /* Same as {test_double}, but for a specific combination
    of {sgn,fudge_0,fudge_1} options. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    do_tests(0.6);
    do_tests(0.06);
    do_tests(0.00000006);
     
    return 0;
  }
   
void do_tests(double prec)
  {
    /* Write a bunch of data to "out/file0_{prec}.txt": */
    
    char *fname0 = NULL;
    asprintf(&fname0, "out/file0_%.12f.txt", prec);
    
    FILE *wr0 = open_write(fname0, TRUE);
    test_write_read_compare(NULL, wr0, prec, 0);
    fclose(wr0);
    
    /* Read it back, compare, and write it again to "out/file1_{prec}.txt": */
    
    char *fname1 = NULL;
    asprintf(&fname1, "out/file1_%.12f.txt", prec);
    
    FILE *rd0 = open_read(fname0, TRUE);
    FILE *wr1 = open_write(fname1, TRUE);
    test_write_read_compare(rd0, wr1, prec, 1);
    fclose(rd0);
    fclose(wr1);
    
    free(fname0);
    free(fname1);
  }
    
void test_write_read_compare(FILE *rd, FILE *wr, double prec, int32_t phase)
  {
    test_int64(rd, wr, phase, -12345678910111213L);
    test_int64(rd, wr, phase, -123456);
    test_int64(rd, wr, phase, -10);
    test_int64(rd, wr, phase, -2);
    test_int64(rd, wr, phase, 0);
    test_int64(rd, wr, phase, +2);
    test_int64(rd, wr, phase, +10);
    test_int64(rd, wr, phase, +123456);
    test_int64(rd, wr, phase, +12345678910111213L);
    
    test_double(rd, wr, phase, NAN,                     prec);
    test_double(rd, wr, phase, -INF,                    prec);
    test_double(rd, wr, phase, -10.0,                   prec);
    test_double(rd, wr, phase, -2.0,                    prec);
    test_double(rd, wr, phase, -10.0/7.0,               prec);
    test_double(rd, wr, phase, -1.0000000000001,        prec);
    test_double(rd, wr, phase, -1.0000000000000,        prec);
    test_double(rd, wr, phase, -0.9999999999999,        prec);
    test_double(rd, wr, phase, -0.5,                    prec);
    test_double(rd, wr, phase, -1.0/700.0,              prec);
    test_double(rd, wr, phase, -0.00012,                prec);
    test_double(rd, wr, phase, -0.00000123456789,       prec);
    test_double(rd, wr, phase, -0.0000000000001,        prec);
    test_double(rd, wr, phase, 0.0,                     prec);
    test_double(rd, wr, phase, +0.0000000000001,        prec);
    test_double(rd, wr, phase, +0.00000123456789,       prec);
    test_double(rd, wr, phase, +0.00012,                prec);
    test_double(rd, wr, phase, +1.0/700.0,              prec);
    test_double(rd, wr, phase, +0.5,                    prec);
    test_double(rd, wr, phase, +0.9999999999999,        prec);
    test_double(rd, wr, phase, +1.0000000000000,        prec);
    test_double(rd, wr, phase, +1.0000000000001,        prec);
    test_double(rd, wr, phase, +10.0/7.0,               prec);
    test_double(rd, wr, phase, +2.0,                    prec);
    test_double(rd, wr, phase, +10.0,                   prec);
    test_double(rd, wr, phase, +INF,                    prec);
  }
  
void test_int64(FILE *rd, FILE *wr, int32_t phase, int64_t v)
  {
    char *lab = NULL; /* Label for compare (the value as in {printf}). */
    asprintf(&lab, "%ld", v);
    if (phase == 0)
      { fprintf(wr, "   ");
        fprintf(wr, "%ld", v);
        fprintf(wr, "\n\n");
      }
    else if (phase == 1)
      { int64_t v_read = nmsim_read_int64_value(rd, lab, INT64_MIN, INT64_MAX);
        fget_eol(rd); fget_eol(rd);
        nmsim_compare_int64_param(lab, v_read, v);
        fprintf(wr, "   ");
        fprintf(wr, "%ld", v_read);
        fprintf(wr, "\n\n");
      }
    free(lab);
  }

void test_double(FILE *rd, FILE *wr, int32_t phase, double v, double prec)
  {
    for (int32_t sgn = 0; sgn <= 1; sgn++)
      { for (int32_t fudge_0 = 0; fudge_0 <= 1; fudge_0++)
          { for (int32_t fudge_1 = 0; fudge_1 <= 1; fudge_1++)
              { test_double_case(rd, wr, phase, v, prec, (sgn != 0), (fudge_0 != 0), (fudge_1 != 0)); }
          }
      }
      if (phase == 0)
        { 
          fprintf(wr, "\n");
        }
      else if (phase == 1) 
        { fget_eol(rd);
          fprintf(wr, "\n");
        }
  }

void test_double_case
  ( FILE *rd, FILE *wr, 
    int32_t phase,
    double v,
    double prec,
    bool_t sgn,
    bool_t fudge_0,
    bool_t fudge_1
  )
  {
    char *lab = NULL; /* Label for compare (the value in "E" format). */
    asprintf(&lab, "%.17e(%d,%.12f,%c,%c,%c)", v, phase, prec, "FT"[sgn], "FT"[fudge_0], "FT"[fudge_1]);

    if (phase == 0)
      { fprintf(wr, "   ");
        nmsim_write_double_value(wr, v, prec, sgn, fudge_0, fudge_1);
        fprintf(wr, "\n");
      }
    else if (phase == 1)
      { double_t v_read = nmsim_read_double_value(rd, lab, -INF, +INF);
        fget_eol(rd); 
        /* Error should be less than {prec}, but fudging may double it: */
        double cmp_prec = prec; /* Precision to use in comparisons. */
        if (fudge_0 && (fabs(v) <= prec)) 
          { if (v == 0.0) 
              { cmp_prec = 0; }
            else
              { cmp_prec = 2*prec; }
          }
        if (fudge_1 && (fabs(1.0 - fabs(v)) <= prec)) 
          { if (fabs(v) == 1.0) 
              { cmp_prec = 0.0; }
            else
              { cmp_prec = 2*prec; }
          }
        nmsim_compare_double_param(lab, v_read, v, cmp_prec, fudge_0, fudge_1);
        fprintf(wr, "   ");
        nmsim_write_double_value(wr, v_read, prec, sgn, fudge_0, fudge_1);
        fprintf(wr, "\n");
      }
    free(lab);
  }
