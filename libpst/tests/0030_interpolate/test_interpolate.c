#define PROG_NAME "test_interpolate"
#define PROG_DESC "tests the {pst_interpolate} routines"
#define PROG_VERS "1.0"

/* Copyright © 2025 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2025-03-16 01:21:20 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -outDir {OUT_DIR}" " \\\n" \
  "    " argparser_help_info_HELP
  
#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Generates a list of points with random weights, integer abscissas, and ordinates that follow linear, quadratic, and cubic polynomials, with undefined gaps.  Then uses {pst_interpolate_four_values} to interpolate those points at half-integer abscissas.  Writes everything out to a text file \"{OUT_DIR}/points.txt\" in a gnuplot-friendly format.\n" \
  "AUTHOR\n" \
  "  Created 2025-02-19 by Jorge Stolfi, UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  By J. Stolfi, IC-UNICAMP unless otherwise noted.\n" \
  "\n" \
  "    2025-02-19 Created."

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

#include <argparser.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <jsrandom.h>
#include <jsmath.h>

#include <pst_interpolate.h>

typedef struct options_t
  { char* outDir;            /* Output directory. */
  } options_t;

/* INTERNAL PROTOTYPES */
    
void tint_define_data_and_exp(int32_t N, double vP[], double wP[], double vH_exp[]);
  /* Fills the data points {vP[k],wP[k]} for {k} in {0..N-1}. In most
    cases, the weights {wP[k]} will be randomly chosen in {[0_1]}, and
    the value {vP[k]} will be a frequency-modulated sinusoid.
    Occasionally there will be gaps of one or more points. where {wP[k]}
    is 0 and {vP[k]} is {NAN}.

    Also fills {vH_exp[k]}, for {k} in {0..N}, with the same function
    evaluated at half a step before {vP[k]} (but witout {NAN}s). */

void tint_write_data(char *outDir, int32_t N, double vP[], double wP[]);
  /* Writes {k}, {vP[k]} and {wP[k]} for {k} in {0..N-1}
    to file "{outDir}/data.txt". */

void tint_write_iterp(char *outDir, int32_t N, double vH_exp[], double vH_cmp[], double wH_cmp[]);
  /* Writes to file "{outDir}/interp.txt" {N+1} lines with the interpolated data.
    Each line has "{k-0.5} {vH_exp[k]} {vH_cmp[k]} {wH_cmp[k]} {vE[k]}"
    where {k} ranges in {0..N} and {vE[k]} is {vH_cmp[k] - vH_exp[k]}.
    If {wH_cmp[k]} is 0 (and thus {vH_cmp[k]} is {NAN}), omits the last three values. */

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int32_t main(int32_t argc, char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    int32_t N = 160;      /* Num of test points. */
    double vP[N], wP[N];  /* Function values and weights at integer args. */
    double vH_exp[N+1];   /* Ideal function values at half-integer args. */
    
    tint_define_data_and_exp(N, vP, wP, vH_exp);
    tint_write_data(o->outDir, N, vP, wP);

    double vH_cmp[N+1];   /* Interpolated values at half-integer args. */
    double wH_cmp[N+1];   /* Weights of interpolated values at half-integer args. */
    
    auto void fetch(int32_t j, double *vR_P, double *wR_P);
      /* Returns {vP[j]} and {wP[j]}, if {j} is in {0..N-1};
        otherwise returns {NAN,0}. */

    for (int32_t j = 0; j <= N; j++)
      { double vR, wR;

        int32_t ja, jb, n, m;
        pst_interpolate_select_data(j-1, fetch, &ja, &jb, &n, &m);
        assert((n >= 0) && (n <= 4));
        assert((m >= 0) && (m <= n));

        /* Collect values to interpolate: */
        double vS[n], wS[n];
        for (int32_t k = 0; k < n; k++)
          { vS[k] = vP[ja + k];
            wS[k] = wP[ja + k];
          }
        
        pst_interpolate_values((uint32_t)n, vS, wS, (uint32_t)m, &vR, &wR);
        vH_cmp[j] = vR;
        wH_cmp[j] = wR;
        
        /* Decide whether weight {wR} is expected to be zero: */
        bool_t expZero = (n == 0);

        bool_t bug = FALSE;
        if (wR == 0)
          { if (! expZero) 
              { fprintf(stderr, "** got zero weight out, expected nonzero\n"); bug = TRUE; }
          }
        else
          { if (expZero)
              { fprintf(stderr, "** got nonzero weight out, expected zero\n"); 
                fprintf(stderr, "  vR = %24.16e  wR = %24.16e\n", vR, wR);
                bug = TRUE;
              }
          }
        if (bug)
          { fprintf(stderr, "  n = %d  m= %d\n", n, m);
            for (int32_t k = 0; k <= n; k++)
              { if (k == m) { fprintf(stderr, "  .....\n"); }
                if (k < n) { fprintf(stderr, "  vS[%d] = %24.16e  wS[%d] = %24.16e\n", k, vS[k], k, wS[k]); }
              }
          }
      }
      
    tint_write_iterp(o->outDir, N, vH_exp, vH_cmp, wH_cmp);
    
    return 0;
    
    void fetch(int32_t j, double *vR_P, double *wR_P)
      { double vR, wR;
        if ((j >= 0) && (j < N))
          { vR = vP[j]; wR = wP[j]; }
        else
          { vR = NAN; wR = 0; }
        assert((isfinite(vR)) == ((isfinite(wR)) && (wR > 0)));
        (*vR_P) = vR;
        (*wR_P) = wR;
      }
  }

void tint_define_data_and_exp(int32_t N, double vP[], double wP[], double vH_exp[])
  { 
    auto void plop(int32_t i);
    
    auto void skip(int32_t i);
    
    auto double func(double z);
      /* Computes a frequency-varying sinusoid. Assumes the 
        {z} argument varies in {[0 _ 1]}. */
    
    int32_t k = 0;
    for (int32_t j = 0; j < 6; j++) { plop(k); k++; }
    int32_t nr = 5;
    for (int32_t b = 0; b < 6; b++)
      { for (int32_t r = 1; (r <= nr) && (k + r + 2 <= N); r++)
          { /* Write a block with 1 gap, {r} points, 1 gap */
            skip(k); k++;
            for (int32_t i = 0; i < r; i++) { plop(k); k++; }
            skip(k); k++;
          }
      }
    while (k < N) { plop(k); k++; }

    for (int32_t k = 0; k <= N; k++)
      { double zHk = ((double)k - 0.5)/((double)N - 1);
        double vHk = func(zHk);
        vH_exp[k] = vHk;
      }

    return;
    
    void plop(int32_t i)
      { double zPk = ((double)i)/((double)N - 1);
        double vPk = func(zPk); 
        double wPk = dabrandom(0, 1);
        if (wPk == 0) { vPk = NAN; }
        vP[i] = vPk; wP[i] = wPk;
      }
    
    void skip(int32_t i)
      { double vPk = NAN; 
        double wPk = 0;
        vP[i] = vPk; wP[i] = wPk;
      }
    
    double func(double z)
      { double t = 3*z*(1 + 1.0*z); 
        return 3.5 + 2.0*sin(2*M_PI*t);  
      }
  }
    
void tint_write_data(char *outDir, int32_t N, double vP[], double wP[])
  { 
    char *fname = jsprintf("%s/data.txt", outDir);
    FILE *wr = open_write(fname, TRUE);

    for (int32_t k = 0; k < N; k++)
      { double xk = k;
        double vk = vP[k];
        double wk = wP[k];
        
        if (wk == 0)
          { fprintf(wr, "\n"); }
        else
          { fprintf(wr, "%8.4f  %8.4f %8.6f\n", xk, vk, wk); }
      }
    fclose(wr);
    free(fname);
  }
    
void tint_write_iterp(char *outDir, int32_t N, double vH_exp[], double vH_cmp[], double wH_cmp[])
  {
    char *fname = jsprintf("%s/interp.txt", outDir);
    FILE *wr = open_write(fname, TRUE);
    for (int32_t j = 0; j <= N; j++)
      { double xR = ((double)j) - 0.5;
        double vF = vH_exp[j];
        double vR = vH_cmp[j];
        double wR = wH_cmp[j];
        if (wR == 0)
          { fprintf(wr, "%10.4f  %12.8f\n", xR, vF); }
        else
          { fprintf(wr, "%10.4f  %12.8f  %12.8f %12.8f  %12.8f\n", xR, vF, vR, wR, vR - vF); }
      }
    fclose(wr);
    free(fname);
  }

options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
   
    options_t *o = (options_t *)malloc(sizeof(options_t));
    
    argparser_get_keyword(pp, "-outDir"); 
    o->outDir = argparser_get_next(pp);

    argparser_finish(pp);
    
    return o;
  }
  

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 2005 by the State University of Campinas (UNICAMP).
**
** Created on 2005-08-15 by Jorge Stolfi, IC-UNICAMP.       
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
