#define PROG_NAME "test_interpolate"
#define PROG_DESC "tests the {pst_interpolate} routines"
#define PROG_VERS "1.0"

/* Copyright © 2025 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2025-02-20 08:42:26 by stolfi */

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

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int32_t main(int32_t argc, char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    int32_t N = 200;      /* Num of test points. */
    double vP[N], wP[N];  /* Data values and weights. */
    double vH[N+1];       /* Expected interpolated values. */
    
    auto void fill_points(int32_t ip0, int32_t np);
      /* Fills the data points {vP[ip],wP[ip]}, for {ip = ip0+k} and {k}
        in {0..np-1}. In most cases, the weights {wP[ip]} will be
        randomly chosen in {[0_1]}, and the value {vP[ip]} will be
        a frequency-modulated sinusoid. Occasionally there
        will be gaps of one or more points. where {wP[ip]} is 0 and
        {vP[ip]} is {NAN}.
        
        Also fills {vH[ip]}, for {ip = ip0+k} and {k} in {0..np},
        with the same function evaluated at half a step before {vP[ip]}. */
        
    auto double func(double z);
      /* The function that defines the point value {vPk[ip]} for the value
        of {z = ip/(np-1)}. */
    
    fill_points(0, N);
    
    char *fname_data = jsprintf("%s/data.txt", o->outDir);
    FILE *wr_data = open_write(fname_data, TRUE);

    for (int32_t k = 0; k < N; k++)
      { double xk = k;
        double vk = vP[k];
        double wk = wP[k];
        
        if (wk == 0)
          { fprintf(wr_data, "\n"); }
        else
          { fprintf(wr_data, "%8.4f  %8.4f %8.6f\n", xk, vk, wk); }
      }
    fclose(wr_data);

    char *fname_interp = jsprintf("%s/interp.txt", o->outDir);
    FILE *wr_interp = open_write(fname_interp, TRUE);
    
    for (int32_t k = 0; k <= N; k++)
      { double vm = (k <= 1 ? NAN : vP[k-2]);
        double wm = (k <= 1 ? 0.0 : wP[k-2]);
        
        double x0 = k-1;
        double v0 = (k == 0 ? NAN : vP[k-1]);
        double w0 = (k == 0 ? 0.0 : wP[k-1]);
        
        double x1 = k;
        double v1 = (k == N ? NAN : vP[k]);
        double w1 = (k == N ? 0.0 : wP[k]);
        
        double vp = (k >= N-1 ? NAN : vP[k+1]);
        double wp = (k >= N-1 ? 0.0 : wP[k+1]);
        
        double xR = (x0 + x1)/2;
        double vR, wR;
        pst_interpolate_four_values(vm, wm, v0, w0, v1, w1, vp, wp, &vR, &wR);
        
        if (wR == 0)
          { fprintf(wr_interp, "\n"); }
        else
          { fprintf(wr_interp, "%10.4f  %12.8f %12.8f  %12.8f %12.8f\n", xR, vR, wR, vH[k], vR - vH[k]); }
      }
      
    fclose(wr_interp);
    return 0;
    
    double func(double z)
      { double t = 3*z*(1 + 1.0*z); 
        return 3.5 + 2.0*sin(2*M_PI*t);  
      }
    
    void fill_points(int32_t ip0, int32_t np)
      { int32_t nskip = 0; /* Points to skip. */
        int32_t ixgap = 0; /* Gaps generated so far, modulo 3. */
        for (int32_t k = 0; k < np; k++)
          { double zPk = ((double)k)/((double)np -1);
            double vPk, wPk;
            if ((k <= 1) || (k >= np-2) || (nskip > 0))
              { vPk = NAN; wPk = 0;
                if (nskip > 0) { nskip--; }
              }
            else
              { vPk = func(zPk);
                wPk = dabrandom(0, 1);
              }
            /* Just in case: */
            if (wPk == 0) { vPk = NAN; }
            /* Save point: */
            int32_t ip = ip0 + k;
            vP[ip] = vPk; wP[ip] = wPk;
            /* Define random bursts of faults: */
            if ((nskip == 0) && (drandom() < 0.10))
              { nskip = ixgap + 1;
                ixgap = (ixgap + 1) % 3;
              }
          }

        for (int32_t k = 0; k <= np; k++)
          { double zHk = ((double)k-0.5)/((double)np -1);
            double vHk = func(zHk);
            /* Save point: */
            int32_t ip = ip0 + k;
            vH[ip] = vHk;
          }
      }
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
