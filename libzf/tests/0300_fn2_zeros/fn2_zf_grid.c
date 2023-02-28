/* See fn2_zf_grid.h */
/* Last edited on 2023-02-18 10:00:38 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <malloc.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <epswr.h>

#include <fn2_functions.h>
#include <fn2_zf_grid.h>
#include <fn2_zf_flt.h>
#include <fn2_zf_plot.h>

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

char *fn2_zf_grid_format_parms(
    Interval xd,
    Interval yd,
    int32_t n
  );
  /* Formats the arguments into a sctring */
  
void fn2_zf_grid_report_stats(
    epswr_figure_t *eps,
    int32_t epsformat, 
    char *arith,    /* "IA" or "AA" */
    int32_t nevals, 
    int32_t ncells
  );
  /* Prints evaluation report, and adds caption lines if not $epsformat$ */
  
void fn2_zf_grid_plot_aa_zeros
  ( epswr_figure_t *eps,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int32_t n
  );

void fn2_zf_grid_plot_ia_zeros(
    epswr_figure_t *eps,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int32_t n
  );

/*** IMPLEMENTATIONS ***/

void fn2_zf_grid_plots(
    char *prefix,
    char *title,
    eval_fp_t *eval_fp,
    eval_ia_t *eval_ia,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int32_t n,
    int32_t m
  )
  { 
    /* Pack the parameters as a string: */
    char *parm_string = fn2_zf_plot_format_parms(xd, yd, n);
    
    /* Ordinary interval arithmetic plot */
    { epswr_figure_t *eps = fn2_zf_plot_new_figure(prefix, "ia", xd, yd, 7, title, parm_string, "Ordinary interval arithmetic");
      epswr_set_pen(eps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
      fn2_zf_grid_plot_ia_zeros(eps, eval_ia, xd, yd, n);
      epswr_set_pen(eps, 0.25,0.25,0.25, 0.10, 0.0, 0.0);
      epswr_grid_lines(eps, n, n);
      fn2_zf_plot_end_figure(eps, eval_fp, xd, yd, m, n, n);
    }
     
    /* Affine arithmetic plot */
    { epswr_figure_t *eps = fn2_zf_plot_new_figure(prefix, "aa", xd, yd, 7, title, parm_string, "Affine arithmetic");
      epswr_set_pen(eps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
      fn2_zf_grid_plot_aa_zeros(eps, eval_aa, xd, yd, n);
      epswr_set_pen(eps, 0.25,0.25,0.25, 0.10, 0.0, 0.0);
      epswr_grid_lines(eps, n, n);
      fn2_zf_plot_end_figure(eps, eval_fp, xd, yd, m, n, n);
    }
  }

void fn2_zf_grid_plot_aa_zeros(
    epswr_figure_t *eps,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int32_t n
  )
  { Interval xv, yv, zv;
    AAP xf, yf, zf;
    double gray = 0.75;
    int32_t ncells = 0;
    int32_t nevals = 0;

    fprintf(stderr, "\nenter fn2_zf_grid_plot_aa_zeros...\n");

    epswr_comment(eps, "Plot of AA-arithmetic zeros");

    for (int32_t yi=0; yi<n; yi++)
      { for (int32_t xi=0; xi<n; xi++)
          { MemP frame = aa_top();

            ROUND_DOWN;
            xv.lo = xd.lo + ((xd.hi - xd.lo)*(Float)xi)/(Float)n;
            yv.lo = yd.lo + ((yd.hi - yd.lo)*(Float)yi)/(Float)n;

            ROUND_UP;
            xv.hi = xd.lo + ((xd.hi - xd.lo)*(Float)(xi+1))/(Float)n;
            yv.hi = yd.lo + ((yd.hi - yd.lo)*(Float)(yi+1))/(Float)n;

            xf = aa_from_interval(xv);
            yf = aa_from_interval(yv);

            zf = eval_aa(xf, yf);
            zv = aa_range(zf);
            nevals++;

            ROUND_NEAR;
            if ((zv.lo <= Zero) && (zv.hi >= Zero))
              { epswr_set_fill_color(eps, gray, gray, gray);
                epswr_grid_cell(eps, xi, n, yi, n, TRUE, FALSE);
                ncells++;
              }

            aa_flush(frame);
          }
      }
    
    fn2_zf_plot_stat_caption(eps, "%d %s function evaluations", nevals, "AA");
    fn2_zf_plot_stat_caption(eps, "%d unresolved cells", ncells, NULL);

    fprintf(stderr, "exit fn2_zf_grid_plot_aa_zeros.\n");
  }

void fn2_zf_grid_plot_ia_zeros(
    epswr_figure_t *eps,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int32_t n
  )
  { Interval xv, yv, zv;
    double gray = 0.75;
    int32_t ncells = 0;
    int32_t nevals = 0;

    fprintf(stderr, "\nenter fn2_zf_grid_plot_ia_zeros...\n");

    epswr_comment(eps, "Plot of interval-arithmetic zeros");

    for (int32_t yi=0; yi<n; yi++)
      { for (int32_t xi=0; xi<n; xi++)
          { ROUND_DOWN;
            xv.lo = xd.lo + ((xd.hi - xd.lo)*(Float)xi)/(Float)n;
            yv.lo = yd.lo + ((yd.hi - yd.lo)*(Float)yi)/(Float)n;

            ROUND_UP;
            xv.hi = xd.lo + ((xd.hi - xd.lo)*(Float)(xi+1))/(Float)n;
            yv.hi = yd.lo + ((yd.hi - yd.lo)*(Float)(yi+1))/(Float)n;

            zv = eval_ia(xv, yv);
            nevals++;

            ROUND_NEAR;
            if ((zv.lo <= Zero) && (zv.hi >= Zero))
              { epswr_set_fill_color(eps, gray, gray, gray);
                epswr_grid_cell(eps, xi, n, yi, n, TRUE, FALSE);
                ncells++;
              }
          }
      }
    
    fn2_zf_plot_stat_caption(eps, "%d %s function evaluations", nevals, "IA");
    fn2_zf_plot_stat_caption(eps, "%d unresolved cells", ncells, NULL);

    fprintf(stderr, "exit fn2_zf_grid_plot_ia_zeros.\n");
  }

