/* See fn2_zf_grid.h */
/* Last edited on 2009-01-06 03:13:57 by stolfi */

#include <fn2_functions.h>
#include <fn2_zf_grid.h>
#include <fn2_zf_flt.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <pswr.h>

#include <math.h>
#include <stdio.h>
#include <malloc.h>

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

char *fn2_zf_grid_format_parms(
    Interval xd,
    Interval yd,
    int n
  );
  /* Formats the arguments into a sctring */
  
void fn2_zf_grid_report_stats(
    PSStream *ps,
    int epsformat, 
    char *arith,    /* "IA" or "AA" */
    int nevals, 
    int ncells
  );
  /* Prints evaluation report, and adds caption lines if not $epsformat$ */
  
void fn2_zf_grid_plot_aa_zeros(PSStream *ps,
    int epsformat,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int n
  );

void fn2_zf_grid_plot_ia_zeros(
    PSStream *ps,
    int epsformat,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int n
  );

/*** IMPLEMENTATIONS ***/

void fn2_zf_grid_plots(
    char *fileprefix,
    int epsformat,
    char *title,
    eval_fp_t *eval_fp,
    eval_ia_t *eval_ia,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int n,
    int m
  )
  { /* Compute the plot scale so that the largest dimension is 6 inches: */
    Float wx = xd.hi - xd.lo;
    Float wy = yd.hi - yd.lo;
    Float wm = (wx > wy ? wx : wy);
    double scale = 6.0 * 72.0 / wm;
    
    /* Compute the canvas dimensions (excluding the margins): */
    double hsize = scale*wx; /* Canvas H size. */
    double vsize = scale*wy; /* Canvas V size. */
    double mrg = 4.0; /* In pt. */

    /* Create the Postscript stream: */
    PSStream *ps = pswr_new_stream(fileprefix, NULL, epsformat, "doc", "letter", hsize+2*mrg, vsize + 2*mrg);
    
    /* Pack the parameters as a string: */
    char *parmstr = fn2_zf_grid_format_parms(xd, yd, n);
    
    /*** Ordinary interval arithmetic plot ***/
    pswr_new_canvas(ps, "ia");
    pswr_new_picture(ps, xd.lo, xd.hi, yd.lo, yd.hi);
    pswr_set_grid(ps, n, n);

    if (! epsformat)
      { pswr_add_caption(ps, title, 0.5);
        pswr_add_caption(ps, "", 0.5);
        pswr_add_caption(ps, parmstr, 0.5);
        pswr_add_caption(ps, "", 0.5);
        pswr_add_caption(ps, "Ordinary interval arithmetic", 0.5);
      }

    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    fn2_zf_grid_plot_ia_zeros(ps, epsformat, eval_ia, xd, yd, n);
    
    pswr_set_pen(ps, 0.5,0.0,0.0, 0.25, 0.0, 0.0);
    fn2_zf_flt_plot(ps, eval_fp, xd, yd, m);
    
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    pswr_grid_lines(ps);

    pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0.0, 0.0);
    pswr_frame(ps);

    /*** Affine arithmetic plot ***/
    pswr_new_canvas(ps, "aa");
    pswr_new_picture(ps, xd.lo, xd.hi, yd.lo, yd.hi);
    pswr_set_grid(ps, n, n);

    if (! epsformat)
      { pswr_add_caption(ps, title, 0.0);
	pswr_add_caption(ps, "", 0.0);
	pswr_add_caption(ps, parmstr, 0.0);
	pswr_add_caption(ps, "", 0.0);
	pswr_add_caption(ps, "Affine arithmetic", 0.0);
      }

    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    fn2_zf_grid_plot_aa_zeros(ps, epsformat, eval_aa, xd, yd, n);
    
    pswr_set_pen(ps, 0.5,0.0,0.0, 0.25, 0.0, 0.0);
    fn2_zf_flt_plot(ps, eval_fp, xd, yd, m);
    
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    pswr_grid_lines(ps);

    pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0.0, 0.0);
    pswr_frame(ps);

    pswr_close_stream(ps);
  }

char *fn2_zf_grid_format_parms(
    Interval xd,
    Interval yd,
    int n
  )
  { char *s = (char *) malloc(100);
    sprintf(s, 
      "x in [%f _ %f]  y in [%f _ %f]  %d intervals",
      xd.lo, xd.hi, yd.lo, yd.hi, n
    );
    return(s);
  }

void fn2_zf_grid_plot_aa_zeros(
    PSStream *ps,
    int epsformat,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int n
  )
  { Interval xv, yv, zv;
    int xi, yi;
    AAP xf, yf, zf;
    double gray = 0.75;
    int ncells = 0;
    char sevals[80];

    fprintf(stderr, "\nenter fn2_zf_grid_plot_aa_zeros...\n");

    pswr_comment(ps, "Plot of AA-arithmetic zeros");

    for (yi=0; yi<n; yi++)
      { for (xi=0; xi<n; xi++)
          { MemP frame = aa_top();

            ROUND_DOWN;
            xv.lo = xd.lo + ((xd.hi - xd.lo)*xi)/n;
            yv.lo = yd.lo + ((yd.hi - yd.lo)*yi)/n;

            ROUND_UP;
            xv.hi = xd.lo + ((xd.hi - xd.lo)*(xi+1))/n;
            yv.hi = yd.lo + ((yd.hi - yd.lo)*(yi+1))/n;

            xf = aa_from_interval(xv);
            yf = aa_from_interval(yv);

            zf = eval_aa(xf, yf);
            zv = aa_range(zf);

            ROUND_NEAR;
            if ((zv.lo <= Zero) && (zv.hi >= Zero))
              { pswr_set_fill_color(ps, gray, gray, gray);
                pswr_grid_cell(ps, xi, yi, TRUE, FALSE);
                ncells++;
              }

            aa_flush(frame);
          }
      }

    fprintf(stderr, " done\n");
    fprintf(stderr, "%d unresolved cells\n", ncells);
    if (! epsformat)
      { sprintf(sevals, "%d unresolved cells", ncells);
        pswr_add_caption(ps, sevals, 0.0);
      }
    fprintf(stderr, "exit fn2_zf_grid_plot_aa_zeros.\n");
  }

void fn2_zf_grid_plot_ia_zeros(
    PSStream *ps,
    int epsformat,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int n
  )
  { Interval xv, yv, zv;
    int xi, yi;
    double gray = 0.75;
    int ncells = 0;
    char sevals[80];

    fprintf(stderr, "\nenter fn2_zf_grid_plot_ia_zeros...\n");

    pswr_comment(ps, "Plot of interval-arithmetic zeros");

    for (yi=0; yi<n; yi++)
      { for (xi=0; xi<n; xi++)
          { ROUND_DOWN;
            xv.lo = xd.lo + ((xd.hi - xd.lo)*xi)/n;
            yv.lo = yd.lo + ((yd.hi - yd.lo)*yi)/n;

            ROUND_UP;
            xv.hi = xd.lo + ((xd.hi - xd.lo)*(xi+1))/n;
            yv.hi = yd.lo + ((yd.hi - yd.lo)*(yi+1))/n;

            zv = eval_ia(xv, yv);

            ROUND_NEAR;
            if ((zv.lo <= Zero) && (zv.hi >= Zero))
              { pswr_set_fill_color(ps, gray, gray, gray);
                pswr_grid_cell(ps, xi, yi, TRUE, FALSE);
                ncells++;
              }
          }
      }

    fprintf(stderr, " done\n");
    fprintf(stderr, "%d unresolved cells\n", ncells);
    if (! epsformat)
      { sprintf(sevals, "%d unresolved cells", ncells);
        pswr_add_caption(ps, sevals, 0.0);
      }
    fprintf(stderr, "exit fn2_zf_grid_plot_ia_zeros.\n");
  }

