/* See fn2_zf_quad.h */
/* Last edited on 2009-01-06 03:14:20 by stolfi */

#include "fn2_zf_quad.h"
#include "fn2_zf_flt.h"

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <pswr.h>
#include <bool.h>
#include <jsfile.h>
#include <jsstring.h>

#include <math.h>
#include <stdio.h>
#include <malloc.h>

/*** PROTOTYPES FOR INTERNAL ROUTINES ***/

void fn2_zf_quad_plot_aa_zeros(
    PSStream *ps,
    int epsformat,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int n
  );
  /* 
    Plots roots of {f(xf,yf) = 0}, using {eval_aa} to 
    evaluate it with affine arithmetic.
    The input intervals are {n}×{n} cells in {xd}×{yd}. */

char *fn2_zf_quad_format_parms(
    Interval xd,
    Interval yd,
    int n
  );
  /* Formats the arguments into a string */
  
void fn2_zf_quad_report_stats(
    PSStream *ps,
    int epsformat, 
    char *arith,    /* "IA" or "AA" */
    int nevals, 
    int ncells
  );
  /* Prints evaluation report, and adds caption lines if not $epsformat$ */
  
void fn2_zf_quad_aa_zeros_aux(
    PSStream *ps,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int n,
    int *nevals,
    int *ncells,
    int ixlo, int mx,
    int iylo, int my
  );
  /* Called by fn2_zf_quad_plot_aa_zeros */

void fn2_zf_quad_plot_ia_zeros(
    PSStream *ps,
    int epsformat,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int n
  );
  /* Plots zeros of fv(xv,yv) = 0, using ordinary interval arithmetic.
    The input inervals are nxn cells in xd x yd. */

void fn2_zf_quad_ia_zeros_aux(
    PSStream *ps,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int n,
    int *nevals,
    int *ncells,
    int ixlo, int mx,
    int iylo, int my
  );
  /* Called by fn2_zf_quad_plot_ia_zeros */

/*** IMPLEMENTATIONS ***/

void fn2_zf_quad_plots(
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
  { 
    /* Compute the plot scale so that the largest dimension is 6 inches: */
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
    char *parmstr = fn2_zf_quad_format_parms(xd, yd, n);
    
    /*** Ordinary interval arithmetic plot ***/
    pswr_new_canvas(ps, "ia");
    pswr_new_picture(ps, xd.lo, xd.hi, yd.lo, yd.hi);
    pswr_set_grid(ps, n, n);
    
    if (! epsformat)
      { pswr_add_caption(ps, title, 0.0);
        pswr_add_caption(ps, "", 0.0);
        pswr_add_caption(ps, parmstr, 0.0);
        pswr_add_caption(ps, "", 0.0);
        pswr_add_caption(ps, "Ordinary interval arithmetic", 0.0);
      }

    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    fn2_zf_quad_plot_ia_zeros(ps, epsformat, eval_ia, xd, yd, n);

    pswr_set_pen(ps, 0.5,0.0,0.0, 0.10, 0.0, 0.0);
    fn2_zf_flt_plot(ps, eval_fp, xd, yd, m);

    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
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
    fn2_zf_quad_plot_aa_zeros(ps, epsformat, eval_aa, xd, yd, n);

    pswr_set_pen(ps, 0.5,0.0,0.0, 0.10, 0.0, 0.0);
    fn2_zf_flt_plot(ps, eval_fp, xd, yd, m);

    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    pswr_frame(ps);

    pswr_close_stream(ps);
  }

char *fn2_zf_quad_format_parms(
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

void fn2_zf_quad_plot_aa_zeros(
    PSStream *ps,
    int epsformat,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int n
  )
  { int nevals = 0;
    int ncells = 0;

    fprintf(stderr, "\nenter fn2_zf_quad_plot_aa_zeros...\n");

    pswr_comment(ps, "Plot of AA-arithmetic zeros");

    fn2_zf_quad_aa_zeros_aux(
      ps, eval_aa, xd, yd, n, &nevals, &ncells,
      0, n, 0, n
    );

    fprintf(stderr, " done\n");
    
    fn2_zf_quad_report_stats(ps, epsformat, "AA", nevals, ncells);

    fprintf(stderr, "exit fn2_zf_quad_plot_aa_zeros.\n");
  }

void fn2_zf_quad_aa_zeros_aux(
    PSStream *ps,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int n,
    int *nevals,
    int *ncells,
    int ixlo, int mx,
    int iylo, int my
  )
  { Interval xv, yv, zv;
    AAP xf, yf, zf;
    double gray = 0.75;

    MemP frame = aa_top();

    fprintf(stderr, "(");

    ROUND_DOWN;
    xv.lo = xd.lo + ((xd.hi - xd.lo)*ixlo)/n;
    yv.lo = yd.lo + ((yd.hi - yd.lo)*iylo)/n;

    ROUND_UP;
    xv.hi = xd.lo + ((xd.hi - xd.lo)*(ixlo+mx))/n;
    yv.hi = yd.lo + ((yd.hi - yd.lo)*(iylo+my))/n;

    xf = aa_from_interval(xv);
    yf = aa_from_interval(yv);

    zf = eval_aa(xf, yf);
    zv = aa_range(zf);
    (*nevals)++;

    aa_flush(frame);

    ROUND_NEAR;
    if ((zv.lo > Zero) || (zv.hi < Zero))
      { pswr_rectangle(ps, xv.lo, xv.hi, yv.lo, yv.hi, FALSE, TRUE); }
    else if ((mx == 1) && (my == 1))
      { (*ncells)++;
        pswr_set_fill_color(ps, gray, gray, gray);
        pswr_rectangle(ps,  xv.lo, xv.hi, yv.lo, yv.hi, TRUE, TRUE); 
      }
    else if (mx >= my)
      { int mx2 = mx / 2;
        fn2_zf_quad_aa_zeros_aux(
          ps, eval_aa, xd, yd, n, nevals, ncells,
          ixlo, mx2, iylo, my
        );
        fn2_zf_quad_aa_zeros_aux(
          ps, eval_aa, xd, yd, n, nevals, ncells,
          ixlo + mx2, mx - mx2, iylo, my
        );
      }
    else /* if (my > mx) */
      { int my2 = my / 2;
        fn2_zf_quad_aa_zeros_aux(
          ps, eval_aa, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo, my2
        );
        fn2_zf_quad_aa_zeros_aux(
          ps, eval_aa, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo + my2, my - my2
        );
      }

    fprintf(stderr, ")");
  }

void fn2_zf_quad_plot_ia_zeros(
    PSStream *ps,
    int epsformat,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int n
  )
  { int nevals = 0;
    int ncells = 0;

    fprintf(stderr, "\nenter fn2_zf_quad_plot_ia_zeros...\n");

    pswr_comment(ps, "Plot of interval-arithmetic zeros");

    fn2_zf_quad_ia_zeros_aux(
      ps, eval_ia, xd, yd, n, &nevals, &ncells,
      0, n, 0, n
    );

    fprintf(stderr, "\n");
    
    fn2_zf_quad_report_stats(ps, epsformat, "IA", nevals, ncells);

    fprintf(stderr, "exit fn2_zf_quad_plot_ia_zeros\n");
  }

void fn2_zf_quad_report_stats(
    PSStream *ps, 
    int epsformat, 
    char *arith, 
    int nevals, 
    int ncells
  )
  { char sevals[80];

    fprintf(stderr, "%d %s function evaluations\n", nevals, arith);
    if (! epsformat)
      { sprintf(sevals, "%d %s function evaluations", nevals, arith);
        pswr_add_caption(ps, sevals, 0.0);
      }

    fprintf(stderr, "%d cells retained\n", ncells);
    if (! epsformat)
      { sprintf(sevals, "%d cells retained", ncells);
        pswr_add_caption(ps, sevals, 0.0);
      }
  }

void fn2_zf_quad_ia_zeros_aux(
    PSStream *ps,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int n,
    int *nevals,
    int *ncells,
    int ixlo, int mx,
    int iylo, int my
  )
  { Interval xv, yv, zv;
    double gray = 0.75;

    fprintf(stderr, "(");

    ROUND_DOWN;
    xv.lo = xd.lo + ((xd.hi - xd.lo)*ixlo)/n;
    yv.lo = yd.lo + ((yd.hi - yd.lo)*iylo)/n;

    ROUND_UP;
    xv.hi = xd.lo + ((xd.hi - xd.lo)*(ixlo+mx))/n;
    yv.hi = yd.lo + ((yd.hi - yd.lo)*(iylo+my))/n;

    zv = eval_ia(xv, yv);
    (*nevals)++;

    ROUND_NEAR;
    if ((zv.lo > Zero) || (zv.hi < Zero))
      { pswr_rectangle(ps, xv.lo, xv.hi, yv.lo, yv.hi, FALSE, TRUE); }
    else if ((mx == 1) && (my == 1))
      { pswr_set_fill_color(ps, gray, gray, gray);
        pswr_rectangle(ps, xv.lo, xv.hi, yv.lo, yv.hi, TRUE, TRUE);
        (*ncells)++;
      }
    else if (mx >= my)
      { int mx2 = mx / 2;
        fn2_zf_quad_ia_zeros_aux(
          ps, eval_ia, xd, yd, n, nevals, ncells,
          ixlo, mx2, iylo, my
        );
        fn2_zf_quad_ia_zeros_aux(
          ps, eval_ia, xd, yd, n, nevals, ncells,
          ixlo + mx2, mx - mx2, iylo, my
        );
      }
    else /* if (my > mx) */
      { int my2 = my / 2;
        fn2_zf_quad_ia_zeros_aux(
          ps, eval_ia, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo, my2
        );
        fn2_zf_quad_ia_zeros_aux(
          ps, eval_ia, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo + my2, my - my2
        );
      }

    fprintf(stderr, ")");
  }

