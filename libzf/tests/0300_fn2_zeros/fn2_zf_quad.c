/* See fn2_zf_quad.h */
/* Last edited on 2024-12-05 10:41:28 by stolfi */

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <malloc.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <epswr.h>
#include <bool.h>
#include <jsfile.h>

#include <fn2_zf_quad.h>
#include <fn2_zf_flt.h>
#include <fn2_zf_plot.h>

/*** PROTOTYPES FOR INTERNAL ROUTINES ***/

void fn2_zf_quad_plot_aa_zeros
  ( epswr_figure_t *eps,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int32_t n
  );
  /* Plots roots of {f(xf,yf) = 0}, using {eval_aa} to evaluate 
    it with affine arithmetic.
    The input intervals are {n}×{n} cells in {xd}×{yd}. */
  
void fn2_zf_quad_aa_zeros_aux
  ( epswr_figure_t *eps,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int32_t n,
    int32_t *nevals,
    int32_t *ncells,
    int32_t ixlo, int32_t mx,
    int32_t iylo, int32_t my
  );
  /* Called by {fn2_zf_quad_plot_aa_zeros} */

void fn2_zf_quad_plot_ia_zeros
  ( epswr_figure_t *eps,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int32_t n
  );
  /* Plots zeros of {fv(xv,yv) = 0}, using ordinary interval arithmetic
    on {n×n} equal boxes in {xd×yd}. */

void fn2_zf_quad_ia_zeros_aux
  ( epswr_figure_t *eps,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int32_t n,
    int32_t *nevals,
    int32_t *ncells,
    int32_t ixlo, int32_t mx,
    int32_t iylo, int32_t my
  );
  /* Called by {fn2_zf_quad_plot_ia_zeros} */

/*** IMPLEMENTATIONS ***/

void fn2_zf_quad_plots(
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
      fn2_zf_quad_plot_ia_zeros(eps, eval_ia, xd, yd, n);
      fn2_zf_plot_end_figure(eps, eval_fp, xd, yd, m, n, n);
    }
    
    /* Affine arithmetic plot */
    { epswr_figure_t *eps = fn2_zf_plot_new_figure(prefix, "aa", xd, yd, 7, title, parm_string, "Affine arithmetic");
      epswr_set_pen(eps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
      fn2_zf_quad_plot_aa_zeros(eps, eval_aa, xd, yd, n);
      fn2_zf_plot_end_figure(eps, eval_fp, xd, yd, m, n, n);
    }
  }

void fn2_zf_quad_plot_aa_zeros(
    epswr_figure_t *eps,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int32_t n
  )
  { int32_t nevals = 0;
    int32_t ncells = 0;

    fprintf(stderr, "\nenter fn2_zf_quad_plot_aa_zeros...\n");

    epswr_comment(eps, "Plot of AA-arithmetic zeros");

    fn2_zf_quad_aa_zeros_aux(
      eps, eval_aa, xd, yd, n, &nevals, &ncells,
      0, n, 0, n
    );

    fprintf(stderr, " done\n");
    
    fn2_zf_plot_stat_caption(eps, "%d %s function evaluations", nevals, "AA");
    fn2_zf_plot_stat_caption(eps, "%d cells retained", ncells, NULL);

    fprintf(stderr, "exit fn2_zf_quad_plot_aa_zeros.\n");
  }

void fn2_zf_quad_aa_zeros_aux(
    epswr_figure_t *eps,
    eval_aa_t *eval_aa,
    Interval xd,
    Interval yd,
    int32_t n,
    int32_t *nevals,
    int32_t *ncells,
    int32_t ixlo, int32_t mx,
    int32_t iylo, int32_t my
  )
  { Interval xv, yv, zv;
    AAP xf, yf, zf;
    double gray = 0.75;

    MemP frame = aa_top();

    fprintf(stderr, "(");

    ROUND_DOWN;
    xv.lo = xd.lo + ((xd.hi - xd.lo)*(Float)ixlo)/(Float)n;
    yv.lo = yd.lo + ((yd.hi - yd.lo)*(Float)iylo)/(Float)n;

    ROUND_UP;
    xv.hi = xd.lo + ((xd.hi - xd.lo)*(Float)(ixlo+mx))/(Float)n;
    yv.hi = yd.lo + ((yd.hi - yd.lo)*(Float)(iylo+my))/(Float)n;

    xf = aa_from_interval(xv);
    yf = aa_from_interval(yv);

    zf = eval_aa(xf, yf);
    zv = aa_range(zf);
    (*nevals)++;

    aa_flush(frame);

    ROUND_NEAR;
    if ((zv.lo > Zero) || (zv.hi < Zero))
      { epswr_rectangle(eps, xv.lo,xv.hi, yv.lo,yv.hi, FALSE, TRUE); }
    else if ((mx == 1) && (my == 1))
      { (*ncells)++;
        epswr_set_fill_color(eps, gray, gray, gray);
        epswr_rectangle(eps,  xv.lo,xv.hi, yv.lo,yv.hi, TRUE, TRUE); 
      }
    else if (mx >= my)
      { int32_t mx2 = mx / 2;
        fn2_zf_quad_aa_zeros_aux(
          eps, eval_aa, xd, yd, n, nevals, ncells,
          ixlo, mx2, iylo, my
        );
        fn2_zf_quad_aa_zeros_aux(
          eps, eval_aa, xd, yd, n, nevals, ncells,
          ixlo + mx2, mx - mx2, iylo, my
        );
      }
    else /* if (my > mx) */
      { int32_t my2 = my / 2;
        fn2_zf_quad_aa_zeros_aux(
          eps, eval_aa, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo, my2
        );
        fn2_zf_quad_aa_zeros_aux(
          eps, eval_aa, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo + my2, my - my2
        );
      }

    fprintf(stderr, ")");
  }

void fn2_zf_quad_plot_ia_zeros
  ( epswr_figure_t *eps,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int32_t n
  )
  { int32_t nevals = 0;
    int32_t ncells = 0;

    fprintf(stderr, "\nenter fn2_zf_quad_plot_ia_zeros...\n");

    epswr_comment(eps, "Plot of interval-arithmetic zeros");

    fn2_zf_quad_ia_zeros_aux(
      eps, eval_ia, xd, yd, n, &nevals, &ncells,
      0, n, 0, n
    );

    fprintf(stderr, "\n");
    
    fn2_zf_plot_stat_caption(eps, "%d %s function evaluations", nevals, "IA");
    fn2_zf_plot_stat_caption(eps, "%d cells retained", ncells, NULL);

    fprintf(stderr, "exit fn2_zf_quad_plot_ia_zeros\n");
  }


void fn2_zf_quad_ia_zeros_aux(
    epswr_figure_t *eps,
    eval_ia_t *eval_ia,
    Interval xd,
    Interval yd,
    int32_t n,
    int32_t *nevals,
    int32_t *ncells,
    int32_t ixlo, int32_t mx,
    int32_t iylo, int32_t my
  )
  { Interval xv, yv, zv;
    double gray = 0.75;

    fprintf(stderr, "(");

    ROUND_DOWN;
    xv.lo = xd.lo + ((xd.hi - xd.lo)*(Float)ixlo)/(Float)n;
    yv.lo = yd.lo + ((yd.hi - yd.lo)*(Float)iylo)/(Float)n;

    ROUND_UP;
    xv.hi = xd.lo + ((xd.hi - xd.lo)*(Float)(ixlo+mx))/(Float)n;
    yv.hi = yd.lo + ((yd.hi - yd.lo)*(Float)(iylo+my))/(Float)n;

    zv = eval_ia(xv, yv);
    (*nevals)++;

    ROUND_NEAR;
    if ((zv.lo > Zero) || (zv.hi < Zero))
      { epswr_rectangle(eps, xv.lo,xv.hi, yv.lo,yv.hi, FALSE, TRUE); }
    else if ((mx == 1) && (my == 1))
      { epswr_set_fill_color(eps, gray, gray, gray);
        epswr_rectangle(eps, xv.lo,xv.hi, yv.lo,yv.hi, TRUE, TRUE);
        (*ncells)++;
      }
    else if (mx >= my)
      { int32_t mx2 = mx / 2;
        fn2_zf_quad_ia_zeros_aux(
          eps, eval_ia, xd, yd, n, nevals, ncells,
          ixlo, mx2, iylo, my
        );
        fn2_zf_quad_ia_zeros_aux(
          eps, eval_ia, xd, yd, n, nevals, ncells,
          ixlo + mx2, mx - mx2, iylo, my
        );
      }
    else /* if (my > mx) */
      { int32_t my2 = my / 2;
        fn2_zf_quad_ia_zeros_aux(
          eps, eval_ia, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo, my2
        );
        fn2_zf_quad_ia_zeros_aux(
          eps, eval_ia, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo + my2, my - my2
        );
      }

    fprintf(stderr, ")");
  }

