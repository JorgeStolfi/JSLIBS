/* See fn2_zf_plot.h */
/* Last edited on 2023-02-18 09:57:32 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <fn2_functions.h>
#include <epswr.h>
#include <bool.h>

#include "fn2_zf_plot.h"
#include "fn2_zf_flt.h"

epswr_figure_t *fn2_zf_plot_new_figure
  ( char *prefix, 
    char *arith_tag, 
    Interval xd,
    Interval yd,
    int32_t ncap,
    char *title, 
    char *parm_string, 
    char *arith_title
  )
  {
    fprintf(stderr, "=== %s (%s) plot ===\n", arith_title, arith_tag);
    
    /* Count extra lines in {title,parm_string}: */
    char *q = title; 
    while ((*q) != 0) { if ((*q) == '\n') { ncap++; } q++; }
    q = parm_string; 
    while ((*q) != 0) { if ((*q) == '\n') { ncap++; } q++; }

    double capFontSize = 10.0; /* Caption font nominal height in pt. */

    /* Compute the plot scale so that the largest dimension is 6 inches: */
    Float wx = xd.hi - xd.lo;
    Float wy = yd.hi - yd.lo;
    Float wm = (wx > wy ? wx : wy);
    double scale = 6.0 * 72.0 / wm;
    
    /* Compute the canvas dimensions (excluding the margins): */
    double hsize = scale*wx; /* Canvas H size. */
    double vsize = scale*wy; /* Canvas V size. */
    double mrg = 4.0; /* In pt. */

    double capHeight = ncap*capFontSize;
    double botMrg = mrg + capHeight + mrg;

    /* Create the EPS figure: */
    bool_t verbose = FALSE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( "out", prefix, arith_tag, -1, NULL, 
        hsize, vsize, mrg, mrg, botMrg, mrg, 
        verbose
      );

    epswr_set_client_window(eps, xd.lo, xd.hi, yd.lo, yd.hi);
    
    epswr_set_text_geometry(eps, FALSE, 0, hsize, -capHeight-mrg, -mrg, 0.0);
    epswr_set_text_font(eps, "Courier", capFontSize);
    epswr_set_fill_color(eps, 0.0,0.0,0.0);
    
    epswr_text(eps, title, FALSE, 0.0, TRUE, FALSE);
    epswr_text(eps, "", FALSE, 0.0, TRUE, FALSE);
    epswr_text(eps, parm_string, FALSE, 0.0, TRUE, FALSE);
    epswr_text(eps, "", FALSE, 0.0, TRUE, FALSE);
    epswr_text(eps, arith_title, FALSE, 0.0, TRUE, FALSE);
    
    return eps;
 }
 
void fn2_zf_plot_stat_caption
  ( epswr_figure_t *eps, 
    char *fmt, 
    int32_t count,
    char *arith
  )
  { char *str = NULL;
    if (arith != NULL)
      { char *str = jsprintf(fmt, count, arith); }
    else
      { char *str = jsprintf(fmt, count); }
    fprintf(stderr, "%s\n", str);
    epswr_set_fill_color(eps, 0.0,0.0,0.0);
    epswr_text(eps, str, FALSE, 0.0, TRUE, FALSE);
    free(str);
  }
 
void fn2_zf_plot_end_figure
  ( epswr_figure_t *eps,
    eval_fp_t *eval_fp,
    Interval xd,
    Interval yd,
    int32_t m, 
    int32_t cols, 
    int32_t rows
  )
  {
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    fn2_zf_flt_plot(eps, eval_fp, xd, yd, m);
    epswr_frame(eps);
  }

char *fn2_zf_plot_format_parms(
    Interval xd,
    Interval yd,
    int32_t n
  )
  { char *s = (char *) malloc(100);
    sprintf(s, 
      "x in [%f _ %f]  y in [%f _ %f]  %d intervals",
      xd.lo, xd.hi, yd.lo, yd.hi, n
    );
    return(s);
  }

