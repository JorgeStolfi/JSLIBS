/* See {fn1_zf.h} */
/* Last edited on 2023-02-18 09:32:49 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <malloc.h>

#include <fltgraph.h>
#include <zf.h>
#include <aa_trapez.h>
#include <aa.h>
#include <ia.h>
#include <ia_butfly.h>
#include <flt.h>
#include <affirm.h>
#include <jsstring.h>
#include <jsfile.h>
#include <epswr.h>

#include <fn1_zf.h>
#include <fn1_functions.h>

#define DEBUG 1

/*** CLOSURE TYPES ***/

typedef struct fn1_zf_problem_t
  { char *title;
    Interval xd;
    Interval yd;
    eval_fp_t *eval_fp;  
    eval_ia_t *eval_ia;
    diff_ia_t *diff_ia;
    eval_aa_t *eval_aa;    
    double epsilon;
    double delta;
    int32_t m;
  } fn1_zf_problem_t;
  
typedef ia_butfly_t eval_butt_t (Interval xv, fn1_zf_problem_t *p);

typedef struct fn1_zf_closure_t
  { FILE *pplot;           /* Progress plot. */
    FILE *splot;           /* Positive/zero/negative plot. */
    eval_butt_t *eval;     /* Function evaluator, returns butterfly enclosure. */
    fn1_zf_problem_t *p; /* Problem description. */
    int32_t nevals;            /* Count of evals performed. */
    int32_t nroots;            /* Number of roots found. */
  } fn1_zf_closure_t;
  

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

void fn1_zf_do_find_and_plot_zeros
  ( char *prefix,
    char *arith_tag,
    char *arith_title,
    fn1_zf_problem_t *p,
    eval_butt_t eval
  );
  /* The {eval} function should compute a butterfly containing the 
    graph of the function over the interval {xv}, without any further
    processing. */

epswr_figure_t *fn1_zf_create_figure
  ( char *prefix, 
    char *arith_tag,
    char *suffix,
    char *arith_title,
    char *subtitle,
    fn1_zf_problem_t *p
  );
  /* Opens an Encapsulated Postscript (EPS) file named 
    "out/zf1-{arith_tag}-{suffix}.eps" and writes the captions
    {arith_title} and {subtitle}. */
  
void fn1_zf_finish_figure
  ( epswr_figure_t *eps,
    int32_t nevals,      /* Number of function evaluations (for caption). */
    int32_t nroots,      /* Number of roots found (for caption). */
    fn1_zf_problem_t *p
  );
  /* Writes the postamble of a plot file {psfile}, 
    including the function's plot, axes, frame, captions, etc.. */
  
char *fn1_zf_format_parms
  ( Interval xd,
    Interval yd,
    double epsilon,
    double delta
  );
  /* Formats the arguments into a string, for figure captions. */
  
void fn1_zf_report_stats
  ( FILE *pplot, FILE *splot,
    int32_t epsformat, 
    char *arith,    /* "IA" or "AA" */
    int32_t nevals, 
    int32_t nroots
  );
  /* Prints evaluation report, and adds caption lines if not {epsformat} */
  
ia_butfly_t fn1_zf_ia_eval (Interval xv, fn1_zf_problem_t *p);
  /* Computes a rectangular enclosing box for the graph of {F} over the
    interval {xv}. The result has {xmd=xv.lo} and {yxlo=yxmd=yxhi=yv}
    where {yv = p->eval_ia(xv)}. */

ia_butfly_t fn1_zf_ar_eval (Interval xv, fn1_zf_problem_t *p);
  /* Computes a  rectangular enclosing box for the graph of {F} over the
    interval {xv}, from the affine form {p->eval_aa(from_interval(xv))}. */

ia_butfly_t fn1_zf_id_eval (Interval xv, fn1_zf_problem_t *p);
  /* Computes a butterfly-shaped enclosing region for the graph of
    {F}, over the X range {xv}. Uses the interval-slope model with
    central interval {yxmd = p->eval_ia(xm)}, where {xm} is a
    singleton interval at the center {xmd} of {xv}. The slope range is
    obtained by {dv = p->diff_ia(xv)}. */

ia_butfly_t fn1_zf_aa_eval (Interval xv, fn1_zf_problem_t *p);
  /* Computes a parallelogram-shaped enclosure for the graph of {F}, 
    from the affine form {p->eval_aa(from_interval(xv))}. */

/*** IMPLEMENTATIONS ***/

void fn1_zf_find_and_plot_zeros
  ( char *prefix,
    eval_fp_t eval_fp,
    eval_ia_t eval_ia,
    diff_ia_t diff_ia,
    eval_aa_t eval_aa,
    char *title,
    Interval xd,
    Interval yd,
    double epsilon,
    double delta,
    int32_t m
  )
  { fn1_zf_problem_t p;
    
    p.eval_fp = eval_fp;
    p.eval_ia = eval_ia;
    p.diff_ia = diff_ia;
    p.eval_aa = eval_aa;
    p.title = title;
    p.xd = xd;
    p.yd = yd;
    p.epsilon = epsilon;
    p.delta = delta;
    p.m = m;

    /*** Box enclosures (Interval Arith.) ***/
    fn1_zf_do_find_and_plot_zeros
      ( prefix, "ia", "Box enclosures (Interval Arith.)",
        &p, fn1_zf_ia_eval 
      );

    /*** Box enclosures (Affine Arith.) ***/
    fn1_zf_do_find_and_plot_zeros
      ( prefix, "ar", "Box enclosures (Affine Arith.)",
        &p, fn1_zf_ar_eval 
      );

    /*** Interval-slope enclosures (Interval Arith.) ***/
    fn1_zf_do_find_and_plot_zeros
      ( prefix, "id", "Interval-slope enclosures (Interval Arith.)",
        &p, fn1_zf_id_eval 
      );

    /*** Slab enclosures (Affine Arith.) ***/
    fn1_zf_do_find_and_plot_zeros
      ( prefix, "aa", "Slab enclosures (Affine Arith.)",
        &p, fn1_zf_aa_eval 
      );
  }
  
void fn1_zf_do_find_and_plot_zeros
  ( char *prefix,
    char *arith_tag,
    char *arith_title,
    fn1_zf_problem_t *p,
    eval_butt_t eval
  )
  { fprintf(stderr, "\nenter fn1_zf_do_find_and_plot_zeros...\n");

    int32_t nevals = 0;
    int32_t nroots = 0;
    
    auto void fn1_zf_eval(Interval *x, Interval *y, ia_butfly_t *bt);
      /* Called by {zf_enum_zeros} to bound the graph of {F} and check for
        roots. This procedure calls {closure->eval(xv, closure->p)} and
        plots the resulting butterfly. The {data} argument should be
        actually an {fn1_zf_closure_t}. */

    auto bool_t fn1_zf_report(Interval *x, Interval *y, zf_kind_t kind);
      /* Called by {zf_enum_zeros} to report a new interval.
        This procedure plots the given root and prints it to {stderr}. */
    /* Progress plot: */
    epswr_figure_t *pplot = fn1_zf_create_figure(prefix, arith_tag, "p", arith_title, "Search progress", p);
    
    /* Sign classification plot: */
    epswr_figure_t *splot = fn1_zf_create_figure(prefix, arith_tag, "s", arith_title, "Function signs", p);
    
    epswr_set_pen(pplot, 0.0, 0.5, 0.0, 0.25, 0.0, 0.0);
    epswr_set_pen(splot, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0);
    zf_enum_zeros
      ( fn1_zf_eval, 
        fn1_zf_report,
        p->xd,
        p->epsilon,
        p->delta
      );
    
    fn1_zf_finish_figure(pplot, nevals, nroots, p);
    fn1_zf_finish_figure(splot, nevals, nroots, p);
    fprintf(stderr, "exit fn1_zf_do_find_and_plot_zeros.\n");
    
    epswr_end_figure(pplot);
    epswr_end_figure(splot);

    return;
    
    /* Internal procedures: */

    void fn1_zf_eval (Interval *x, Interval *y, ia_butfly_t *bt)
      { *bt = eval(*x, p);
        Interval ybt0 = ia_join(bt->tp[0].yxlo, bt->tp[0].yxhi);
        Interval ybt1 = ia_join(bt->tp[1].yxlo, bt->tp[1].yxhi);
        *y = ia_join(ybt0, ybt1);
        nevals++;
        ia_butfly_draw(pplot, &(p->yd), bt);
      }

    bool_t fn1_zf_report(Interval *x, Interval *y, zf_kind_t kind)
      { nroots++;
        fprintf(stderr, "  apparent root in "); 
        ia_print(stderr, *x); 
        fprintf(stderr, "\n");

        { /* Plot result: */
          float cr, cg, cb; /* Polygon color */
          switch (kind)
            { case zf_kind_undefined:
                cr = 0.75f; cg = 0.75f; cb = 0.75f; break;
              case zf_kind_positive:
                cr = 1.00f; cg = 0.25f; cb = 0.25f; break;
              case zf_kind_negative:
                cr = 0.50f; cg = 0.60f; cb = 1.00f; break;
              case zf_kind_root:
                cr = 0.00f; cg = 0.85f; cb = 0.00f; break;
              case zf_kind_mixed:
                cr = 0.00f; cg = 0.00f; cb = 0.00f; break;
            }

          ia_trapez_t tr = (ia_trapez_t){ *x, *y, *y };
          ia_trapez_fill(splot, &(p->yd), &tr, cr,cg,cb);
        }
        return 0;
      }
  }
  
epswr_figure_t *fn1_zf_create_figure
  ( char *prefix, 
    char *arith_tag,
    char *suffix,
    char *arith_title,
    char *subtitle,
    fn1_zf_problem_t *p
  )
  { 
    fprintf(stderr, "starting figure\n");
    fprintf(stderr, "title = %s\n", p->title);
    fprintf(stderr, "arith = %s (%s)\n", arith_title, arith_tag);
    fprintf(stderr, "subtitle = %s\n", subtitle);
    
    /* Count lines in {p->title}: */
    int32_t nl = 1;
    char *q = p->title; 
    while ((*q) != 0) { if ((*q) == '\n') { nl++; } q++; }

    double capFontSize = 10.0; /* Caption font nominal height in pt. */

   /* Pack the parameters as a string: */
    char *parm_string = fn1_zf_format_parms(p->xd, p->yd, p->epsilon, p->delta);
    
    /* Compute the plot scale so that the largest dimension is 6 inches: */
    Float wx = p->xd.hi - p->xd.lo;
    Float wy = p->yd.hi - p->yd.lo;
    
    /* Compute the canvas dimensions (excluding the margins): */
    Float wm = (wx > wy ? wx : wy);
    double scale = 6.0 * 72.0 / wm;
    double hsize = scale*wx; /* Canvas H size. */
    double vsize = scale*wy; /* Canvas V size. */

    /* Define margins and space for caption at bottom: */
    double mrg = 4.0; /* Margin around plot area in pt. */
    double capHeight = (6 + nl)*capFontSize;
    double botMrg = mrg + capHeight + mrg;
    
    bool_t verbose = TRUE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( "out", prefix, arith_tag, -1, suffix, 
        hsize, vsize, mrg, mrg, botMrg, mrg, 
        verbose
      );
   
    /* Set the client window size: */
    epswr_set_client_window(eps, p->xd.lo, p->xd.hi, p->yd.lo, p->yd.hi);
    
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.25, 0.0,0.0);
    
    epswr_set_text_geometry(eps, FALSE, 0, hsize, -capHeight-mrg, -mrg, 0.0);
    epswr_set_text_font(eps, "Courier", capFontSize);
    epswr_set_fill_color(eps, 0.0,0.0,0.0);
    
    epswr_text(eps, p->title, FALSE, 0.0, TRUE, FALSE);
    epswr_text(eps, txtcat("Arithmetic: ", arith_title), FALSE, 0.0, TRUE, FALSE);
    epswr_text(eps, subtitle, FALSE, 0.0, TRUE, FALSE);
    epswr_text(eps, "", FALSE, 0.0, TRUE, FALSE);
    epswr_text(eps, parm_string, FALSE, 0.0, TRUE, FALSE);

    return eps;
  }

void fn1_zf_finish_figure
  ( epswr_figure_t *eps,
    int32_t nevals,
    int32_t nroots,
    fn1_zf_problem_t *p
  )
  { char *sevals = NULL;
    char *sevals = jsprintf(
      "%d function evaluations, %d roots found", 
      nevals, nroots
    );
    fprintf(stderr, "%s\n", sevals);
    epswr_set_fill_color(eps, 0.0,0.0,0.0);
    epswr_text(eps, sevals, FALSE, 0.0, TRUE, FALSE);
    free(sevals);
    
    epswr_set_pen(eps, 0.0, 0.0, 0.0, 0.10, 0.0, 0.0);
    fltgraph_draw_axes(eps, p->xd, p->yd);
     
    epswr_set_pen(eps, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0);
    epswr_frame(eps);
   
    epswr_set_pen(eps, 0.5, 0.0, 0.0, 0.35, 0.0, 0.0);
    fltgraph_plot(eps, p->eval_fp, p->xd, p->yd, p->m);
  }

ia_butfly_t fn1_zf_ia_eval (Interval xv, fn1_zf_problem_t *p)
  { Interval fv = p->eval_ia(xv);
#   if DEBUG
      fprintf(stderr, "  f =  "); ia_print(stderr, fv); fprintf(stderr, "\n");
#   endif
    return ia_butfly_from_box(&xv, &fv);
  }

ia_butfly_t fn1_zf_id_eval (Interval xv, fn1_zf_problem_t *p)
  { /* Compute Y range at central point of {xv}: */
    Float xmd = ia_mid(xv);
    Interval xm = (Interval){xmd, xmd};
    Interval yxmd = p->eval_ia(xm);
    /* Compute slope range over whole range {xv}: */
    Interval dfv = p->diff_ia(xv);
#   if DEBUG
      fprintf(stderr, "  f =  "); ia_print(stderr, yxmd); fprintf(stderr, "\n");
      fprintf(stderr, "  df = "); ia_print(stderr, dfv); fprintf(stderr, "\n");
#   endif
    return ia_butfly_from_ia_diff(&xv, xmd, &yxmd, &dfv);
  }

ia_butfly_t fn1_zf_aa_eval (Interval xv, fn1_zf_problem_t *p)
  { MemP frame = aa_top();
    AAP xf = aa_from_interval(xv);
    AAP yf = p->eval_aa(xf);
#   if DEBUG
      fprintf(stderr, "  f =  "); aa_print(stderr, yf); fprintf(stderr, "\n");
#   endif

    ia_trapez_t tr = aa_trapez_from_pair(&xv, xf, yf);
    aa_flush(frame);
    return ia_butfly_from_trapez(&tr);
  }

ia_butfly_t fn1_zf_ar_eval (Interval xv, fn1_zf_problem_t *p)
  { MemP frame = aa_top();
    AAP xf = aa_from_interval(xv);
    AAP yf = p->eval_aa(xf);
#   if DEBUG
      fprintf(stderr, "  f =  "); aa_print(stderr, yf); fprintf(stderr, "\n");
#   endif
    Interval yv = aa_range(yf);
    ia_butfly_t bt = ia_butfly_from_box(&xv, &yv);
    aa_flush(frame);
    return bt;
  }

char *fn1_zf_format_parms
  ( Interval xd,
    Interval yd,
    double epsilon,
    double delta
  )
  { char *s = (char *) malloc(100);
    sprintf(s, 
      "x in [%f _ %f]  y in [%f _ %f]\nepsilon = %.4e delta = %.4e",
      xd.lo, xd.hi, yd.lo, yd.hi, epsilon, delta
    );
    return s;
  }

