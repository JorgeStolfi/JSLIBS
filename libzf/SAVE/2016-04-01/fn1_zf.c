/* See {fn1_zf.h} */
/* Last edited on 2011-05-29 10:42:04 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <malloc.h>

#include <fn1_zf.h>
#include <fn1_functions.h>

#include <fltgraph.h>
#include <zf.h>
#include <aa_trapez.h>
#include <aa.h>
#include <ia.h>
#include <ia_butfly.h>
#include <flt.h>
#include <affirm.h>
#include <jsstring.h>
#include <pswr.h>

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
    int m;
  } fn1_zf_problem_t;
  
typedef ia_butfly_t eval_butt_t (Interval xv, fn1_zf_problem_t *p);

typedef struct fn1_zf_closure_t
  { FILE *pplot;           /* Progress plot. */
    FILE *splot;           /* Positive/zero/negative plot. */
    eval_butt_t *eval;     /* Function evaluator, returns butterfly enclosure. */
    fn1_zf_problem_t *p; /* Problem description. */
    int nevals;            /* Count of evals performed. */
    int nroots;            /* Number of roots found. */
  } fn1_zf_closure_t;
  

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

void fn1_zf_do_find_and_plot_zeros
  ( PSStream *pplot,
    PSStream *splot,
    char *arith_tag,
    char *arith_title,
    fn1_zf_problem_t *p,
    eval_butt_t eval
  );
  /* The {eval} function should compute a butterfly containing the 
    graph of the function over the interval {xv}, without any further
    processing. */

void fn1_zf_open_ps_files
  ( char *prefix,
    bool_t epsformat,
    PSStream **pplotp,   /* Zero-finding progress plot. */
    PSStream **splotp    /* Function signs plot. */
  );
  /* Opens Postscript files for the progress plots ({*pplotp}) and
    for the sign (positive/negative/zero) classification of intervals 
    ({*splotp}).  If {epsformat} is TRUE, each figure will be written
    as "{prefix}-{suffix}-{tag}.eps", where {tag} is specified
    to {fn1_zf_begin_figure} and the {suffix} is "p" or "s",
    respectively.  If {epsformat} is FALSE, each plot will be written
    as a document called "{prefix}-{suffix}.ps", and the {tag} will
    be used as page number.  */
  
void fn1_zf_start_figure
  ( PSStream *ps,
    char *arith_tag,
    char *arith_title,
    char *subtitle,
    fn1_zf_problem_t *p
  );
  /* Writes the preamble of a plot file {psfile}, sets scales, etc. */
  
void fn1_zf_finish_figure
  ( PSStream *ps,
    int nevals,      /* Number of function evaluations (for caption). */
    int nroots,      /* Number of roots found (for caption). */
    fn1_zf_problem_t *p
  );
  /* Writes the postamble of a plot file {psfile}, 
    including the function's plot, axes, frame, captions, etc.. */
  
char *fn1_zf_format_parms
  ( Interval xd,
    Interval yd,
    Float epsilon,
    Float delta
  );
  /* Formats the arguments into a string, for figure captions. */
  
void fn1_zf_report_stats
  ( FILE *pplot, FILE *splot,
    int epsformat, 
    char *arith,    /* "IA" or "AA" */
    int nevals, 
    int nroots
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
  ( char *fileprefix,
    int epsformat,
    eval_fp_t eval_fp,
    eval_ia_t eval_ia,
    diff_ia_t diff_ia,
    eval_aa_t eval_aa,
    char *title,
    Interval xd,
    Interval yd,
    Float epsilon,
    Float delta,
    int m
  )
  { PSStream *pplot;
    PSStream *splot;
    fn1_zf_problem_t p;
    
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
    
    fn1_zf_open_ps_files(fileprefix, epsformat, &pplot, &splot);

    /*** Box enclosures (Interval Arith.) ***/
    fn1_zf_do_find_and_plot_zeros
      ( pplot, splot,
        "ia", "Box enclosures (Interval Arith.)",
        &p, fn1_zf_ia_eval 
      );

    /*** Box enclosures (Affine Arith.) ***/
    fn1_zf_do_find_and_plot_zeros
      ( pplot, splot,
        "ar", "Box enclosures (Affine Arith.)",
        &p, fn1_zf_ar_eval 
      );

    /*** Interval-slope enclosures (Interval Arith.) ***/
    fn1_zf_do_find_and_plot_zeros
      ( pplot, splot,
        "id", "Interval-slope enclosures (Interval Arith.)",
        &p, fn1_zf_id_eval 
      );

    /*** Slab enclosures (Affine Arith.) ***/
    fn1_zf_do_find_and_plot_zeros
      ( pplot, splot,
        "aa", "Slab enclosures (Affine Arith.)",
        &p, fn1_zf_aa_eval 
      );
    
    pswr_close_stream(pplot);
    pswr_close_stream(splot);
  }
  
void fn1_zf_open_ps_files
  ( char *prefix,
    bool_t epsformat,
    PSStream **pplotp,
    PSStream **splotp
  )
  { 
    /* The actual figure size (if {epsformat=TRUE}) will be set later: */
    (*pplotp) = pswr_new_stream(prefix, NULL, epsformat, "p", "letter", 300, 300);
    (*splotp) = pswr_new_stream(prefix, NULL, epsformat, "s", "letter", 300, 300);
  }
  
void fn1_zf_do_find_and_plot_zeros
  ( PSStream *pplot,
    PSStream *splot,
    char *arith_tag,
    char *arith_title,
    fn1_zf_problem_t *p,
    eval_butt_t eval
  )
  { fprintf(stderr, "\nenter fn1_zf_do_find_and_plot_zeros...\n");

    int nevals = 0;
    int nroots = 0;
    
    auto void fn1_zf_eval(Interval *x, Interval *y, ia_butfly_t *bt);
      /* Called by {zf_enum_zeros} to bound the graph of {F} and check for
        roots. This procedure calls {closure->eval(xv, closure->p)} and
        plots the resulting butterfly. The {data} argument should be
        actually an {fn1_zf_closure_t}. */

    auto bool_t fn1_zf_report(Interval *x, Interval *y, zf_kind_t kind);
      /* Called by {zf_enum_zeros} to report a new interval.
        This procedure plots the given root and prints it to {stderr}. */
    
    fn1_zf_start_figure(pplot, arith_tag, arith_title, "Search progress",  p);
    fn1_zf_start_figure(splot, arith_tag, arith_title, "Function signs",   p);

    pswr_set_pen(pplot, 0.0, 0.5, 0.0, 0.25, 0.0, 0.0);
    pswr_set_pen(splot, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0);
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
                cr = 0.75; cg = 0.75; cb = 0.75; break;
              case zf_kind_positive:
                cr = 1.00; cg = 0.25; cb = 0.25; break;
              case zf_kind_negative:
                cr = 0.50; cg = 0.60; cb = 1.00; break;
              case zf_kind_root:
                cr = 0.00; cg = 0.85; cb = 0.00; break;
              case zf_kind_mixed:
                cr = 0.00; cg = 0.00; cb = 0.00; break;
            }

          ia_trapez_t tr = (ia_trapez_t){ *x, *y, *y };
          ia_trapez_fill(splot, &(p->yd), &tr, cr,cg,cb);
        }
        return 0;
      }
  }
  
void fn1_zf_start_figure
  ( PSStream *ps,
    char *arith_tag,
    char *arith_title,
    char *subtitle,
    fn1_zf_problem_t *p
  )
  { 
    if (ps == NULL)
      { fatalerror ("fn1_zf_do_find_and_plot_zeros: null PostScript stream"); }

    /* Start a new page or new EPS figure: */
    pswr_new_canvas(ps, arith_tag);
    
    /* Pack the parameters as a string: */
    char *parmstr = fn1_zf_format_parms(p->xd, p->yd, p->epsilon, p->delta);
    
    /* Compute the plot scale so that the largest dimension is 6 inches: */
    Float wx = p->xd.hi - p->xd.lo;
    Float wy = p->yd.hi - p->yd.lo;
    Float wm = (wx > wy ? wx : wy);
    double scale = 6.0 * 72.0 / wm;
    
    /* Compute the canvas dimensions (excluding the margins): */
    double hsize = scale*wx; /* Canvas H size. */
    double vsize = scale*wy; /* Canvas V size. */
    double mrg = 2.0; /* In pt. */

    /* Set the canvas size (effective only for EPS figures): */
    pswr_set_canvas_size(ps, hsize+2*mrg, vsize + 2*mrg);
    
    /* Set the client window size: */
    pswr_new_picture(ps, p->xd.lo, p->xd.hi, p->yd.lo, p->yd.hi);
    
    pswr_set_pen(ps, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0);
    if (! pswr_is_eps(ps))
      { pswr_add_caption(ps, p->title, 0.0);
        pswr_add_caption(ps, txtcat("Arithmetic: ", arith_title), 0.0);
        pswr_add_caption(ps, subtitle, 0.0);
        pswr_add_caption(ps, "", 0.0);
        pswr_add_caption(ps, parmstr, 0.0);
      }

    pswr_comment(ps, arith_title);
  }

void fn1_zf_finish_figure
  ( PSStream *ps,
    int nevals,
    int nroots,
    fn1_zf_problem_t *p
  )
  { pswr_set_pen(ps, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0);

    { char *sevals = NULL;
      asprintf(&sevals, 
        "%d function evaluations, %d roots found", 
        nevals, nroots
      );
      fprintf(stderr, "%s\n", sevals);
      if (!  pswr_is_eps(ps)) { pswr_add_caption(ps, sevals, 0.0); }
      free(sevals);
    }
    
    pswr_set_pen(ps, 0.0, 0.0, 0.0, 0.10, 0.0, 0.0);
    fltgraph_draw_axes(ps, p->xd, p->yd);
    
    pswr_set_pen(ps, 0.5, 0.0, 0.0, 0.35, 0.0, 0.0);
    fltgraph_plot(ps, p->eval_fp, p->xd, p->yd, p->m);
    
    pswr_set_pen(ps, 0.0, 0.0, 0.0, 0.25, 0.0, 0.0);
    pswr_frame(ps);
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
    Float epsilon,
    Float delta
  )
  { char *s = (char *) malloc(100);
    sprintf(s, 
      "x in [%f _ %f]  y in [%f _ %f]\nepsilon = %.4e delta = %.4e",
      xd.lo, xd.hi, yd.lo, yd.hi, epsilon, delta
    );
    return s;
  }

