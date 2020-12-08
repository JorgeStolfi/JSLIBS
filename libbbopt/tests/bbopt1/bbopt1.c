/* Univariate nonlinear branch-and-bound optimization with interval estimators. */
/* Last edited on 2009-01-06 04:51:45 by stolfi */

#include <bbopt.h>
#include <bbgoal.h>

#include <fbox.h>
#include <fboxheap.h>
#include <fboxlist.h>

#include <aa.h>
#include <ia.h>
#include <fltgraph.h>
#include <pswr.h>
#include <affirm.h>
#include <jsstring.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* INTERNAL PROTOS */

int main(int argc, char **argv);

PSStream *bb1_new_plot_stream
  ( bool_t epsfmt, 
    char *prefix,
    char *title,
    Interval xr, 
    Interval fr,
    double hsize,
    double vsize, 
    double *ysc
  );
  /* Opens the output plot file and initializes scales, captions, etc.
    If {epsfmt} is true the file is encapsulated Postscript
    ("{outname}.eps") otherwise it is a plain Postscripr
    ("{outname}.ps").
    
    The procedure also stores in {*ysc} a scale factor} to be applied to
    the function values before plotting. The client plot window will
    contain the rectangle {xr × ((*ysc)*fr)}. */

void bb1_print_range_scale
  ( char *which,
    Interval r,
    double psize,
    double scale
  );
  /* Prints the client range {r}, the device window size {psize} (in pt),
    and the client-to-pswr scale {scale}, for axis {which} ("X" or "Y"). */
    
void bb1_print_interval(Interval *xr, Interval *fr, bool_t final);
  /* A reporting function that prints the domain interval {*xr} and 
    its image {*fr} to {stderr}. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  { 
    /* Get function tag from command line: */
    demand(argc == 2, "bad arguments");
    char *ftag = argv[1];
    fprintf(stderr, "ftag = %s\n", ftag);
    
    /* Get function attributes: */
    bbgoal_data_t f = bbgoal_from_tag(ftag);
    demand(f.dim == 1, "Goal function has wrong domain dim");
    affirm(strcmp(f.tag, ftag) == 0, "Goal function has wrong tag");
    
    /* Postscript file name prefix: */
    char *prefix = NULL;
    asprintf(&prefix, "bb1-out-%s-", f.tag);

    /* These variables will be set later: */
    PSStream *ps = NULL;
    double yscale;

    auto void draw_interval(Interval *xr, Interval fr, bool_t final);

    void draw_interval(Interval *xr, Interval fr, bool_t final)
    { bb1_print_interval(&(xr[0]), &fr, final);
        Interval x = xr[0];
        Interval y = ia_scale(fr, yscale, 1.0);
        if (final)
          { pswr_set_fill_color(ps, 0.5,0.5,0.5);
            pswr_rectangle(ps, x.lo, x.hi, y.lo, y.hi, TRUE, TRUE);
          }
        else
          { pswr_rectangle(ps, x.lo, x.hi, y.lo, y.hi, FALSE, TRUE); }
      }

    ia_init();
    
    /* Choose the initial domain {xr} and get its width {xw}: */
    Interval xr = (Interval){0.0, 1.0};
    Float xw = xr.hi - xr.lo;
    Float tol = 0.01 * xw;
    
    /* Get the approximate function value range {fr}: */
    Interval fr = f.eval_ia(&xr);

    /* Ensure that the function value range {fr} includes the Y=0 axis: */
    if (fr.lo > 0.0) { fr.lo = 0.0; }
    if (fr.hi < 0.0) { fr.hi = 0.0; }

    /* Get the true global minimum point {sr} and value {mr}: */
    Interval sr;
    f.true_opt(&xr, &sr);
    Interval mr = f.eval_ia(&sr);

    /* Postscript figure size: */
    double hsize = 432.0;
    double vsize = 288.0;
    
    /* Generate EPS and PS output: */
    int epsfmt;
    for (epsfmt = 0; epsfmt < 2; epsfmt++)
      { 
        ps = bb1_new_plot_stream
          ( epsfmt, prefix, f.descr, 
            xr, fr, hsize, vsize, &yscale
          );
        bb1_print_range_scale("X", xr, hsize, 1.00);
        bb1_print_range_scale("F", fr, vsize, yscale);
        
        /* Draw axes and plot of {F}: */
        auto Float fplot(Float x);
        Float fplot(Float x) { return yscale * f.eval_fp(&x); }
        Interval xp = xr;
        Interval yp = ia_scale(fr, yscale, 1.0);
        pswr_set_pen(ps, 0.0,0.2,0.5, 0.25, 0.0,0.0);
        fltgraph_plot(ps, fplot, xp, yp, 100);
        
        /* Plot the true minimum: */
        Float xsol = (sr.lo + sr.hi)/2;
        Float ysol = fplot(xsol);
        pswr_set_pen(ps, 0.5,0.0,0.0, 0.25, 0.0,0.0);
        pswr_dot(ps, xsol, ysol, 2.0, FALSE, TRUE);
        
        /* Optimize and plot node boxes: */
        pswr_set_pen(ps, 0.0,0.5,0.0, 0.10, 0.0,0.0);
        Interval gr;
        FBoxList R = bb_optimize(1, f.eval_ia, &xr, &tol, &gr, draw_interval); 
        affirm(R != NULL, "bb_optimize returned no boxes");
        fprintf(stderr, "minimum found = [%8.4f _ %8.4f]\n", gr.lo, gr.hi);
        
        /* Print the true minimum: */
        fprintf(stderr, "true minimum:\n");
        bb1_print_interval(&sr, &mr, FALSE);

        /* Draw frame: */
        pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0.0,0.0);
        pswr_frame(ps);
        
        /* Terminate plot: */
        pswr_close_stream(ps);
      }
    return 0;
  }

PSStream *bb1_new_plot_stream
  ( bool_t epsfmt, 
    char *prefix,
    char *title,
    Interval xr, 
    Interval fr,
    double hsize,
    double vsize, 
    double *ysc
  )
  { /* Compute widths of {xr} and {yr}: */
    Interval yr = fr;
    double xw = xr.hi - xr.lo;
    double yw = yr.hi - yr.lo;
    
    /* Add a safety margin around {xr,yr}: */
    double xmrg = 0.03*xw;
    double ymrg = 0.03*yw;
    xr.lo -= xmrg; xr.hi += xmrg;
    yr.lo -= ymrg; yr.hi += ymrg;

    /* Recompute {xw,yw} for expanded intervals: */
    xw = xr.hi - xr.lo;
    yw = yr.hi - yr.lo;
    
    /* Compute {yscale} to use the whole Y range: */
    double yscale = (xw*vsize)/(yw*hsize);

    /* Create the stream: */
    double mrg = 4.0; /* Figure margin, if EPS format. */
    PSStream *ps = pswr_new_stream(prefix, NULL, epsfmt, "doc", "letter", FALSE, hsize + 2*mrg, vsize + 2*mrg);
    pswr_new_picture(ps, xr.lo, xr.hi, yscale*yr.lo, yscale*yr.hi);
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0,0.0);
    if (! epsfmt) { pswr_add_caption(ps, title, 0.0); }
    *ysc = yscale;
    return ps;
  }


void bb1_print_range_scale
  ( char *which,
    Interval r,
    double psize,
    double scale
  )
  {
    fprintf(stderr, "%s range  = ", which);
    ia_print(stderr, r); 
    fprintf(stderr, "  size  = %6.1f", psize);
    fprintf(stderr, "  scale  = %9.4f", scale);
    fprintf(stderr, "\n");
  }

void bb1_print_interval(Interval *xr, Interval *fr, bool_t final)
  { ia_print(stderr, *xr); 
    fprintf(stderr, " × ");
    ia_print(stderr, *fr); 
    if (final) { fprintf(stderr, " *"); }
    fprintf(stderr, "\n");
  }
