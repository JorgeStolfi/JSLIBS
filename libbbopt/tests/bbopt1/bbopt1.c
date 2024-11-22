/* Univariate nonlinear branch-and-bound optimization with interval estimators. */
/* Last edited on 2023-02-20 10:03:36 by stolfi */

#define _GNU_SOURCE
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <jsstring.h>
#include <affirm.h>
#include <epswr.h>
#include <fltgraph.h>
#include <ia.h>
#include <aa.h>

#include <fboxlist.h>
#include <fboxheap.h>
#include <fbox.h>

#include <bbgoal.h>
#include <bbopt.h>

/* INTERNAL PROTOS */

int32_t main(int32_t argc, char **argv);

epswr_figure_t *bb1_new_figure
  ( char *outname,
    char *title,
    Interval xr, 
    Interval fr,
    double *yscale_P
  );
  /* Opens the output Encapsulated Postscript file, called "{outname}.eps".
    Writes the {title} as caption. under the plot.
    
    The procedure also stores in {*yscale_P} a scale factor {yscale} to be applied to
    the function values before plotting. The client plot window will
    contain the rectangle {xr × yscale*fr}. */

void bb1_print_range_scale
  ( epswr_figure_t *eps,
    char *which,
    Interval r,
    double scale
  );
  /* Prints the actual plot range {r} and scale factor
    {scale} applied before plotting, for axis {which} ("X" or "Y"). */
    
void bb1_print_interval(Interval *xr, Interval *fr, bool_t final);
  /* A reporting function that prints the domain interval {*xr} and 
    its image {*fr} to {stderr}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { 
    /* Get function tag from command line: */
    demand(argc == 2, "bad arguments");
    char *ftag = argv[1];
    fprintf(stderr, "ftag = %s\n", ftag);
    
    /* Get function attributes: */
    bbgoal_data_t f = bbgoal_from_tag(ftag);
    demand(f.dim == 1, "Goal function has wrong domain dim");
    affirm(strcmp(f.tag, ftag) == 0, "Goal function has wrong tag");

    ia_init();
    
    /* Choose the initial domain {xr} and get its width {xw}: */
    Interval xr = (Interval){0.0, 1.0};
    Float xw = xr.hi - xr.lo;
    Float tol = (Float)(0.01 * xw);
    
    /* Get the approximate function value range {fr}: */
    Interval fr = f.eval_ia(&xr);

    /* Ensure that the function value range {fr} includes the Y=0 axis: */
    if (fr.lo > 0.0) { fr.lo = 0.0; }
    if (fr.hi < 0.0) { fr.hi = 0.0; }
    
    /* Postscript file: */
    char *outname = jsprintf("bb1_%s", f.tag);
    double yscale;
    epswr_figure_t *eps = bb1_new_figure(outname, f.descr, xr, fr, &yscale);

    auto void draw_interval(Interval *xr, Interval fr, bool_t final);

    void draw_interval(Interval *xr, Interval fr, bool_t final)
      { bb1_print_interval(&(xr[0]), &fr, final);
        Interval x = xr[0];
        Interval y = ia_scale(fr, (Float)yscale, 1.0);
        if (final)
          { epswr_set_fill_color(eps, 0.5,0.5,0.5);
            epswr_rectangle(eps, x.lo, x.hi, y.lo, y.hi, TRUE, TRUE);
          }
        else
          { epswr_rectangle(eps, x.lo, x.hi, y.lo, y.hi, FALSE, TRUE); }
      }

    /* Get the true global minimum point {sr} and value {mr}: */
    Interval sr;
    f.true_opt(&xr, &sr);
    Interval mr = f.eval_ia(&sr);
   
    /* Generate EPS output: */
    bb1_print_range_scale(eps, "X", xr, 1.00);
    bb1_print_range_scale(eps, "F", fr, yscale);

    /* Draw axes and plot of {F}: */
    auto Float fplot(Float x);
    
    Float fplot(Float x) { return (Float)(yscale * f.eval_fp(&x)); }
    Interval xp = xr;
    Interval yp = ia_scale(fr, (Float)yscale, 1.0);
    epswr_set_pen(eps, 0.0,0.2,0.5, 0.25, 0.0,0.0);
    fltgraph_plot(eps, fplot, xp, yp, 100);

    /* Plot the true minimum: */
    Float xsol = (sr.lo + sr.hi)/2;
    Float ysol = fplot(xsol);
    epswr_set_pen(eps, 0.5,0.0,0.0, 0.25, 0.0,0.0);
    epswr_dot(eps, xsol, ysol, 2.0, FALSE, TRUE);

    /* Optimize and plot node boxes: */
    epswr_set_pen(eps, 0.0,0.5,0.0, 0.10, 0.0,0.0);
    Interval gr;
    FBoxList R = bb_optimize(1, f.eval_ia, &xr, &tol, &gr, draw_interval); 
    affirm(R != NULL, "bb_optimize returned no boxes");
    fprintf(stderr, "minimum found = [%8.4f _ %8.4f]\n", gr.lo, gr.hi);

    /* Print the true minimum: */
    fprintf(stderr, "true minimum:\n");
    bb1_print_interval(&sr, &mr, FALSE);

    /* Draw frame: */
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.25, 0.0,0.0);
    epswr_frame(eps);

    /* Terminate plot: */
    epswr_end_figure(eps);

    return 0;
  }

epswr_figure_t *bb1_new_figure
  ( char *outname,
    char *title,
    Interval xr, 
    Interval fr,
    double *yscale_P
  )
  { /* Compute plot ranges from {xr,fr} with some skosh: */
    double xw = xr.hi - xr.lo;
    double xmrg = 0.03*xw;
    double xmin = xr.lo - xmrg, xmax = xr.hi + xmrg;

    double fw = fr.hi - fr.lo;
    double yscale = xw/fw;
    
    double yw = xw;
    double ymrg = 0.03*yw;
    double ymin = yscale*fr.lo - ymrg, ymax = yscale*fr.hi + ymrg;
   
    double maxSize = 150*epswr_pt_per_mm;
    int32_t capLines = 1;
    double fontHeight = 10.0;
    bool_t eps_verbose = FALSE;
    epswr_figure_t *eps = epswr_new_captioned_figure
      ( "out", outname, NULL, -1, NULL,
        xmin,xmax, ymin,ymax, maxSize, maxSize,
        capLines, fontHeight, eps_verbose
      );
    epswr_text(eps, title, FALSE, 0.5, TRUE, FALSE);
    (*yscale_P) = yscale;
    return eps;
  }

void bb1_print_range_scale
  ( epswr_figure_t *eps,
    char *which,
    Interval r,
    double scale
  )
  {
    fprintf(stderr, "%s range  = ", which);
    ia_print(stderr, r); 
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
