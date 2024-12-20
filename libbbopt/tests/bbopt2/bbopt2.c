/* Bivariate nonlinear branch-and-bound optimization with interval estimators. */
/* Last edited on 2024-12-05 10:22:09 by stolfi */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <jsstring.h>
#include <affirm.h>
#include <epswr.h>
#include <ia.h>
#include <aa.h>

#include <fboxlist.h>
#include <fboxheap.h>
#include <fbox.h>
#include <bbgoal.h>
#include <bbopt.h>

#define D 2

/* INTERNAL PROTOS */

int main(int argc, char **argv);

epswr_figure_t *bb2_new_figure
  ( char *outname,
    char *title,
    Interval xr, 
    Interval yr
  );
  /* Opens the output Encapsulated Postscript file, called "{outname}.eps".
    Writes the {title} as caption. under the plot. The client plot window will
    contain the rectangle {xr × yr}. */
    
void bb2_print_range_scale(char *which, Interval r);
  /* Prints the client range {r} for axis {which} ("X" or "Y"). */
    
void bb2_print_box(Interval *xr, Interval *fr, bool_t final);
  /* A reporting function that prints the domain box {xr[0..D-1]} and 
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
    demand(f.dim == 2, "Goal function has wrong domain dim");
    affirm(strcmp(f.tag, ftag) == 0, "Goal function has wrong tag");
    
    /* File name: */
    char *outname = jsprintf("bb2_%s", f.tag);

    /* These variables will be set later: */
    epswr_figure_t *eps = NULL;

    auto void draw_box(Interval *xr, Interval fr, bool_t final);

    void draw_box(Interval *xr, Interval fr, bool_t final)
      { bb2_print_box(xr, &fr, final);
        Interval x = xr[0];
        Interval y = xr[1];
        if (final)
          { epswr_set_fill_color(eps, 0.5,0.5,0.5);;
            epswr_rectangle(eps, x.lo, x.hi, y.lo, y.hi, TRUE, TRUE);
            
          }
        else
          { epswr_rectangle(eps, x.lo, x.hi, y.lo, y.hi, FALSE, TRUE); }
      }

    ia_init();

    /* Choose the initial domain {xr[0..D-1]} and get its width {xw[0..D-1]}: */
    Interval xr[D];
    Float xw[D], tol[D];
    int i;
    for(i = 0; i < D; i++)
      { xr[i] = (Interval){0.0, 1.0};
        xw[i] = xr[i].hi - xr[i].lo;
        tol[i] = (Float)(0.01 * xw[i]);
      }

    /* Get the approximate function value range {fr}: */
    /* Interval fr = f.ia(xr); */
    
    /* Get the true global minimum point {sr} and value {mr}: */
    Interval sr[D];
    f.true_opt(xr, sr);
    Interval mr = f.eval_ia(sr);
    
    /* Generate EPS output: */
    eps = bb2_new_figure(outname, f.descr, xr[0], xr[1]);
    bb2_print_range_scale("X", xr[0]);
    bb2_print_range_scale("Y", xr[1]);

    /* Draw axes and plot of {F}: */
    /* (TO BE WRITTEN) */

    /* Plot the true minimum: */
    Float xsol = (sr[0].lo + sr[0].hi)/2;
    Float ysol = (sr[1].lo + sr[1].hi)/2;
    epswr_set_pen(eps, 0.5,0.0,0.0, 0.25, 0.0,0.0);
    epswr_dot(eps, xsol, ysol, 2.0, FALSE, TRUE);

    /* Optimize and plot node boxes: */
    epswr_set_pen(eps, 0.0,0.5,0.0, 0.10, 0.0,0.0);
    Interval gr;
    FBoxList R = bb_optimize(2, f.eval_ia, &(xr[0]), &(tol[0]), &gr, draw_box); 
    affirm(R != NULL, "bb_optimize returned no boxes");
    fprintf(stderr, "minimum found = [%8.4f _ %8.4f]\n", gr.lo, gr.hi);

    /* Print the true minimum: */
    fprintf(stderr, "true minimum:\n");
    bb2_print_box(sr, &mr, FALSE);

     /* Terminate plot: */
    epswr_end_figure(eps);
    return 0;
  }

epswr_figure_t *bb2_new_figure
  ( char *outname,
    char *title,
    Interval xr, 
    Interval yr
  )
  { /* Compute plot ranges from {xr,yr} with some skosh: */
    double xw = xr.hi - xr.lo;
    double xmrg = 0.03*xw;
    double xmin = xr.lo - xmrg, xmax = xr.hi + xmrg;
    
    double yw = yr.hi - yr.lo;
    double ymrg = 0.03*yw;
    double ymin = yr.lo - ymrg, ymax = yr.hi + ymrg;
   
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
    return eps;
  }

void bb2_print_range_scale(char *which, Interval r)
  {
    fprintf(stderr, "%s range  = ", which);
    ia_print(stderr, r); 
    fprintf(stderr, "\n");
  }

void bb2_print_box(Interval *xr, Interval *fr, bool_t final)
  {
    ia_print(stderr, xr[0]); 
    fprintf(stderr, " × ");
    ia_print(stderr, xr[1]); 
    fprintf(stderr, " = ");
    ia_print(stderr, *fr); 
    if (final) { fprintf(stderr, " *"); }
    fprintf(stderr, "\n");
  }
