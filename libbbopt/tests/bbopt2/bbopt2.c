/* Bivariate nonlinear branch-and-bound optimization with interval estimators. */
/* Last edited on 2009-01-06 04:58:35 by stolfi */

#include <bbopt.h>
#include <bbgoal.h>
#include <fbox.h>
#include <fboxheap.h>
#include <fboxlist.h>

#include <aa.h>
#include <ia.h>
#include <pswr.h>
#include <affirm.h>
#include <jsstring.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define D 2

/* INTERNAL PROTOS */

int main(int argc, char **argv);

PSStream *bb2_new_plot_stream
  ( bool_t epsfmt, 
    char *prefix,
    char *title,
    Interval xr, 
    Interval yr,
    double hsize,
    double vsize
  );
  /* Opens the output plot file and initializes scales, captions, etc.
    If {epsfmt} is true, the stream consists of encapsulated
    Postscript figures called "{outname}-{NNNNNN}.eps". If {epsfmt} is
    false, the stream is a plain Postscript document called
    "{outname}.ps", with "letter" paper size. The effective plot area
    will be {hsize} by {vsize} (in pt) plus a margin of 4 pt all around. */
    
void bb2_print_range_scale
  ( char *which,
    Interval r,
    double psize,
    double scale
  );
  /* Prints the client range {r}, the device window size {psize} (in pt),
    and the client-to-pswr scale {scale}, for axis {which} ("X" or "Y"). */
    
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
    char *prefix = NULL;
    asprintf(&prefix, "bb2-out-%s-", f.tag);

    /* These variables will be set later: */
    PSStream *ps = NULL;

    auto void draw_box(Interval *xr, Interval fr, bool_t final);

    void draw_box(Interval *xr, Interval fr, bool_t final)
      { bb2_print_box(xr, &fr, final);
        Interval x = xr[0];
        Interval y = xr[1];
        if (final)
          { pswr_set_fill_color(ps, 0.5,0.5,0.5);;
            pswr_rectangle(ps, x.lo, x.hi, y.lo, y.hi, TRUE, TRUE);
            
          }
        else
          { pswr_rectangle(ps, x.lo, x.hi, y.lo, y.hi, FALSE, TRUE); }
      }

    ia_init();

    /* Choose the initial domain {xr[0..D-1]} and get its width {xw[0..D-1]}: */
    Interval xr[D];
    Float xw[D], tol[D];
    int i;
    for(i = 0; i < D; i++)
      { xr[i] = (Interval){0.0, 1.0};
        xw[i] = xr[i].hi - xr[i].lo;
        tol[i] = 0.01 * xw[i];
      }

    /* Get the approximate function value range {fr}: */
    /* Interval fr = f.ia(xr); */
    
    /* Get the true global minimum point {sr} and value {mr}: */
    Interval sr[D];
    f.true_opt(xr, sr);
    Interval mr = f.eval_ia(sr);
    
    /* Postscript figure size: */
    double hsize = 432.0;
    double vsize = 432.0;
    
    /* Generate EPS and PS output: */
    int epsfmt;
    for (epsfmt = 0; epsfmt < 2; epsfmt++)
      { 
        ps = bb2_new_plot_stream
          ( epsfmt, prefix, f.descr, 
            xr[0], xr[1], 432.0, 432.0
          );
        bb2_print_range_scale("X", xr[0], hsize, 1.0);
        bb2_print_range_scale("Y", xr[1], vsize, 1.0);
        
        /* Draw axes and plot of {F}: */
        /* (TO BE WRITTEN) */

        /* Plot the true minimum: */
        Float xsol = (sr[0].lo + sr[0].hi)/2;
        Float ysol = (sr[1].lo + sr[1].hi)/2;
        pswr_set_pen(ps, 0.5,0.0,0.0, 0.25, 0.0,0.0);
        pswr_dot(ps, xsol, ysol, 2.0, FALSE, TRUE);
        
        /* Optimize and plot node boxes: */
        pswr_set_pen(ps, 0.0,0.5,0.0, 0.10, 0.0,0.0);
        Interval gr;
        FBoxList R = bb_optimize(2, f.eval_ia, &(xr[0]), &(tol[0]), &gr, draw_box); 
        affirm(R != NULL, "bb_optimize returned no boxes");
        fprintf(stderr, "minimum found = [%8.4f _ %8.4f]\n", gr.lo, gr.hi);
        
        /* Print the true minimum: */
        fprintf(stderr, "true minimum:\n");
        bb2_print_box(sr, &mr, FALSE);

         /* Terminate plot: */
        pswr_close_stream(ps);
      }
    return 0;
  }

PSStream *bb2_new_plot_stream
  ( bool_t epsfmt, 
    char *prefix,
    char *title,
    Interval xr, 
    Interval yr,
    double hsize,
    double vsize
  )
  { /* Compute widths of {xr} and {yr}: */
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
    
    double mrg = 4.0; /* Figure margin, if EPS format. */
    PSStream *ps = pswr_new_stream(prefix, NULL, epsfmt, "doc", "letter", FALSE, hsize + 2*mrg, vsize + 2*mrg);
    pswr_new_picture(ps, xr.lo, xr.hi, yr.lo, yr.hi);
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0,0.0);
    if (! epsfmt) { pswr_add_caption(ps, title, 0.0); }
    return ps;
  }

void bb2_print_range_scale
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
