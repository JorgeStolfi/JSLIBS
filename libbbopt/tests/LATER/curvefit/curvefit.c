#define PROG_NAME "curvefit"
#define PROG_DESC "Finds the best Euclidean fit of two fragments by branch-and-bound"
#define PROG_VERS "1.0"

#define curvefit_C_COPYRIGHT "Copyright © 2004 by the State University of Campinas (UNICAMP)"
/* Last edited on 2024-12-05 10:22:01 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -curves {NAME1} {NAME2} \\\n" \
  "    -tol {TOL} -sigma {SIGMA} \\\n" \
  "    [ -maxEvals {MAX_EVALS} ] \\\n" \
  "    -outName NAME \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program bla bla bla" \
  " bla bla bla bla {X+Y} bla bla" \
  " bla {INFILE} bla \"foobar.ppm\" bla bla bla\n" \
  "\n" \
  "  Beware that bla bla bla BLEBBLE BLOB bla" \
  " bla bla bla bla.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -blabber {AMOUNT}\n" \
  "    Blabbers for that {AMOUNT}. May also bla" \
  " bla bla bla bla bla bla bla bla bla bla bla bla" \
  " bla bla bla bla bla bla bla.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  stuffa(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2004-11-02 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Option bla bla added by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " curvefit_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <fbox.h>
#include <fboxheap.h>
#include <fboxlist.h>
#include <bbopt.h>
#include <bbgoal.h>

#include <aa.h>
#include <ia.h>
#include <affirm.h>
#include <epswr.h>
#include <jsstring.h>
#include <argparser.h>
#include <r2.h>

#define D 3

typedef struct options_t 
  { char *curve1;         /* Name of first curve. */
    char *curve2;         /* Name of second curve. */
    char *outName;        /* Prefix for output files. */
    /* Plotting options: */
    bool_t epsfmt;        /* TRUE for encapsulated postscript output. */
    /* Parameters for branch-and-bound optimization: */
    double tol;           /* Positional accuracy required. */
    double sigma;         /* SD of noise in curve samples. */
    int32_t maxEvals;     /* Max number of function evaluations. */
  } options_t;
  /* Command line options passed to the program. */

/* INTERNAL PROTOS */

int32_t main(int32_t argc, char **argv);

options_t *cf_parse_options(int32_t argc, char **argv);

epswr_figure_t *cf_new_figure
  ( char *outName,
    Interval *xplot,
    Interval *yplot,
    int32_t nevals
  );
  /* Opens the output Encapsulated Postscript file, called "{outname}_{NNNNN}.eps"
    where {NNNNN} is the value {nevals} formatted with "%05d".
    Writes a caption with {nevals} under the plot. The client plot window will
    contain the rectangle {xplot × yplot}. */

void cf_print_state(Interval *xr, Interval *fr, bool_t final);

void cf_print_box(Interval *xr, Interval *fr, bool_t final);
  /* A reporting function that prints the domain box {xr[0..D-1]} and 
    its image {*fr} to {stderr}. */

r2_vec_t cf_read_curve(char *filename);
r2_t cf_barycenter(r2_vec_t cv);
double cf_radius(r2_vec_t cv, r2_t bar);

void cf_draw_curves
  ( epswr_figure_t *eps,
    r2_vec_t c1,
    r2_t b1,
    double R1,
    r2_vec_t c2,
    r2_t b2,
    double R2,
    double d,
    double t1,
    double t2
  );
  
/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { options_t *o = cf_parse_options(argc, argv);
    
    ia_init();

    /* Read the two curves, get barycenters and radii: */
    r2_vec_t c1 = cf_read_curve(o->curve1);
    r2_t b1 = cf_barycenter(c1);
    double R1 = cf_radius(c1, b1);

    r2_vec_t c2 = cf_read_curve(o->curve2);
    r2_t b2 = cf_barycenter(c2);
    double R2 = cf_radius(c2, b2);
    
    /* Choose the initial domain {xr[0..D-1]} and get its width {xw[0..D-1]}: */
    /*   xr[0] = relative X displacement of center of curve 2 rel to curve 1. */
    /*   xr[1] = rotation angle of curve 1. */
    /*   xr[2] = rotation angle of curve 2. */
    /* The angles xr[1],xr[2] are in radians times the mean radii. */
    Interval xr[D];
    Float tol[D];
    xr[0] = (Interval){0, (Float)(R1+R2)};       tol[0] = (Float)o->tol;
    xr[1] = (Interval){0.0, (Float)(2*M_PI*R1)}; tol[1] = (Float)o->tol;
    xr[2] = (Interval){0.0, (Float)(2*M_PI*R2)}; tol[2] = (Float)o->tol;

    /* Plot intervals: */
    double Rmax = (R1 > R2 ? R1 : R2);
    Float hw = (Float)(1.25*Rmax);
    Interval xplot = (Interval){ -hw, +hw };
    Interval yplot = (Interval){ -hw, +hw };
    
    /* Evaluation count: */
    int32_t nevals = 0;
    
    /* Optimize and plot node boxes: */

    auto void draw_state(Interval *xr, Interval fr, bool_t final);
    
    void draw_state(Interval *xr, Interval fr, bool_t final)
      { cf_print_state(xr, &fr, final);
        Interval d = xr[0];
        Interval t1 = ia_scale(xr[1], 1.0, (Float)R1);
        Interval t2 = ia_scale(xr[2], 1.0, (Float)R2);
        nevals++;
        if ((nevals < 20) || (nevals % 20 == 0) || final)
          { epswr_figure_t *eps = cf_new_figure(o->outName, &xplot, &yplot, nevals);
            cf_draw_curves(eps, c1, b1, R1, c2, b2, R2, d.lo, t1.lo, t2.lo);
            cf_draw_curves(eps, c1, b1, R1, c2, b2, R2, d.lo, t1.lo, t2.hi);
            cf_draw_curves(eps, c1, b1, R1, c2, b2, R2, d.lo, t1.hi, t2.lo);
            cf_draw_curves(eps, c1, b1, R1, c2, b2, R2, d.lo, t1.hi, t2.hi);
            cf_draw_curves(eps, c1, b1, R1, c2, b2, R2, d.hi, t1.lo, t2.lo);
            cf_draw_curves(eps, c1, b1, R1, c2, b2, R2, d.hi, t1.lo, t2.hi);
            cf_draw_curves(eps, c1, b1, R1, c2, b2, R2, d.hi, t1.hi, t2.lo);
            cf_draw_curves(eps, c1, b1, R1, c2, b2, R2, d.hi, t1.hi, t2.hi);
            epswr_end_figure(eps);
          }
      }

    Interval gr;
    FBoxList R = bb_optimize(3, cf_eval_fit, &(xr[0]), &(tol[0]), &gr, draw_state); 
    affirm(R != NULL, "bb_optimize returned no boxes");
    fprintf(stderr, "minimum found = [%8.4f _ %8.4f]\n", gr.lo, gr.hi);

    return 0;
  }

epswr_figure_t *cf_new_figure
  ( char *outName,
    Interval *xplot,
    Interval *yplot,
    int32_t nevals
  )
  { /* Compute plot ranges from {xplot,yr} with some skosh: */
    double xw = xplot->hi - xplot->lo;
    double xmrg = 0.03*xw;
    double xmin = xplot->lo - xmrg, xmax = xplot->hi + xmrg;
    
    double yw = yplot->hi - yplot->lo;
    double ymrg = 0.03*yw;
    double ymin = yplot->lo - ymrg, ymax = yplot->hi + ymrg;
   
    double maxSize = 150*epswr_pt_per_mm;
    int32_t capLines = 1;
    double fontHeight = 10.0;
    bool_t eps_verbose = FALSE;
    epswr_figure_t *eps = epswr_new_captioned_figure
      ( NULL, outname, NULL, nevals, NULL,
        xmin,xmax, ymin,ymax, maxSize, maxSize,
        capLines, fontHeight, eps_verbose
      );
      
    char *fnum = jsprintf("%06d evaluations", nevals);
    epswr_text(eps, fnum, FALSE, 0.5, TRUE, FALSE);
    free(fnum);
    return eps;
  }

void cf_draw_curves
  ( epswr_figure_t *eps,
    r2_vec_t c1,
    r2_t b1,
    double R1,
    r2_vec_t c2,
    r2_t b2,
    double R2,
    double d,
    double t1,
    double t2
  )
  { /* Draw curve 1 in red: */
    double x1 = -R1*d/(R1+R2);
    epswr_set_pen(eps, 0.667,0.000,0.000, 0.10, 0.0,0.0);
    cf_draw_curve(c1, b1, x1, t1);
    /* Draw curve 2 in blue: */
    double x2 = R2*d/(R1+R2);
    epswr_set_pen(eps, 0.000,0.167,1.000, 0.10, 0.0,0.0);
    cf_draw_curve(c2, b2, x2, t2);
  }

void cf_print_state(Interval *xr, Interval *fr, bool_t final)
  {
    ia_print(stderr, xr[0]); 
    fprintf(stderr, " × ");
    ia_print(stderr, xr[1]); 
    fprintf(stderr, " × ");
    ia_print(stderr, xr[2]); 
    fprintf(stderr, " = ");
    ia_print(stderr, *fr); 
    if (final) { fprintf(stderr, " *"); }
    fprintf(stderr, "\n");
  }
  
options_t *cf_parse_options(int32_t argc, char **argv)
   /*
     Parses the command line options, returns a record with their options. */
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)malloc(sizeof(options_t));

    argparser_get_keyword(pp, "-curves");
    o->curve1 = argparser_get_next(pp);
    o->curve2 = argparser_get_next(pp);

    argparser_get_keyword(pp, "-tol");
    o->tol = argparser_get_next_double(pp, 1.0e-10, 1.0e+10);

    argparser_get_keyword(pp, "-sigma");
    o->sigma = argparser_get_next_double(pp, 1.0e-10, 1.0e+10);
   
    if (argparser_keyword_present(pp, "-maxEvals"))
      { o->maxEvals = argparser_get_next_int(pp, 0, 10000000); }
    else
      { o->maxEvals = 1000; }
    
    argparser_get_keyword(pp, "-outName");
    o->outName = argparser_get_next(pp);
    
    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

