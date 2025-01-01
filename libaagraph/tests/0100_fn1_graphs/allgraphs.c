/* See {allgraphs.h} */
/* Last edited on 2024-12-31 01:22:42 by stolfi */

#include <allgraphs.h>

#include <fltgraph.h>
#include <iagraph.h>
#include <aagraph.h>
#include <flt.h>
#include <ia.h>
#include <aa.h>
#include <epswr.h>
#include <bool.h>
#include <jsstring.h>
#include <jsfile.h>
#include <jsprintf.h>

#include <math.h>
#include <stdio.h>
#include <malloc.h>

/*** INTERNAL PROTOTYPES ***/
  
char *allgraphs_format_parms(
    Float xmin, Float xmax,
    Float ymin, Float ymax,
    int nsub
  );
  
epswr_figure_t *allgraphs_open_stream(bool_t epsformat, char *fileprefix);

epswr_figure_t *allgraphs_new_plot
  ( char *fileprefix, 
    char *filetag,
    char *method,
    char *function,
    char *parmstr,
    double xmin, 
    double xmax, 
    double ymin, 
    double ymax
  );

/* IMPLEMENTATIONS */

char *allgraphs_format_parms(
    Float xmin, Float xmax,
    Float ymin, Float ymax,
    int nsub
  )
  { char *s = jsprintf(
      "x in [%f _ %f]  y in [%f _ %f]\n%d intervals",
      xmin, xmax, ymin, ymax, nsub
    );
    return(s);
  }

void allgraphs_plot(
    char *fileprefix,
    int epsformat,
    char *function,
    eval_fp_t eval_fp,
    eval_ia_t eval_ia,
    diff_ia_t diff_ia,
    eval_aa_t eval_aa,
    Float xmin, Float xmax,
    Float ymin, Float ymax,
    int nsub,
    int nsteps
  )
  { 
    char *parmstr = allgraphs_format_parms(xmin, xmax, ymin, ymax, nsub);
    Interval xd = (Interval){xmin, xmax};
    Interval yd = (Interval){ymin, ymax};

    /*** IA BOX ENCLOSURES ***/

    char *method = NULL;
    epswr_figure_t *fig = NULL;
    
    /* Start a new EPS figure: */
    method = "Box enclosures (Interval Arith.)";
    fig = allgraphs_new_plot
      ( fileprefix, "ia", 
        method, function, parmstr,
        xmin, xmax,  ymin, ymax
      );
    
    epswr_set_pen(fig, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    iagraph_plot_boxes(fig, eval_ia, xd, yd, nsub);
    
    epswr_set_pen(fig, 0.5,0.0,0.0, 0.10, 0.0, 0.0);
    fltgraph_plot(fig, eval_fp, xd, yd, nsteps);
    
    epswr_set_pen(fig, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    fltgraph_draw_axes(fig, xd, yd);
    
    epswr_set_pen(fig, 0.0,0.0,0.0, 0.25, 0.0, 0.0);
    epswr_frame(fig);
    epswr_end_figure(fig); fig = NULL;
      
    /*** AA BOX ENCLOSURES ***/
    method = "Box enclosures (Affine Arith.)";
    fig = allgraphs_new_plot
      ( fileprefix, "ab",
        method, function,  parmstr,
        xmin, xmax,  ymin, ymax
      );

    aagraph_plot_boxes(fig, eval_aa, xd, yd, nsub);
    fltgraph_draw_axes(fig, xd, yd);
    fltgraph_plot(fig, eval_fp, xd, yd, nsteps);
    epswr_frame(fig);
    epswr_end_figure(fig); fig = NULL;
      
    /*** IA INTERVAL-SLOPE ENCLOSURES ***/
    if (diff_ia != NULL)
      { method = "Interval-slope enclosures (Interval Arith.)";
        fig = allgraphs_new_plot
          ( fileprefix, "id", 
            method, function,  parmstr,
            xmin, xmax,  ymin, ymax
          );

        epswr_set_pen(fig, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
        iagraph_plot_butterflies(fig, eval_ia, diff_ia, xd, yd, nsub);

        epswr_set_pen(fig, 0.5,0.0,0.0, 0.10, 0.0, 0.0);
        fltgraph_plot(fig, eval_fp, xd, yd, nsteps);

        epswr_set_pen(fig, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
        fltgraph_draw_axes(fig, xd, yd);

        epswr_set_pen(fig, 0.0,0.0,0.0, 0.25, 0.0, 0.0);
        epswr_frame(fig);
        epswr_end_figure(fig); fig = NULL;
      }

    /*** AA BAND ENCLOSURES ***/
    method = "Band enclosures (Affine Arith.)";
    fig = allgraphs_new_plot
      ( fileprefix, "aa",
        method, function,  parmstr,
        xmin, xmax,  ymin, ymax
      );
    
    aagraph_plot_paralelograms(fig, eval_aa, xd, yd, nsub);
    fltgraph_draw_axes(fig, xd, yd);
    fltgraph_plot(fig, eval_fp, xd, yd, nsteps);
    epswr_frame(fig);
    epswr_end_figure(fig); fig = NULL;
  }

epswr_figure_t *allgraphs_new_plot
  ( char *fileprefix, 
    char *filetag,
    char *method,
    char *function,
    char *parmstr,
    double xmin, 
    double xmax, 
    double ymin, 
    double ymax
  )
  { 
    /* Compute the plot scale so that the largest dimension is 6 inches: */
    double wx = xmax - xmin;
    double wy = ymax - ymin;
    double wm = (wx > wy ? wx : wy);
    double scale = 6.0 * 72.0 / wm; /* Plot scale (pt per client unit). */
    
    /* Compute the canvas dimensions (excluding the margins): */
    double hsize = scale*wx; /* Figure H size (pt). */
    double vsize = scale*wy; /* Figure V size (pt). */
    double fontsize = 10.0;

    char *fname = jsprintf("%s%s.eps", fileprefix, filetag);
    FILE *wr = open_write(fname, TRUE);
    
    double smrg = 5;
    double tmrg = 2*smrg + 5*fontsize + smrg;
    epswr_figure_t *fig = epswr_new_figure(wr, hsize, vsize, smrg,smrg,smrg,tmrg, TRUE);
    epswr_set_text_font(fig, "Courier", fontsize);
    
    /* Set the client window size: */
    epswr_set_client_window(fig, xmin, xmax, ymin, ymax);

    /* Define text margins: */
    double vtext = vsize + 2*smrg + 5*fontsize; 
    epswr_set_text_geometry(fig, FALSE, 0, hsize, vsize, vtext, 0);
    
    epswr_text(fig, function, FALSE, 0.0, TRUE, FALSE);
    epswr_text(fig, method,   FALSE, 0.0, TRUE, FALSE);
    epswr_text(fig, parmstr,  FALSE, 0.0, TRUE, FALSE);

    free(fname);
    return fig;
  }
