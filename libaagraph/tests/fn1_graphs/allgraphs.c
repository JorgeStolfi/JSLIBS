/* See {allgraphs.h} */
/* Last edited on 2009-01-06 03:10:25 by stolfi */

#define _GNU_SOURCE

#include <allgraphs.h>

#include <fltgraph.h>
#include <iagraph.h>
#include <aagraph.h>
#include <flt.h>
#include <ia.h>
#include <aa.h>
#include <pswr.h>
#include <bool.h>
#include <jsstring.h>
#include <jsfile.h>

#include <math.h>
#include <stdio.h>
#include <malloc.h>

/*** INTERNAL PROTOTYPES ***/
  
char *allgraphs_format_parms(
    Float xmin, Float xmax,
    Float ymin, Float ymax,
    int nsub
  );
  
PSStream *allgraphs_open_stream(bool_t epsformat, char *fileprefix);

void allgraphs_new_plot
  ( PSStream *ps,
    char *tag,
    char *pagenum,
    char *method,
    char *function,
    char *parmstr,
    double xmin, double xmax,
    double ymin, double ymax,
    int nsub
  );

/* IMPLEMENTATIONS */

char *allgraphs_format_parms(
    Float xmin, Float xmax,
    Float ymin, Float ymax,
    int nsub
  )
  { char *s = NULL;
    asprintf(&s, 
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
  { char *parmstr = allgraphs_format_parms(xmin, xmax, ymin, ymax, nsub);
    PSStream *ps = allgraphs_open_stream(epsformat, fileprefix);
    
    Interval xd = (Interval){xmin, xmax};
    Interval yd = (Interval){ymin, ymax};

    /*** IA BOX ENCLOSURES ***/
    allgraphs_new_plot
      ( ps, "ia", "ia", "Box enclosures (Interval Arith.)",
        function,  parmstr,
        xmin, xmax,  ymin, ymax,
        nsub
      );
    
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    iagraph_plot_boxes(ps, eval_ia, xd, yd, nsub);
    
    pswr_set_pen(ps, 0.5,0.0,0.0, 0.10, 0.0, 0.0);
    fltgraph_plot(ps, eval_fp, xd, yd, nsteps);
    
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    fltgraph_draw_axes(ps, xd, yd);
    
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0.0, 0.0);
    pswr_frame(ps);
      
    /*** AA BOX ENCLOSURES ***/
    allgraphs_new_plot
      ( ps, "ar", "ar", "Box enclosures (Affine Arith.)",
        function,  parmstr,
        xmin, xmax,  ymin, ymax,
        nsub
      );

    aagraph_plot_boxes(ps, eval_aa, xd, yd, nsub);
    fltgraph_draw_axes(ps, xd, yd);
    fltgraph_plot(ps, eval_fp, xd, yd, nsteps);
    pswr_frame(ps);
      
    /*** IA INTERVAL-SLOPE ENCLOSURES ***/
    if (diff_ia != NULL)
      { allgraphs_new_plot
          ( ps, "id", "ad", "Interval-slope enclosures (Interval Arith.)",
            function,  parmstr,
            xmin, xmax,  ymin, ymax,
            nsub
          );

        pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
        iagraph_plot_butterflies(ps, eval_ia, diff_ia, xd, yd, nsub);

        pswr_set_pen(ps, 0.5,0.0,0.0, 0.10, 0.0, 0.0);
        fltgraph_plot(ps, eval_fp, xd, yd, nsteps);

        pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
        fltgraph_draw_axes(ps, xd, yd);

        pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0.0, 0.0);
        pswr_frame(ps);
      }

    /*** AA BAND ENCLOSURES ***/
    allgraphs_new_plot
      ( ps, "aa", "aa", "Band enclosures (Affine Arith.)",
        function,  parmstr,
        xmin, xmax,  ymin, ymax,
        nsub
      );
    
    aagraph_plot_paralelograms(ps, eval_aa, xd, yd, nsub);
    fltgraph_draw_axes(ps, xd, yd);
    fltgraph_plot(ps, eval_fp, xd, yd, nsteps);
    pswr_frame(ps);
    
    /* FINISH */
    
    pswr_close_stream(ps);
  }

PSStream *allgraphs_open_stream(bool_t epsformat, char *fileprefix)
  { /* Actual EPS {hCanvasSize, vCanvasSize} will be set for each figure: */
    PSStream *ps = pswr_new_stream(fileprefix, NULL, epsformat, "doc", "letter", FALSE, 300.0, 300.0);
    return ps;
  }

void allgraphs_new_plot
  ( PSStream *ps,
    char *tag,
    char *pagenum,
    char *method,
    char *function,
    char *parmstr,
    double xmin, double xmax,
    double ymin, double ymax,
    int nsub
  )
  {
    /* Start a new page or new EPS figure: */
    pswr_new_canvas(ps, tag);
    
    /* Compute the plot scale so that the largest dimension is 6 inches: */
    double wx = xmax - xmin;
    double wy = ymax - ymin;
    double wm = (wx > wy ? wx : wy);
    double scale = 6.0 * 72.0 / wm; /* Plot scale (pt per client unit). */
    
    /* Compute the canvas dimensions (excluding the margins): */
    double hsize = scale*wx; /* Canvas H size. */
    double vsize = scale*wy; /* Canvas V size. */
    double mrg = 2.0; /* In pt. */

    /* Set the canvas size (effective only for EPS figures): */
    pswr_set_canvas_size(ps, hsize+2*mrg, vsize + 2*mrg);
    
    /* Set the client window size: */
    pswr_new_picture(ps, xmin, xmax, ymin, ymax);
    
    /* Set the grid size: */
    pswr_set_grid(ps, nsub, 1);
    
    if (! pswr_is_eps(ps))
      { /* pswr_add_caption(ps, "", 0.0); */
        pswr_add_caption(ps, function, 0.0);
        /* pswr_add_caption(ps, "", 0.0); */
        pswr_add_caption(ps, method, 0.0);
        /* pswr_add_caption(ps, "", 0.0); */
        pswr_add_caption(ps, parmstr, 0.0);
      }
  }
