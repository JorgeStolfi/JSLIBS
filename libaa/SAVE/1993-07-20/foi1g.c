/* See foi1g.h */

#include "foi1g.h"
#include "fgraph.h"
#include "foi.h"
#include "foifloat.h"
#include "foimisc.h"
#include "interval.h"
#include "iomisc.h"
#include "plt0.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>

/*** INTERNAL PROTOTYPES ***/

void foi1g_plot_foi_graph(
    FILE *psfile,
    FOIP ff (FOIP x),
    Interval xd,
    Interval yd,
    int n
  );

void foi1g_plot_interval_graph(
    FILE *psfile,
    Interval fv (Interval x),
    Interval xd,
    Interval yd,
    int n
  );
  
char *foi1g_format_parms(
    Interval xd,
    Interval yd,
    int n
  );
  
/*** IMPLEMENTATIONS ***/

void foi1g_plot_foi_graph(
    FILE *psfile,
    FOIP ff (FOIP x),
    Interval xd,
    Interval yd,
    int n
  )
  {
    Interval xv, yv;
    int xi;
    FOIP xf, yf;
    double gray = 0.75;

    plt0_begin_section(psfile, "Plot of FOI-arithmetic graph");

    for (xi=0; xi<n; xi++)
      {
        MemP frame = foi_top();

        ROUND_DOWN;
        xv.lo = xd.lo + ((xd.hi - xd.lo)*xi)/n;

        ROUND_UP;
        xv.hi = xd.lo + ((xd.hi - xd.lo)*(xi+1))/n;

        xf = foi_from_interval(xv);

        yf = ff(xf);
        yv = foi_range(yf);

        if (IV_ISFULL(yv)) yv = yd;

        ROUND_NEAR;
        plt0_fill_and_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi, gray);

        foi_flush(frame);
      }

    plt0_end_section(psfile);
  }

void foi1g_plot_interval_graph(
    FILE *psfile,
    Interval fv (Interval x),
    Interval xd,
    Interval yd,
    int n
  )
  {
    Interval xv, yv;
    int xi;
    double gray = 0.75;

    plt0_begin_section(psfile, "Plot of interval-arithmetic graph");

    for (xi=0; xi<n; xi++)
      {
        ROUND_DOWN;
        xv.lo = xd.lo + ((xd.hi - xd.lo)*xi)/n;

        ROUND_UP;
        xv.hi = xd.lo + ((xd.hi - xd.lo)*(xi+1))/n;

        yv = fv(xv);

        if (IV_ISFULL(yv)) yv = yd;

        ROUND_NEAR;
        plt0_fill_and_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi, gray);

      }

    plt0_end_section(psfile);
  }
  
char *foi1g_format_parms(
    Interval xd,
    Interval yd,
    int n
  )
  {
    char *s = (char *) malloc(100);
    sprintf(s, 
      "x in [%f _ %f]  y in [%f _ %f]  %d intervals",
      xd.lo, xd.hi, yd.lo, yd.hi, n
    );
    return(s);
  }

void foi1g_plots(
    char *filename,
    char *title,
    Float f(Float x),
    Interval fv(Interval x),
    FOIP ff(FOIP x),
    Interval xd,
    Interval yd,
    int n,
    int m
  )
  {
    FILE *psfile = fopen(filename, "w");
    char *parmstr = foi1g_format_parms(xd, yd, n);

    if (psfile == NULL)
      { error ("foi1g_plots: can't open PostScript file"); }

    plt0_begin_file(psfile);

    plt0_begin_page(psfile, 1, xd.lo, xd.hi, yd.lo, yd.hi, n, 1);
    plt0_add_caption(psfile, title);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, parmstr);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, "Ordinary interval arithmetic");

    foi1g_plot_interval_graph(psfile, fv, xd, yd, n);
    fgraph_draw_axes(psfile, xd, yd);
    fgraph_plot(psfile, f, xd, yd, m);
    plt0_draw_frame(psfile);

    plt0_end_page(psfile);

    plt0_begin_page(psfile, 2, xd.lo, xd.hi, yd.lo, yd.hi, n, 1);
    plt0_add_caption(psfile, title);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, parmstr);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, "FOI arithmetic");

    foi1g_plot_foi_graph(psfile, ff, xd, yd, n);
    fgraph_draw_axes(psfile, xd, yd);
    fgraph_plot(psfile, f, xd, yd, m);
    plt0_draw_frame(psfile);

    plt0_end_page(psfile);

    plt0_end_file(psfile);
    fclose(psfile);
  }
