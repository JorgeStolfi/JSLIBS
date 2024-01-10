/* See foi2z.h */

#include "foi2z.h"
#include "foi.h"
#include "foifloat.h"
#include "foimisc.h"
#include "interval.h"
#include "iomisc.h"
#include "zeros2.h"
#include "plt0.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

char *foi2z_format_parms(
    Interval xd,
    Interval yd,
    int n
  );
  /* Formats the arguments into a sctring */
  
void foi2z_plot_foi_zeros(
    FILE *psfile,
    FOIP ff (FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n
  );

void foi2z_plot_interval_zeros(
    FILE *psfile,
    Interval fv (Interval x, Interval y),
    Interval xd,
    Interval yd,
    int n
  );

/*** IMPLEMENTATIONS ***/

void foi2z_plots(
    char *filename,
    char *title,
    Float f(Float x, Float y),
    Interval fv(Interval x, Interval y),
    FOIP ff(FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n,
    int m
  )
  {
    FILE *psfile = fopen(filename, "w");
    char *parmstr = foi2z_format_parms(xd, yd, n);

    if (psfile == NULL)
      { error ("foi2z_plots: can't open PostScript file"); }

    plt0_begin_file(psfile);

    plt0_begin_page(psfile, 1, xd.lo, xd.hi, yd.lo, yd.hi, n, n);
    plt0_add_caption(psfile, title);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, parmstr);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, "Ordinary interval arithmetic");

    foi2z_plot_interval_zeros(psfile, fv, xd, yd, n);
    zeros2_plot(psfile, f, xd, yd, m);
    plt0_draw_grid_lines(psfile);
    plt0_draw_frame(psfile);

    plt0_end_page(psfile);

    plt0_begin_page(psfile, 2, xd.lo, xd.hi, yd.lo, yd.hi, n, n);
    plt0_add_caption(psfile, title);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, parmstr);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, "FOI arithmetic");

    foi2z_plot_foi_zeros(psfile, ff, xd, yd, n);
    zeros2_plot(psfile, f, xd, yd, m);
    plt0_draw_grid_lines(psfile);
    plt0_draw_frame(psfile);

    plt0_end_page(psfile);

    plt0_end_file(psfile);
    fclose(psfile);
  }

char *foi2z_format_parms(
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

void foi2z_plot_foi_zeros(
    FILE *psfile,
    FOIP ff (FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n
  )
  {
    Interval xv, yv, zv;
    int xi, yi;
    FOIP xf, yf, zf;
    double gray = 0.75;
    int ncells = 0;
    char sevals[80];

    fprintf(stderr, "\nenter foi2z_plot_foi_zeros...\n");

    plt0_begin_section(psfile, "Plot of FOI-arithmetic zeros");

    for (yi=0; yi<n; yi++)
      {
        for (xi=0; xi<n; xi++)
          {
            MemP frame = foi_top();

            ROUND_DOWN;
            xv.lo = xd.lo + ((xd.hi - xd.lo)*xi)/n;
            yv.lo = yd.lo + ((yd.hi - yd.lo)*yi)/n;

            ROUND_UP;
            xv.hi = xd.lo + ((xd.hi - xd.lo)*(xi+1))/n;
            yv.hi = yd.lo + ((yd.hi - yd.lo)*(yi+1))/n;

            xf = foi_from_interval(xv);
            yf = foi_from_interval(yv);

            zf = ff(xf, yf);
            zv = foi_range(zf);

            ROUND_NEAR;
            if ((zv.lo <= Zero) && (zv.hi >= Zero))
              { plt0_fill_grid_cell(psfile, xi, yi, gray);
                ncells++;
              }

            foi_flush(frame);
          }
      }

    plt0_end_section(psfile);

    fprintf(stderr, " done\n");
    fprintf(stderr, "%d non-empty cells\n", ncells);
    sprintf(sevals, "%d non-empty cells", ncells);
    plt0_add_caption(psfile, sevals);
    fprintf(stderr, "exit foi2z_plot_foi_zeros.\n");
  }

void foi2z_plot_interval_zeros(
    FILE *psfile,
    Interval fv (Interval x, Interval y),
    Interval xd,
    Interval yd,
    int n
  )
  {
    Interval xv, yv, zv;
    int xi, yi;
    double gray = 0.75;
    int ncells = 0;
    char sevals[80];

    fprintf(stderr, "\nenter foi2z_plot_interval_zeros...\n");

    plt0_begin_section(psfile, "Plot of interval-arithmetic zeros");

    for (yi=0; yi<n; yi++)
      {
        for (xi=0; xi<n; xi++)
          {
            ROUND_DOWN;
            xv.lo = xd.lo + ((xd.hi - xd.lo)*xi)/n;
            yv.lo = yd.lo + ((yd.hi - yd.lo)*yi)/n;

            ROUND_UP;
            xv.hi = xd.lo + ((xd.hi - xd.lo)*(xi+1))/n;
            yv.hi = yd.lo + ((yd.hi - yd.lo)*(yi+1))/n;

            zv = fv(xv, yv);

            ROUND_NEAR;
            if ((zv.lo <= Zero) && (zv.hi >= Zero))
              { plt0_fill_grid_cell(psfile, xi, yi, gray);
                ncells++;
              }

          }
      }

    plt0_end_section(psfile);

    fprintf(stderr, " done\n");
    fprintf(stderr, "%d non-empty cells\n", ncells);
    sprintf(sevals, "%d non-empty cells", ncells);
    plt0_add_caption(psfile, sevals);
    fprintf(stderr, "exit foi2z_plot_interval_zeros.\n");
  }

