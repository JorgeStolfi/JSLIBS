/* See foi2qe.h */

#include "foi2qe.h"
#include "foifloat.h"
#include "interval.h"
#include "foi.h"
#include "foimisc.h"
#include "iomisc.h"
#include "zeros2.h"
#include "plte.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>

/*** PROTOTYPES FOR INTERNAL ROUTINES ***/

void foi2qe_plot_foi_zeros(
    FILE *psfile,
    FOIP ff (FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n
  );
  /* Plots roots of ff(xf,yf) = 0, using FOI arithmetic. */
  /* The input intervals are nxn cells in xd x yd. */

char *foi2qe_format_parms(
    Interval xd,
    Interval yd,
    int n
  );
  /* Formats the arguments into a sctring */
  
void foi2qe_foi_zeros_aux(
    FILE *psfile,
    FOIP ff (FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n,
    int *nevals,
    int *ncells,
    int ixlo, int mx,
    int iylo, int my
  );
  /* Called by foi2qe_plot_foi_zeros */

void foi2qe_plot_interval_zeros(
    FILE *psfile,
    Interval fv (Interval x, Interval y),
    Interval xd,
    Interval yd,
    int n
  );
  /* Plots zeros of fv(xv,yv) = 0, using ordinary interval arithmetic. */
  /* The input inervals are nxn cells in xd x yd. */

void foi2qe_interval_zeros_aux(
    FILE *psfile,
    Interval fv (Interval x, Interval y),
    Interval xd,
    Interval yd,
    int n,
    int *nevals,
    int *ncells,
    int ixlo, int mx,
    int iylo, int my
  );
  /* Called by foi2qe_plot_interval_zeros */

/*** IMPLEMENTATIONS ***/

void foi2qe_plots(
    char *filename1, char *filename2,
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
    {
      FILE *psfile = fopen(filename1, "w");
      if (psfile == NULL)
	{ error ("foi2qe_plots: can't open 1st PostScript file"); }

      plte_begin_figure(psfile, xd.lo, xd.hi, yd.lo, yd.hi, n, n);
      foi2qe_plot_interval_zeros(psfile, fv, xd, yd, n);
      zeros2_plot(psfile, f, xd, yd, m);
      plte_draw_frame(psfile);
      plte_end_figure(psfile);
      fclose(psfile);
    }
    
    {
      FILE *psfile = fopen(filename2, "w");
      if (psfile == NULL)
	{ error ("foi2qe_plots: can't open 2nd PostScript file"); }

      plte_begin_figure(psfile, xd.lo, xd.hi, yd.lo, yd.hi, n, n);
      foi2qe_plot_foi_zeros(psfile, ff, xd, yd, n);
      zeros2_plot(psfile, f, xd, yd, m);
      plte_draw_frame(psfile);
      plte_end_figure(psfile);
      fclose(psfile);
    }
  }

char *foi2qe_format_parms(
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

void foi2qe_plot_foi_zeros(
    FILE *psfile,
    FOIP ff (FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n
  )
  { int nevals = 0;
    int ncells = 0;

    fprintf(stderr, "\nenter foi2qe_plot_foi_zeros...\n");

    plte_begin_section(psfile, "Plot of FOI-arithmetic zeros");

    foi2qe_foi_zeros_aux(
      psfile, ff, xd, yd, n, &nevals, &ncells,
      0, n, 0, n
    );

    plte_end_section(psfile);

    fprintf(stderr, " done\n");
    fprintf(stderr, "%d FOI function evaluations\n", nevals);
    fprintf(stderr, "%d non-empty cells\n", ncells);
    fprintf(stderr, "exit foi2qe_plot_foi_zeros.\n");
  }

void foi2qe_foi_zeros_aux(
    FILE *psfile,
    FOIP ff (FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n,
    int *nevals,
    int *ncells,
    int ixlo, int mx,
    int iylo, int my
  )
  {
    Interval xv, yv, zv;
    FOIP xf, yf, zf;
    double gray = 0.75;

    MemP frame = foi_top();

    fprintf(stderr, "(");

    ROUND_DOWN;
    xv.lo = xd.lo + ((xd.hi - xd.lo)*ixlo)/n;
    yv.lo = yd.lo + ((yd.hi - yd.lo)*iylo)/n;

    ROUND_UP;
    xv.hi = xd.lo + ((xd.hi - xd.lo)*(ixlo+mx))/n;
    yv.hi = yd.lo + ((yd.hi - yd.lo)*(iylo+my))/n;

    xf = foi_from_interval(xv);
    yf = foi_from_interval(yv);

    zf = ff(xf, yf);
    zv = foi_range(zf);
    (*nevals)++;

    foi_flush(frame);

    ROUND_NEAR;
    if ((zv.lo > Zero) || (zv.hi < Zero))
      { plte_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi); }
    else if ((mx == 1) && (my == 1))
      { (*ncells)++;
        plte_fill_and_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi, gray);
      }
    else if (mx >= my)
      { int mx2 = mx / 2;
        foi2qe_foi_zeros_aux(
          psfile, ff, xd, yd, n, nevals, ncells,
          ixlo, mx2, iylo, my
        );
        foi2qe_foi_zeros_aux(
          psfile, ff, xd, yd, n, nevals, ncells,
          ixlo + mx2, mx - mx2, iylo, my
        );
      }
    else /* if (my > mx) */
      { int my2 = my / 2;
        foi2qe_foi_zeros_aux(
          psfile, ff, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo, my2
        );
        foi2qe_foi_zeros_aux(
          psfile, ff, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo + my2, my - my2
        );
      }

    fprintf(stderr, ")");
  }

void foi2qe_plot_interval_zeros(
    FILE *psfile,
    Interval fv (Interval x, Interval y),
    Interval xd,
    Interval yd,
    int n
  )
  { int nevals = 0;
    int ncells = 0;

    fprintf(stderr, "\nenter foi2qe_plot_interval_zeros...\n");

    plte_begin_section(psfile, "Plot of interval-arithmetic zeros");

    foi2qe_interval_zeros_aux(
      psfile, fv, xd, yd, n, &nevals, &ncells,
      0, n, 0, n
    );

    plte_end_section(psfile);

    fprintf(stderr, "\n");
    fprintf(stderr, "%d interval function evaluations\n", nevals);
    fprintf(stderr, "%d non-empty cells\n", ncells);
    fprintf(stderr, "exit foi2qe_plot_interval_zeros\n");
  }

void foi2qe_interval_zeros_aux(
    FILE *psfile,
    Interval fv (Interval x, Interval y),
    Interval xd,
    Interval yd,
    int n,
    int *nevals,
    int *ncells,
    int ixlo, int mx,
    int iylo, int my
  )
  {
    Interval xv, yv, zv;
    double gray = 0.75;

    fprintf(stderr, "(");

    ROUND_DOWN;
    xv.lo = xd.lo + ((xd.hi - xd.lo)*ixlo)/n;
    yv.lo = yd.lo + ((yd.hi - yd.lo)*iylo)/n;

    ROUND_UP;
    xv.hi = xd.lo + ((xd.hi - xd.lo)*(ixlo+mx))/n;
    yv.hi = yd.lo + ((yd.hi - yd.lo)*(iylo+my))/n;

    zv = fv(xv, yv);
    (*nevals)++;

    ROUND_NEAR;
    if ((zv.lo > Zero) || (zv.hi < Zero))
      { plte_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi); }
    else if ((mx == 1) && (my == 1))
      { plte_fill_and_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi, gray);
        (*ncells)++;
      }
    else if (mx >= my)
      { int mx2 = mx / 2;
        foi2qe_interval_zeros_aux(
          psfile, fv, xd, yd, n, nevals, ncells,
          ixlo, mx2, iylo, my
        );
        foi2qe_interval_zeros_aux(
          psfile, fv, xd, yd, n, nevals, ncells,
          ixlo + mx2, mx - mx2, iylo, my
        );
      }
    else /* if (my > mx) */
      { int my2 = my / 2;
        foi2qe_interval_zeros_aux(
          psfile, fv, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo, my2
        );
        foi2qe_interval_zeros_aux(
          psfile, fv, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo + my2, my - my2
        );
      }

    fprintf(stderr, ")");
  }

