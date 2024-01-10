/* See foi2q.h */

#include "foi2q.h"
#include "foifloat.h"
#include "interval.h"
#include "foi.h"
#include "foimisc.h"
#include "iomisc.h"
#include "zeros2.h"
#include "plt0.h"
#include <math.h>
#include <stdio.h>
#include <malloc.h>

/*** PROTOTYPES FOR INTERNAL ROUTINES ***/

void foi2q_plot_foi_zeros(
    FILE *psfile,
    FOIP ff (FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n
  );
  /* Plots roots of ff(xf,yf) = 0, using FOI arithmetic. */
  /* The input intervals are nxn cells in xd x yd. */

char *foi2q_format_parms(
    Interval xd,
    Interval yd,
    int n
  );
  /* Formats the arguments into a sctring */
  
void foi2q_foi_zeros_aux(
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
  /* Called by foi2q_plot_foi_zeros */

void foi2q_plot_interval_zeros(
    FILE *psfile,
    Interval fv (Interval x, Interval y),
    Interval xd,
    Interval yd,
    int n
  );
  /* Plots zeros of fv(xv,yv) = 0, using ordinary interval arithmetic. */
  /* The input inervals are nxn cells in xd x yd. */

void foi2q_interval_zeros_aux(
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
  /* Called by foi2q_plot_interval_zeros */

/*** IMPLEMENTATIONS ***/

void foi2q_plots(
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
    char *parmstr = foi2q_format_parms(xd, yd, n);

    if (psfile == NULL)
      { error ("foi2q_plots: can't open PostScript file"); }

    plt0_begin_file(psfile);

    plt0_begin_page(psfile, 1, xd.lo, xd.hi, yd.lo, yd.hi, n, n);
    plt0_add_caption(psfile, title);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, parmstr);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, "Ordinary interval arithmetic");

    foi2q_plot_interval_zeros(psfile, fv, xd, yd, n);
    zeros2_plot(psfile, f, xd, yd, m);
    plt0_draw_frame(psfile);

    plt0_end_page(psfile);

    plt0_begin_page(psfile, 2, xd.lo, xd.hi, yd.lo, yd.hi, n, n);
    plt0_add_caption(psfile, title);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, parmstr);
    plt0_add_caption(psfile, "");
    plt0_add_caption(psfile, "FOI arithmetic");

    foi2q_plot_foi_zeros(psfile, ff, xd, yd, n);
    zeros2_plot(psfile, f, xd, yd, m);
    plt0_draw_frame(psfile);

    plt0_end_page(psfile);

    plt0_end_file(psfile);
    fclose(psfile);
  }

char *foi2q_format_parms(
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

void foi2q_plot_foi_zeros(
    FILE *psfile,
    FOIP ff (FOIP x, FOIP y),
    Interval xd,
    Interval yd,
    int n
  )
  { int nevals = 0;
    int ncells = 0;
    char sevals[80];

    fprintf(stderr, "\nenter foi2q_plot_foi_zeros...\n");

    plt0_begin_section(psfile, "Plot of FOI-arithmetic zeros");

    foi2q_foi_zeros_aux(
      psfile, ff, xd, yd, n, &nevals, &ncells,
      0, n, 0, n
    );

    plt0_end_section(psfile);

    fprintf(stderr, " done\n");
    fprintf(stderr, "%d FOI function evaluations\n", nevals);
    sprintf(sevals, "%d FOI function evaluations", nevals);
    plt0_add_caption(psfile, sevals);
    fprintf(stderr, "%d non-empty cells\n", ncells);
    sprintf(sevals, "%d non-empty cells", ncells);
    plt0_add_caption(psfile, sevals);
    fprintf(stderr, "exit foi2q_plot_foi_zeros.\n");
  }

void foi2q_foi_zeros_aux(
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
      { plt0_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi); }
    else if ((mx == 1) && (my == 1))
      { (*ncells)++;
        plt0_fill_and_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi, gray);
      }
    else if (mx >= my)
      { int mx2 = mx / 2;
        foi2q_foi_zeros_aux(
          psfile, ff, xd, yd, n, nevals, ncells,
          ixlo, mx2, iylo, my
        );
        foi2q_foi_zeros_aux(
          psfile, ff, xd, yd, n, nevals, ncells,
          ixlo + mx2, mx - mx2, iylo, my
        );
      }
    else /* if (my > mx) */
      { int my2 = my / 2;
        foi2q_foi_zeros_aux(
          psfile, ff, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo, my2
        );
        foi2q_foi_zeros_aux(
          psfile, ff, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo + my2, my - my2
        );
      }

    fprintf(stderr, ")");
  }

void foi2q_plot_interval_zeros(
    FILE *psfile,
    Interval fv (Interval x, Interval y),
    Interval xd,
    Interval yd,
    int n
  )
  { int nevals = 0;
    int ncells = 0;
    char sevals[80];

    fprintf(stderr, "\nenter foi2q_plot_interval_zeros...\n");

    plt0_begin_section(psfile, "Plot of interval-arithmetic zeros");

    foi2q_interval_zeros_aux(
      psfile, fv, xd, yd, n, &nevals, &ncells,
      0, n, 0, n
    );

    plt0_end_section(psfile);

    fprintf(stderr, "\n");
    fprintf(stderr, "%d interval function evaluations\n", nevals);
    sprintf(sevals, "%d interval function evaluations", nevals);
    plt0_add_caption(psfile, sevals);
    fprintf(stderr, "%d non-empty cells\n", ncells);
    sprintf(sevals, "%d non-empty cells", ncells);
    plt0_add_caption(psfile, sevals);
    fprintf(stderr, "exit foi2q_plot_interval_zeros\n");
  }

void foi2q_interval_zeros_aux(
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
      { plt0_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi); }
    else if ((mx == 1) && (my == 1))
      { plt0_fill_and_draw_rectangle(psfile, xv.lo, xv.hi, yv.lo, yv.hi, gray);
        (*ncells)++;
      }
    else if (mx >= my)
      { int mx2 = mx / 2;
        foi2q_interval_zeros_aux(
          psfile, fv, xd, yd, n, nevals, ncells,
          ixlo, mx2, iylo, my
        );
        foi2q_interval_zeros_aux(
          psfile, fv, xd, yd, n, nevals, ncells,
          ixlo + mx2, mx - mx2, iylo, my
        );
      }
    else /* if (my > mx) */
      { int my2 = my / 2;
        foi2q_interval_zeros_aux(
          psfile, fv, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo, my2
        );
        foi2q_interval_zeros_aux(
          psfile, fv, xd, yd, n, nevals, ncells,
          ixlo, mx, iylo + my2, my - my2
        );
      }

    fprintf(stderr, ")");
  }

