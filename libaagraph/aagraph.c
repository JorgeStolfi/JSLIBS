/* See aagraph.h */
/* Last edited on 2021-06-26 02:30:59 by jstolfi */

#include <aagraph.h>
#include <affirm.h>
#include <flt.h>
#include <ia.h>
#include <aa.h>
#include <aarange.h>
#include <epswr.h>
#include <math.h>
#include <stdio.h>

void aagraph_plot_paralelograms(
    epswr_figure_t *fig,
    AAP f (AAP x),
    Interval xd,
    Interval yd,
    int n
  )
  {
    Interval xv, yvlo, yvhi;
    int xi;
    AAP xf, yf;
    double xp[4], yp[4];
    AATerm eps[1];
    double gray = 0.75;

    epswr_comment(fig, "Plot of AA graph");

    for (xi=0; xi<n; xi++)
      {
        MemP frame = aa_top();

        ROUND_DOWN;
        xv.lo = (Float)(xd.lo + ((xd.hi - xd.lo)*(double)xi)/(double)n);

        ROUND_UP;
        xv.hi = (Float)(xd.lo + ((xd.hi - xd.lo)*(double)(xi+1))/(double)n);

        xf = aa_from_interval(xv);
        affirm((xf->nterms == 1), "aagraph_plot_paralelograms: xf->nterms != 1");
        eps[0].id = ((AATermP)(xf + 1))->id;

        yf = f(xf);
        
        eps[0].coef = -One;
        yvlo = aa_range(aa_fix_eps(yf, 1, eps));
        eps[0].coef = One;
        yvhi = aa_range(aa_fix_eps(yf, 1, eps));

        if (ia_is_full(&yvlo) || ia_is_full(&yvhi)) 
          { yvlo = (Interval){xd.lo, xd.hi}; yvhi = yvlo; }

        ROUND_NEAR;
        xp[0] = xv.lo;  yp[0] = yvlo.lo;
        xp[1] = xv.hi;  yp[1] = yvhi.lo;
        xp[2] = xv.hi;  yp[2] = yvhi.hi;
        xp[3] = xv.lo;  yp[3] = yvlo.hi;
        
        epswr_set_fill_color(fig, gray,gray,gray); 
        epswr_polygon(fig, TRUE, xp, yp, 4, TRUE, TRUE, TRUE);

        aa_flush(frame);
      }

    fprintf(stderr, "\n");
  }

void aagraph_plot_boxes(
    epswr_figure_t *fig,
    AAP f (AAP x),
    Interval xd,
    Interval yd,
    int n
  )
  {
    Interval xv, yv;
    int xi;
    AAP xf, yf;
    double gray = 0.75;

    epswr_comment(fig, "Plot of AA range graph");

    for (xi=0; xi<n; xi++)
      {
        MemP frame = aa_top();

        ROUND_DOWN;
        xv.lo = (Float)(xd.lo + ((xd.hi - xd.lo)*(double)xi)/(double)n);

        ROUND_UP;
        xv.hi = (Float)(xd.lo + ((xd.hi - xd.lo)*(double)(xi+1))/(double)n);

        xf = aa_from_interval(xv);

        yf = f(xf);
        yv = aa_range(yf);

        if (ia_is_full(&yv)) { yv = (Interval){xd.lo, xd.hi}; }

        ROUND_NEAR;
        epswr_set_fill_color(fig, gray,gray,gray);
        epswr_rectangle(fig, xv.lo, xv.hi, yv.lo, yv.hi, TRUE, TRUE);

        aa_flush(frame);
      }

    fprintf(stderr, "\n");
  }

void aagraph_fill_and_draw_2d_range 
  ( epswr_figure_t *fig, 
    AAP x, 
    AAP y, 
    double R, double G, double B
  )
  {
    AATermCount xn = (x->nterms);
    AATermCount yn = (y->nterms);
    AATermCount nvmax = 2*(xn+yn);
    double xv[nvmax], yv[nvmax];
    AATermCount nv;

    epswr_comment(fig, "AA joint range");

    aa_2d_range(x, y, &nv, xv, yv);
    if (nv == 0)
      { epswr_dot(fig, (double)(x->center), (double)(y->center), 0.5, FALSE, TRUE); }
    else
      { epswr_set_fill_color(fig, R,G,B);
        epswr_polygon(fig, TRUE, xv, yv, (int)nv, TRUE, TRUE, TRUE);
      }
    
    fprintf(stderr, "\n");
  }
