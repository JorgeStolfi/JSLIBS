/* See {iagraph.h} */
/* Last edited on 2016-04-01 01:17:11 by stolfilocal */

#include <iagraph.h>

#include <flt.h>
#include <ia.h>
#include <ia_butfly.h>
#include <pswr.h>
#include <affirm.h>

#include <math.h>
#include <stdio.h>

#define DEBUG 0

void iagraph_plot_boxes
  ( PSStream *ps,
    Interval f (Interval x),
    Interval xd,
    Interval yd,
    int n
  )
  { Interval xv, yv;
    int xi;
    double gray = 0.75;

    pswr_comment(ps, "Function plot with IA box enclosures");

    for (xi=0; xi<n; xi++)
      {
        ROUND_DOWN;
        xv.lo = xd.lo + ((xd.hi - xd.lo)*((float)xi))/((float)n);

        ROUND_UP;
        xv.hi = xd.lo + ((xd.hi - xd.lo)*((float)(xi+1)))/((float)n);

        yv = f(xv);

        ROUND_NEAR;
        
        iagraph_fill_and_draw_box(ps, xv, yv, xd, yd, gray,gray,gray);
      }

    fprintf(stderr, "\n");
  }

void iagraph_plot_butterflies
  ( PSStream *ps,
    Interval f (Interval x),
    Interval df (Interval x),
    Interval xd,
    Interval yd,
    int n
  )
  { Interval xv;
    int xi;
    double gray = 0.75;

    pswr_comment(ps, "Function plot IA interval-slope enclosures");

    for (xi=0; xi<n; xi++)
      {
        ROUND_DOWN;
        xv.lo = xd.lo + ((xd.hi - xd.lo)*((float)xi))/((float)n);

        ROUND_UP;
        xv.hi = xd.lo + ((xd.hi - xd.lo)*((float)(xi+1)))/((float)n);
        
        ROUND_NEAR;

        Float xmd = ia_mid(xv);
        ia_butfly_t bt;
        iagraph_compute_butterfly(xv, xmd, f, df, &bt);
        
        iagraph_fill_and_draw_butterfly(ps, &bt, xd, yd, gray,gray,gray);
      }

    fprintf(stderr, "\n");
  }

void iagraph_compute_butterfly
  ( Interval xv, 
    Float xm, 
    Interval f (Interval x),
    Interval df (Interval x),
    ia_butfly_t *bt
  )
  { demand((xm >= xv.lo) && (xm <= xv.hi), "xm outside xv");
    
    /* Set up X ranges: */
    bt->tp[0].x = (Interval){ xv.lo, xm };
    bt->tp[1].x = (Interval){ xm, xv.hi };
    
    /* Compute slope range over whole range {xv}: */
    Interval dfv = df(xv);

    /* Compute Y range {fxm} at the goven point {xm} of {xv}: */
    Interval xmv = (Interval){xm, xm};
    Interval fxmv = f(xmv);

    #if DEBUG
      fprintf(stderr, "  f(xm) = "); ia_print(stderr, fxmv); fprintf(stderr, "\n");
      fprintf(stderr, "  df(x) = "); ia_print(stderr, dfv); fprintf(stderr, "\n");
    #endif

    if (ia_is_full(&dfv) || ia_is_full(&fxmv))
      { /* Center interval or slope is unbounded; butterfly is unbounded too: */
        bt->tp[0].yxlo = bt->tp[0].yxhi = ia_full();
        bt->tp[1].yxlo = bt->tp[1].yxhi = ia_full();
        return;
      }
    
    /* Compute trapezoids from central value range {yxm} and slope range {dfv}: */
    bt->tp[0].yxhi = bt->tp[1].yxlo = fxmv;
    ROUND_UP;
    double dxlo = xm - xv.lo;
    double dxhi = xv.hi - xm;
    bt->tp[0].yxlo.hi = (float)(fxmv.hi + dfv.lo * (-dxlo));
    bt->tp[1].yxhi.hi = (float)(fxmv.hi + dfv.hi * (+dxhi));
    ROUND_DOWN;
    bt->tp[0].yxlo.lo = (float)(fxmv.lo + dfv.hi * (-dxlo));
    bt->tp[1].yxhi.lo = (float)(fxmv.lo + dfv.lo * (+dxhi));
  }

void iagraph_fill_and_draw_box
  ( PSStream *ps, 
    Interval xv, 
    Interval yv, 
    Interval xd,
    Interval yd,
    double R, double G, double B
  )
  { /* Clip box vertically to {yd}: */
    yv = ia_meet(yv, yd);
    /* Plot box: */
    pswr_set_fill_color(ps, R,G,B);
    pswr_rectangle(ps, xv.lo, xv.hi, yv.lo, yv.hi, TRUE, TRUE);
  }

void iagraph_fill_and_draw_trapezoid
  ( PSStream *ps, 
    ia_trapez_t *tp,
    Interval xd,
    Interval yd,
    double R, double G, double B
  )
  {
    Interval x = tp->x;
    if (x.hi < x.lo) { /* Empty X range, nothing to plot: */ return; }
    
    Interval yxlo = tp->yxlo;
    Interval yxhi = tp->yxhi; 
    
    /* Clip trapezoid vertically against {xd}. May result in hexagon. */
    if (ia_is_full(&yxlo) || ia_is_full(&yxhi)) { yxlo = yxhi = yd; }

    if 
      ( ((yxlo.lo < yd.hi) || (yxhi.lo < yd.hi)) &&
        ((yxlo.hi > yd.lo) || (yxhi.hi > yd.lo))
      ) 
      {
        /* Some part of the trapezoid is visible, plot it: */
        double xp[6], yp[6];
        int np = 0;  /* Number of plot corners. */
        xp[0] = x.lo;  yp[0] = yxlo.lo;
        xp[1] = x.hi;  yp[1] = yxhi.lo;
        xp[2] = x.hi;  yp[2] = yxhi.hi;
        xp[3] = x.lo;  yp[3] = yxlo.hi;
        np = 4;
        /* Should clip polygon against {xd × yd}: */
        pswr_set_fill_color(ps, R,G,B);  
        pswr_polygon(ps, TRUE, xp, yp, np, TRUE, TRUE, TRUE); 
      }
  }

void iagraph_fill_and_draw_butterfly
  ( PSStream *ps, 
    ia_butfly_t *bt, 
    Interval xd,
    Interval yd,
    double R, double G, double B
  )
  {
    int i;
    for (i = 0; i < 2; i++)
      { ia_trapez_t *ti = &(bt->tp[i]); 
        iagraph_fill_and_draw_trapezoid(ps, ti, xd, yd, R,G,B);
      }
  }
