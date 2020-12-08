/* See {ia_trapez.h} */
/* Last edited on 2012-12-08 23:35:19 by stolfilocal */

#include <ia_trapez.h>

#include <affirm.h>
#include <assert.h>
#include <flt.h>
#include <ia.h>
#include <pswr.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

ia_trapez_t ia_trapez_from_box(Interval *xr, Interval *yr)
  { return (ia_trapez_t) { *xr, *yr, *yr }; }

ia_trapez_t ia_trapez_from_ia_diff(Interval *xr, Interval *yr, Interval *dyr, int which)
  { Float dx = xr->hi - xr->lo;
    if (which == 1) { dx = -dx; }
    Interval yo = ia_add(*yr, ia_scale(*dyr, dx, One));
    if (which == 0)
      { return (ia_trapez_t){ *xr, *yr, yo }; }
    else
      { return (ia_trapez_t){ *xr, yo, *yr }; }
  }

void ia_trapez_print(FILE *wr, ia_trapez_t *tr)
  { int full = ia_is_full(&(tr->x));
    putc('[', wr);
    putc(' ', wr);
    ia_print_bound(wr, tr->x.lo, 0, full);
    fputs(" : ", wr);
    ia_print(wr, tr->yxlo);
    fputs(" __ ", wr);
    ia_print_bound(wr, tr->x.hi, 1, full);
    fputs(" : ", wr);
    ia_print(wr, tr->yxhi);
    putc(' ', wr);
    putc(']', wr);
  }

void ia_trapez_fill(PSStream *ps, Interval *yr, ia_trapez_t *tr, float cr, float cg, float cb)
  { int is_full = ia_is_full(&(tr->yxlo)) || ia_is_full(&(tr->yxhi));
    ROUND_NEAR;
    double xp[4], yp[4];
    xp[0] = tr->x.lo;   yp[0] = (is_full ? yr->lo : tr->yxlo.lo);
    xp[1] = tr->x.hi;   yp[1] = (is_full ? yr->lo : tr->yxhi.lo);
    xp[2] = tr->x.hi;   yp[2] = (is_full ? yr->hi : tr->yxhi.hi);
    xp[3] = tr->x.lo;   yp[3] = (is_full ? yr->hi : tr->yxlo.hi);
    pswr_set_fill_color(ps, cr,cg,cb); 
    pswr_polygon(ps, TRUE, xp, yp, 4, TRUE, FALSE, FALSE); 
  }    

void ia_trapez_draw(PSStream *ps, Interval *yr, ia_trapez_t *tr)
  { int is_full = ia_is_full(&(tr->yxlo)) || ia_is_full(&(tr->yxhi));
    ROUND_NEAR;
    double xp[4], yp[4];
    xp[0] = tr->x.lo;   yp[0] = (is_full ? yr->lo : tr->yxlo.lo);
    xp[1] = tr->x.hi;   yp[1] = (is_full ? yr->lo : tr->yxhi.lo);
    xp[2] = tr->x.hi;   yp[2] = (is_full ? yr->hi : tr->yxhi.hi);
    xp[3] = tr->x.lo;   yp[3] = (is_full ? yr->hi : tr->yxlo.hi);
    pswr_polygon(ps, TRUE, xp, yp, 4, FALSE, TRUE, FALSE); 
  }

ia_trapez_t ia_trapez_clip(Interval *xr, ia_trapez_t *tp)
  {
    ia_trapez_t tc;
    if ((tp->x.hi < xr->lo) || (tp->x.lo > xr->hi) || (xr->lo > xr->hi) || (tp->x.lo > tp->x.hi))
      { /* Intersection is empty: */
        tc.x = tc.yxlo = tc.yxhi = (Interval) { +1, -1 }; 
      }
    else if (ia_is_full(xr))
      { /* {xr} is the whole line, clipping is no-op: */
        tc = (*tp);
      }
    else 
      { if (ia_is_full(&(tp->x)))
          { /* The trapezoid {tp} should be an infinite horizontal band: */
            demand((tp->yxlo.lo == tp->yxhi.lo) && (tp->yxlo.hi == tp->yxhi.hi), "bad trapezoid");
            tc.x = (*xr);
            tc.yxlo = tp->yxlo; tc.yxhi = tp->yxhi;
          }
        else
          { /* The trapezoid {tp} is finite. */
            tc.x = ia_meet(*xr, tp->x);
            /* Compute the vertical side at the lo end: */
            assert(tc.x.lo >= tp->x.lo);
            if (tc.x.lo == tp->x.lo)
              { tc.yxlo = tp->yxlo; }
            else
              { tc.yxlo = ia_interp(tp->x.lo, tp->yxlo, tp->x.hi, tp->yxhi, tc.x.lo); }
            /* Compute the vertical side at the hi end: */
            assert(tc.x.hi <= tp->x.hi);
            if (tc.x.hi == tp->x.hi)
              { tc.yxhi = tp->yxhi; }
            else
              { tc.yxhi = ia_interp(tp->x.lo, tp->yxlo, tp->x.hi, tp->yxhi, tc.x.hi); }
          }
      }
    return tc;
  }
