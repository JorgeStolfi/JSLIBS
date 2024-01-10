/* See interval.h */
/* Last edited on 2011-05-14 10:19:04 by stolfi */

#define _ISOC99_SOURCE
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <fenv.h>
#include <fpu_control.h>

#include <jsmath.h>
#include <affirm.h>

#include <interval.h> 

#define interval_IS_EMPTY(X)   ((X).end[0] > (X).end[1])
#define interval_IS_FULL(X)    (((X).end[0] == -INF) || ((X).end[1] == +INF))
#define interval_IS_FINITE(X)  (((X).end[0] > -INF) && ((X).end[1] < +INF))
#define interval_IS_TRIVIAL(X) ((X).end[0] == (X).end[1])

double interval_is_empty(interval_t *X)
  { return interval_IS_EMPTY(*X); }

double interval_is_full(interval_t *X)
  { return interval_IS_FULL(*X); }

double interval_is_finite(interval_t *X)
  { return interval_IS_FINITE(*X); }

double interval_is_trivial(interval_t *X)
  { return interval_IS_TRIVIAL(*X); }

void interval_mid_rad (interval_t *X, double *mid, double *rad)
  { if (interval_IS_FULL(*X))
      { if (mid != NULL) { (*mid) = 0; } if (rad != NULL) { (*rad) = +INF; } }
    else if (interval_IS_EMPTY(*X))
      { if (mid != NULL) { (*mid) = 0; } if (rad != NULL) { (*rad) = -INF; } }
    else if (interval_IS_TRIVIAL(*X))
      { if (mid != NULL) { (*mid) = LO(*X); } if (rad != NULL) { (*rad) = 0; } }
    else if (LO(*X) == -INF)
      { if (mid != NULL) { (*mid) = -INF; } if (rad != NULL) { (*rad) = +INF; } }
    else if (HI(*X) == +INF)
      { if (mid != NULL) { (*mid) = +INF; } if (rad != NULL) { (*rad) = +INF; } }
    else
      { double m;
        int oround = fegetround();
        fesetround(FE_TONEAREST);
        /* We divide before adding, to avoid overflow: */
        m = (LO(*X) * 0.5) + (HI(*X) * 0.5);
        if ((m < LO(*X)) || (m > HI(*X)))
          { /* Underflow must have occurred at both ends, which must have the same sign. */
            /* Adding before dividing should not overflow and give a middling value: */
            m = 0.5 * (LO(*X) + HI(*X));
            affirm((m >= LO(*X)) && (m <= HI(*X)), "rounding failed");
          }
        fesetround(oround);
        if (mid != NULL) { (*mid) = m; }
        if (rad != NULL)
          { fesetround(FE_UPWARD);
            double rlo = m - LO(*X);
            double rhi = HI(*X) - m;
            affirm((rlo >= 0.0) && (rhi >= 0.0), "rounding failed");
            (*rad) = (rlo > rhi ? rlo: rhi);
          }
        fesetround(oround);
      }
  }

double interval_mid (interval_t *X)
  { double m;
    interval_mid_rad(X, &m, NULL);
    return m;
  }

double interval_rad (interval_t *X)
  { double r;
    interval_mid_rad(X, NULL, &r);
    return r;
  }

interval_t interval_from_mid_rad (double mid, double rad)
  { 
    if (rad < 0) 
      { /* Empty interval: */
        return (interval_t) {{ +INF, -INF }};
      }
    else if (rad == 0)
      { /* Trivial interval: */
        return (interval_t) {{ mid, mid }};
      }
    else if (rad == +INF)
      { /* Even if {mid} is {±INF}, the only safe answer is the full interval: */
        return (interval_t) {{ -INF, +INF }};
      }
    else if (mid == -INF)
      { /* Since {rad} is finite, the interval must be the singleton {-INF}: */
        return (interval_t) {{ -INF, -INF }};
      }
    else if (mid == +INF)
      { /* Since {rad} is finite, the interval must be the singleton {+INF}: */
        return (interval_t) {{ +INF, +INF }};
      }
    else
      { /* Finite {mid} and finite, positive {rad}. */
        /* The result is a finite interval (except for possible overflows: */
        int oround = fegetround();
        /* We must round {LO = mid - rad} down and {HI = mid + rad} up, so: */
        fesetround(FE_UPWARD);
        double nlo = rad - mid; /* {-LO}. */
        double phi = rad + mid; /* {+HI}. */
        fesetround(oround);
        return (interval_t){{ -nlo, +phi }};
      }
  }

double interval_width (interval_t *X)
  { int oround = fegetround();
    fesetround(FE_UPWARD);
    double w = HI(*X) - LO(*X);
    fesetround(oround);
    return w;
  }

interval_t interval_split(interval_t *X, interval_side_t dir)
  { 
    double mid = interval_mid(X);
    if (dir == 0)
      { return (interval_t){{ LO(*X), mid }}; }
    else
      { return (interval_t){{ mid, HI(*X) }}; }
  }

interval_t interval_join(interval_t *X, interval_t *Y)
  { double Xlo = LO(*X), Xhi = HI(*X);
    double Ylo = LO(*Y), Yhi = HI(*Y);
    if (Xlo > Xhi) 
      { return *Y; }
    else if (Ylo > Yhi)
      { return *X; }
    else
      { interval_t w;
        LO(w) = (Xlo < Ylo ? Xlo : Ylo);
        HI(w) = (Xhi > Yhi ? Xhi : Yhi);
        return w;
      }
  }

interval_t interval_meet(interval_t *X, interval_t *Y)
  { double Xlo = LO(*X), Xhi = HI(*X);
    double Ylo = LO(*Y), Yhi = HI(*Y);
    if (Xlo > Xhi)
      { return *X; }
    if (Ylo > Yhi)
      { return *Y; }
    else
      { interval_t w;
        LO(w) = (Xlo > Ylo ? Xlo : Ylo);
        HI(w) = (Xhi < Yhi ? Xhi : Yhi);
        return w;
      }
  }

void interval_widen(interval_t *X, double margin)
  { if (interval_IS_EMPTY(*X)) { return; }
    int oround = fegetround();
    fesetround(FE_UPWARD);
    double nlo = margin - LO(*X);
    double phi = margin + HI(*X);
    fesetround(oround);
    if (-nlo > +phi) 
      { LO(*X) = +INF; HI(*X) = -INF; }
    else
      { LO(*X) = -nlo; HI(*X) = +phi; }
  }

void interval_adjust_ratio(interval_t *X, interval_t *Y, double tx, double ty)
  { demand(! interval_IS_EMPTY(*X), "interval {X} must be non-empty");
    demand(! interval_IS_EMPTY(*Y), "interval {Y} must be non-empty");
    demand(tx > 0, "ratio {tx} must be positive");
    demand(ty > 0, "ratio {ty} must be positive");
    double wx = HI(*X) - LO(*X);
    double wy = HI(*Y) - LO(*Y);
    if (ty*wx > tx*wy)
      { double ey = ((ty/tx)*wx - wy)/2;
        LO(*Y) = LO(*Y) - ey;
        HI(*Y) = HI(*Y) + ey;
      }
    else if (ty*wx < tx*wy)
      { double ex = ((tx/ty)*wy - wx)/2;
        LO(*X) = LO(*X) - ex;
        HI(*X) = HI(*X) + ex;
      }
  }

double interval_project(interval_t *X, double y)
  {
    demand(LO(*X) <= HI(*X), "empty interval");
    if (y < LO(*X))
      { return LO(*X); }
    else if (y > HI(*X))
      { return HI(*X); }
    else 
      { return y; }
  }

