/* See interval.h */
/* Last edited on 2025-02-06 23:23:34 by stolfi */

/* We need to set these in order to get {asinh}. What a crock... */
#undef __STRICT_ANSI__
#define _ISOC99_SOURCE 1
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <fenv.h>
#include <fpu_control.h>
#include <assert.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>

#include <interval.h> 

#define interval_IS_NAN(X)     (isnan((X).end[0]) || isnan((X).end[1]))
#define interval_IS_EMPTY(X)   ((X).end[0] > (X).end[1])
#define interval_IS_FULL(X)    (((X).end[0] == -INF) || ((X).end[1] == +INF))
#define interval_IS_FINITE(X)  (((X).end[0] > -INF) && ((X).end[1] < +INF))
#define interval_IS_TRIVIAL(X) ((X).end[0] == (X).end[1])

double interval_is_empty(interval_t *X)
  { return (! interval_IS_NAN(*X)) && interval_IS_EMPTY(*X); }

double interval_is_full(interval_t *X)
  { return (! interval_IS_NAN(*X)) && interval_IS_FULL(*X); }

double interval_is_finite(interval_t *X)
  { return (! interval_IS_NAN(*X)) && interval_IS_FINITE(*X); }

double interval_is_trivial(interval_t *X)
  { return (! interval_IS_NAN(*X)) && interval_IS_TRIVIAL(*X); }

void interval_mid_rad (interval_t *X, double *mid, double *rad)
  { if (interval_IS_NAN(*X))
      { if (mid != NULL) { (*mid) = NAN; } if (rad != NULL) { (*rad) = NAN; } }
    else if (interval_IS_FULL(*X))
      { if (mid != NULL) { (*mid) = 0; } if (rad != NULL) { (*rad) = +INF; } }
    else if (interval_IS_EMPTY(*X))
      { if (mid != NULL) { (*mid) = NAN; } if (rad != NULL) { (*rad) = -INF; } }
    else if (interval_IS_TRIVIAL(*X))
      { if (mid != NULL) { (*mid) = LO(*X); } if (rad != NULL) { (*rad) = 0; } }
    else if (LO(*X) == -INF)
      { if (mid != NULL) { (*mid) = -INF; } if (rad != NULL) { (*rad) = +INF; } }
    else if (HI(*X) == +INF)
      { if (mid != NULL) { (*mid) = +INF; } if (rad != NULL) { (*rad) = +INF; } }
    else
      { double m;
        int32_t oround = fegetround();
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
            if ((rlo < 0) && (rhi < 0)) 
              { fprintf(stderr, "!! [ %24.16e  %24.16e ] m = %24.16e\n", LO(*X), HI(*X), m);
                fprintf(stderr, "!! rlo = %24.16e  rhi = %24.16e\n", rlo, rhi);
              } 
            affirm((rlo >= 0.0) && (rhi >= 0.0), "rounding failed");
            (*rad) = (rlo > rhi ? rlo: rhi);
          }
        fesetround(oround);
      }
  }

double interval_mid (interval_t *X)
  { double m;
    interval_mid_rad(X, &m, NULL);
    demand(! isnan(m), "empty interval");
    return m;
  }

double interval_rad (interval_t *X)
  { double r;
    interval_mid_rad(X, NULL, &r);
    return r;
  }

double interval_width (interval_t *X)
  { if (interval_IS_NAN(*X))
      { return NAN; }
    else if (interval_IS_FULL(*X))
      { return +INF; }
    else if (interval_IS_EMPTY(*X))
      { return -INF; }
    else if (interval_IS_TRIVIAL(*X))
      { return 0.0; }
    else if (LO(*X) == -INF)
      { return +INF; }
    else if (HI(*X) == +INF)
      { return +INF; }
    else
      { int32_t oround = fegetround();
        fesetround(FE_UPWARD);
        double w = HI(*X) - LO(*X);
        fesetround(oround);
        return w;
      }
  }

bool_t interval_equal(interval_t *A, interval_t *B)
  { bool_t na = interval_IS_NAN(*A);
    bool_t nb = interval_IS_NAN(*B);
    if (na && nb) 
      { return TRUE; }
    else if (na || nb)
      { return FALSE; }
    bool_t ea = interval_IS_EMPTY(*A);
    bool_t eb = interval_IS_EMPTY(*B);
    if (ea && eb) 
      { return TRUE; }
    else if (ea || eb)
      { return FALSE; }
    else
      { return ((LO(*A) == LO(*B)) && (HI(*A) == HI(*B))); }
  }

bool_t interval_closed_has_point(interval_t *A, double p)
  { if (interval_IS_NAN(*A) || (LO(*A) > HI(*A)))
      { return FALSE; }
    else
      { return ((LO(*A) <= p) && (p <= HI(*A))); }
  }
  
bool_t interval_open_has_point(interval_t *A, double p)
  { if (interval_IS_NAN(*A) || (LO(*A) >= HI(*A)))
      { return FALSE; }
    else
      { return ((LO(*A) < p) && (p < HI(*A))); }
  }

interval_t interval_from_mid_rad (double mid, double rad)
  { 
    if (isnan(mid) || isnan(rad)) 
      { /* Invalid interval: */
        return (interval_t) {{ NAN, NAN }};
      }
    else if (rad < 0) 
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
        int32_t oround = fegetround();
        /* We must round {LO = mid - rad} down and {HI = mid + rad} up, so: */
        fesetround(FE_UPWARD);
        double nlo = rad - mid; /* {-LO}. */
        double phi = rad + mid; /* {+HI}. */
        fesetround(oround);
        return (interval_t){{ -nlo, +phi }};
      }
  }

interval_t interval_split(interval_t *X, interval_side_t dir)
  { 
    if (interval_IS_NAN(*X)) 
      { return (interval_t){{ NAN, NAN }}; }
    double mid = interval_mid(X);
    if (dir == 0)
      { return (interval_t){{ LO(*X), mid }}; }
    else
      { return (interval_t){{ mid, HI(*X) }}; }
  }

interval_t interval_include(interval_t *X, double z)
  { if (interval_IS_NAN(*X) || isnan(z)) 
      { return (interval_t){{ NAN, NAN }}; }
    double Xlo = LO(*X), Xhi = HI(*X);
    if (Xlo > Xhi) 
      { return (interval_t){{ z, z }}; }
    else
      { interval_t w;
        LO(w) = (z < Xlo ? z : Xlo);
        HI(w) = (z > Xhi ? z : Xhi);
        return w;
      }
  }

interval_t interval_join(interval_t *X, interval_t *Y)
  { if (interval_IS_NAN(*X) || interval_IS_NAN(*Y)) 
      { return (interval_t){{ NAN, NAN }}; }
    double Xlo = LO(*X), Xhi = HI(*X);
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
  { if (interval_IS_NAN(*X) || interval_IS_NAN(*Y)) 
      { return (interval_t){{ NAN, NAN }}; }
    double Xlo = LO(*X), Xhi = HI(*X);
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
  { if (interval_IS_NAN(*X) || isnan(margin)) 
      { (*X) = (interval_t){{ NAN, NAN }}; return; }
    if (interval_IS_EMPTY(*X)) { return; }
    int32_t oround = fegetround();
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
  { if (isnan(tx) || isnan(ty)) 
      { (*X) = (interval_t){{ NAN, NAN }}; (*Y) = (interval_t){{ NAN, NAN }}; return; }
    else if (interval_IS_NAN(*X)) 
      { (*X) = (interval_t){{ NAN, NAN }}; return; }
    if (interval_IS_NAN(*Y)) 
      { (*Y) = (interval_t){{ NAN, NAN }}; return; }
    demand(! interval_IS_EMPTY(*X), "interval {X} must be non-empty");
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
    if (interval_IS_NAN(*X) || isnan(y)) { return NAN; }
    demand(LO(*X) <= HI(*X), "empty interval");
    if (y < LO(*X))
      { return LO(*X); }
    else if (y > HI(*X))
      { return HI(*X); }
    else 
      { return y; }
  }

