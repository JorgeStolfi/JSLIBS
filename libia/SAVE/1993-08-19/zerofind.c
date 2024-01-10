/* See zerofind.h */

#include "zerofind.h"
#include "ia.h"
#include <flt.h>
#include <js.h>
#include <math.h>
#include <stdio.h>

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

Interval zerofind_refine_interval(Interval xv, Interval ya, Interval yb);
  /*
    Returns a sub-interval of $xv$ containing all roots of $y(x) = 0$
    in $xv$, assuming that $y$ is affine, $y(xv.lo)$ lies in $ya$,
    and $y(xv.hi)$ lies in $yb$.
  */

Float zerofind_mix_low (Float x1, Float w1, Float x2, Float w2);
Float zerofind_mix_high (Float x1, Float w1, Float x2, Float w2);
  /* 
    Return lower and upper bounds for the interpolation formula 
    (x1*w1 + x2*w2)/(w1 + w2).  Requires that
    $w1 >= 0$, $w2 >= 0$, and $x1 <= x2$.
  */

float zerofind_bits_gained(Float curw, Float orgw);
  /* 
    Number of information bits gained when reducing an interval
    of width $orgw$ to one of width $curw$.
  */
  
IntervalPair zerofind_join_trapezoids(
    Interval xv1,
    IntervalPair yv1, 
    Interval xv2,
    IntervalPair yv2
  );
  /*
    Given two trapezoids $xvi$, $yvi$ for consecutive segments of 
    the graph of a curve $f$, returns a single trapezoid for the
    same curve over the union of $xv1$ and $xv2$.
  */

/*** IMPLEMENTATIONS ***/

#define IVSTACKSIZE (300)
#define MINSHRINK (0.5)
#define MINMERGE (0.25)
#define DEBUG (1)

void zerofind(
    IntervalPair f (Interval xv, void *closure),
    void report (Interval xv, IntervalPair yv, void *closure),
    void *closure,
    Interval xd,
    Float epsilon,
    Float delta
  )
  {
    Interval xv;
    int ivcount = 0;               /* Count of unresolved intervals */
    Interval ivstack[IVSTACKSIZE]; /* Unresolved intervals */
    int ivlevel[IVSTACKSIZE];      /* Generation number of ivstack[i] */
    int xvlevel;                   /* Generation number of xv */
    Interval xroot;                /* Current root interval */
    IntervalPair yroot;            /* Trapezoid containing f in xroot */
    
    xv = xd;
    xvlevel = 0;
    
    xroot.lo = One; xroot.hi = Zero;   /* No roots yet. */
    
    while (xv.lo <= xv.hi)
      {
        Float minw = delta;
        IntervalPair yv;
        int is_void, is_small, is_zero;

        /* Stack invariant: 
           /\_{i=1}^{ivcount-1} ivstack[i].hi < ivstack[i-1].lo
           /\ (ivcount > 0 => xv.hi <= ivstack[ivcount-1].lo)
           /\_{i=1}^{ivcount-1} ivlevel[i] > ivlevel[i-1]
           /\ (ivcount > 0 => xvlevel >= ivlevel[ivcount-1])
        */
        
        /* Compute minimum interval width for the current xv.lo: */
        if (xv.hi < Zero) 
          { ROUND_DOWN; 
            minw = FMAX(minw, - xv.lo * epsilon); 
          }
        else if (xv.lo > Zero)
          { Float c;
            ROUND_UP;
            c = One - epsilon;
            ROUND_DOWN; 
            minw = FMAX(minw, c * epsilon);
          }

        /* Avoid small consecutive intervals: */
        while ( 
            (ivcount > 0) &&
            (xv.hi >= ivstack[ivcount-1].lo) &&
            (xv.hi < xv.lo + minw)
          )
          { ROUND_DOWN;
            xv.hi = xv.lo + minw;
            if (xv.hi >= ivstack[ivcount-1].hi)
              { xv.hi = ivstack[ivcount-1].hi;
                ivcount--;
              }
            else
              { ivstack[ivcount-1].lo = xv.hi;
                ivlevel[ivcount-1] = xvlevel;
              }
            ROUND_DOWN; /* For the "while" test */
          }
        
        #if DEBUG
          fprintf(stderr, "\n");
          fprintf(stderr, "  x =  "); ia_print(stderr, xv);
          fprintf(stderr, "  level = %3d", xvlevel);
          fprintf(stderr, "  bits = %6.2f\n", 
            zerofind_bits_gained(xv.hi - xv.lo, xd.hi - xd.lo)
          );
        #endif

        yv = f(xv, closure);

        #if DEBUG
          fprintf(stderr, "  y0 = "); ia_print(stderr, yv.a); fprintf(stderr, "\n");
          fprintf(stderr, "  y1 = "); ia_print(stderr, yv.b); fprintf(stderr, "\n");
        #endif

        is_void = (
          (yv.a.lo > yv.a.hi) ||
          (yv.b.lo > yv.b.hi) ||
          ((yv.a.lo > Zero) && (yv.b.lo > Zero)) || 
          ((yv.a.hi < Zero) && (yv.b.hi < Zero))
        );
        
        if ( (xroot.lo <= xroot.hi) && (is_void || ( xv.lo > xroot.hi)) )
          { /* Flush current root: */
            report (xroot, yroot, closure);
            xroot.lo = One;
            xroot.hi = Zero;
          };

        { Float xm, maxx;
          ROUND_NEAR;
          xm = Half * xv.lo + Half * xv.hi;
          ROUND_DOWN;
          maxx = xv.lo + minw;
          is_small = ((xv.hi <= maxx) || (xm <= xv.lo) || (xm >= xv.hi));
        }
        
        is_zero = (
          (! is_void) &&
          (yv.a.lo == Zero) && 
          (yv.a.hi == Zero) && 
          (yv.b.lo == Zero) && 
          (yv.b.hi == Zero)
        );

        if ( is_void || is_small || is_zero )
          {
            if ( ! is_void ) 
              { /* Found another root-like interval.  */
                if (xroot.hi < xv.lo)
                  { assert (xroot.hi = xv.lo, "zerofind: didn't flush root!");
                    xroot.hi = xv.hi; 
                    yroot = zerofind_join_trapezoids(xroot, yroot, xv, yv);
                  }
                else
                  { xroot = xv; yroot = yv; }
              }
            
            /* Unstack another interval: */
            if (ivcount <= 0) 
              { 
                /* All done */
                xv.lo = 1.0; xv.hi = Zero;
              }
            else
              { 
                /* Unstack an unresolved interval: */
                ivcount--;
                xv = ivstack[ivcount];
                xvlevel = ivlevel[ivcount];
              }
          }
        else
          {
            Interval zv = zerofind_refine_interval(xv, yv.a, yv.b);
            Interval nextv;
            Float zw, maxzw;
            
            #if DEBUG
              fprintf(stderr, "  z  = "); ia_print(stderr, zv); fprintf(stderr, "\n");
            #endif

            ROUND_UP;
            maxzw = (MINSHRINK * xv.hi) - (MINSHRINK * xv.lo);
            zw = zv.hi - zv.lo;
            if (zw <= maxzw) 
              { xv = zv; xvlevel++; }
            else
              { Float zm;
                ROUND_NEAR;
                zm = Half * zv.lo + Half * zv.hi;

                if ((zm <= zv.lo) || (zm >= zv.hi))
                  { /* Can't split $zv$;  should stop at next iteration */
                    xv = zv; xvlevel++;
                  }
                else
                  {
                    /* Push the interval [zm __ zv.hi] onto the stack: */
                    xvlevel++; 
                    if (
                      (ivcount > 0) && 
                      (zv.hi >= (nextv = ivstack[ivcount-1]).lo) &&
                      ((zv.hi - zm) > MINMERGE*(nextv.hi - nextv.lo))
                    )
                      { /* Merge with previous interval */
                        ivstack[ivcount-1].lo = zm;
                        ivlevel[ivcount-1] = xvlevel;
                      }
                    else
                      { /* Can't merge, really push: */
                        assert(
                          ivcount < IVSTACKSIZE, 
                          "zerofind: stack overflow"
                        );
                        ivstack[ivcount].lo = zm;
                        ivstack[ivcount].hi = zv.hi;
                        ivlevel[ivcount] = xvlevel;
                        ivcount++;
                      }
                    xv.lo = zv.lo; xv.hi = zm;
                  }
              }
          }
      }
    /* Flush last root: */
    if (xroot.lo <= xroot.hi) { report(xroot, yroot, closure); }
  }
  
IntervalPair zerofind_join_trapezoids(
    Interval xv1,
    IntervalPair yv1, 
    Interval xv2,
    IntervalPair yv2
  )
  { 
    IntervalPair yvr;
    /* Could be improved... */
    yvr.a.lo = FMIN(FMIN(yv1.a.lo, yv1.b.lo), FMIN(yv2.a.lo, yv2.b.lo));
    yvr.a.hi = FMAX(FMAX(yv1.a.hi, yv1.b.hi), FMAX(yv2.a.hi, yv2.b.hi));
    yvr.b = yvr.a;
    return(yvr);
  }
  
Interval zerofind_refine_interval(Interval xv, Interval ya, Interval yb)
  {
    Interval zv;
    /* Compute zv = sub-inderval of xv where y may be zero: */
    if (ya.hi < Zero)
      { assert(yb.hi >= Zero, "zerofind_refine_interval: all negative!\n");
        zv.lo = zerofind_mix_low(xv.lo, yb.hi, xv.hi, -ya.hi);
      }
    else if (ya.lo > Zero)
      { assert(yb.lo <= Zero, "zerofind_refine_interval: all positive!\n");
        zv.lo = zerofind_mix_low(xv.lo, -yb.lo, xv.hi, ya.lo);
      }
    else
      { zv.lo = xv.lo; }

    if (yb.hi < Zero)
      { assert(ya.hi >= Zero, "zerofind_refine_interval: all negative!\n");
        zv.hi = zerofind_mix_high(xv.lo, -yb.hi, xv.hi, ya.hi);
      }
    else if (yb.lo > Zero)
      { assert(ya.lo <= Zero, "zerofind_refine_interval: all positive!\n");
        zv.hi = zerofind_mix_high(xv.lo, yb.lo, xv.hi, -ya.lo);
      }
    else
      { zv.hi = xv.hi; }

    return(zv);
  }
  
Float zerofind_mix_low (Float x1, Float w1, Float x2, Float w2)
  {
    /* The code is abit complicated to avoid overflow: */
    Float hw, hdx, xr;
    assert(x1 < x2, "zerofind_mix_low: x out of order!\n");
    assert((w1 >= Zero) && (w2 >= Zero), "zerofind_mix_low: negative weights!\n");
    ROUND_UP;
    hw = (0.5 * w1) + (0.5 * w2);
    if (hw == Zero) return (x1);
    ROUND_DOWN;
    hdx = ((0.5 * x2) + (0.5 * (-x1))) * ((0.5 * w2) / hw);
    ROUND_DOWN;
    xr = (x1 + hdx) + hdx;
    assert ((xr >= x1) && (xr <= x2), "zerofind_mix_low: bug!\n");
    return (xr);
  }
    
Float zerofind_mix_high (Float x1, Float w1, Float x2, Float w2)
  {
    /* The code is abit complicated to avoid overflow: */
    Float hw, hdx, xr;
    assert(x1 < x2, "zerofind_mix_high: x out of order!\n");
    assert((w1 >= Zero) && (w2 >= Zero), "zerofind_mix_high: negative weights!\n");
    ROUND_UP;
    hw = (0.5 * w1) + (0.5 * w2);
    if (hw == Zero) return (x2);
    ROUND_DOWN;
    hdx = ((0.5 * x2) + (0.5 * (-x1))) * ((0.5 * w1) / hw);
    ROUND_UP;
    xr = (x2 - hdx) - hdx;
    assert ((xr >= x1) && (xr <= x2), "zerofind_mix_low: bug!\n");
    return (xr);
  }

float zerofind_bits_gained(Float curw, Float orgw)
  {
    float curbits, orgbits;
    assert (curw <= orgw, "aazeros_bits_gained: bad arguments");
    ROUND_NEAR;
    if (orgw == Zero) 
      { return(Zero); }
    else
      { orgbits = -(float)log2((double)orgw); }
    if (curw == Zero) 
      { curbits = (float)log2((double)MaxFloat); }
    else
      { curbits = -(float)log2((double)curw); }
    return(curbits - orgbits);
  }

