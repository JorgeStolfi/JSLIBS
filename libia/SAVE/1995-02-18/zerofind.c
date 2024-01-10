/* See zerofind.h */

#include "zerofind.h"
#include <ia.h>
#include <flt.h>
#include <js.h>
#include <ioprotos.h>
#include <math.h>
#include <stdio.h>

typedef struct {
    Interval x;
    Interval y;
    zf_type type;
    int level;
  } zf_stack_entry;

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

Float compute_max_root_hi (Float xlo, Float epsilon, Float delta);
  /*
    Returns the maximum $hi$ endpoint for a root-size interval whose 
    lower endpoint is $xlo$.  
   
    (By definition, an interval $xv$ is root-size if and only if 
      $xv.hi-xv.lo <= max(delta, epsilon*max(|xv.lo|,|xv.hi|))$.
  */

void zerofind_split_interval(
    Interval *xv,
    Interval *yv,
    int lv,
    Interval *ya,
    Interval *yb,
    zf_stack_entry **topp,  /* Top-of-stack entry */
    int *nstp               /* Num of entries on stack */
  );
  /*
    Splits the (complex) interval $xv$ into at least two and at
    most four pieces, and pushes them all onto the stack.
    
    For the split, assumes that the function $F$ is contained in the
    trapezoid whose altitude is $xv$ and whose bases are $ya$ at
    $xv.lo$ and $yb$ at $xv.hi$.  
    
    Ensures that any complex intervals pushed onto the stack are
    either undivisible or no more than MINSHRINK times as wide as $xv$
    itself.
    
    Assumes that $yv$ is $ia_join(ya, yb)$ (could compute it, but since
    it is available...) */
    
void zerofind_refine_interval(Interval *xv, Interval *ya, Interval *yb, Interval *xm);
  /*
    Returns a sub-interval of $xv$ containing all roots of $F(x)$ in
    $xv$, assuming that its graph is contained in the trapezoid whose
    altitude is $xv$ and whose bases are $ya$ at $xv.lo$ and $yb$ at
    $xv.hi$. */

void zerofind_stack_interval(
    Interval *xv,
    Interval *yv,
    zf_type tv,
    int lv,
    zf_stack_entry **topp,
    int *nstp
  );
  /*
    Pushes the interval $xv$ and its properties ($yv$, $tv$, $lv$)
    onto the stack, merging it with the top entry if possible. */
    
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
  
/*** IMPLEMENTATIONS ***/

#define STACKSIZE (300)
#define MINSHRINK (0.50)
#define MINMERGE (0.25)
#define DEBUG (0)

zf_type zerofind(
    IntervalPair f (Interval *xv, void *data),
    int report (Interval *xv, Interval *yv, zf_type tv, void *data),
    void *data,
    Interval xd,
    Float epsilon,
    Float delta
  )
  {
    int nst = 0;                          /* Count of stacked intervals */
    zf_stack_entry stack[STACKSIZE];      /* Stacked intervals */
    zf_stack_entry *top = &(stack[0])-1;  /* The top-of-stack entry */

    Interval x_this;         /* The current interval */
    Interval y_this;         /* Its range of values */
    zf_type t_this;          /* Its type */
    int l_this;              /* Generation number of x_this */
    
    Interval x_prev;      /* Previous interval */
    Interval y_prev;      /* Its range of values */
    zf_type t_prev;       /* Its type */
    
    x_this = xd;
    l_this = 0;
    t_this = zf_type_complex;
    
    x_prev.lo = One; x_prev.hi = Zero;   /* No previous interval yet. */
    
    while (x_this.lo <= x_this.hi)
      {
        Float max_root_hi;
        IntervalPair yp;

        /* Stack invariant: 
           /\ top == &(stack[nst-1])
           /\_{i=1}^{nst-1} stack[i].x.hi == stack[i-1].x.lo
           /\ (nst > 0 => xv.hi == stack[nst-1].x.lo)
           /\_{i=1}^{nst-1} stack[i].level >= stack[i-1].level
           /\ (nst > 0 => l_this >= stack[nst-1].level)
        */
        
        /* filter invariant:
            /\ (x_prev.lo <= x_prev.hi) => 
                  /\ (x_prev.hi == x_this.lo)
                  /\ (t_prev != zf_type_undefined) => (y_prev.lo <= y_prev.hi)
            /\ (x_prev.lo > x_prev.hi) =>
                  /\ (xd.lo == x_this.lo))
        */

        if (t_this == zf_type_complex)
          { /* Try to resolve the type of $x_this$ by evaluating $F$ on it: */
            /* Also recompute $y_this$ in the process. */

	    max_root_hi = compute_max_root_hi (x_this.lo, epsilon, delta);

	    /* Avoid evaluating on excessively small intervals; */
            /*   if $x_this$ is too small, grow it by eating up the stack: */
	    while ((x_this.hi < max_root_hi) && (nst > 0))
	      { assert (top->x.lo == x_this.hi, "zerofind: stack not contiguous");
                if (max_root_hi >= top->x.hi)
		  { x_this.hi = top->x.hi; nst--; top--; }
		else
		  { x_this.hi = max_root_hi;
		    top->x.lo = x_this.hi;
		    top->level = l_this;
		  }
	      }

	    #if DEBUG
	      fprintf(stderr, "\n");
	      fprintf(stderr, "  x_this = "); ia_print(stderr, x_this);
	      fprintf(stderr, "  l_this = %3d", l_this);
	      fprintf(stderr, "  bits = %6.2f\n", 
		zerofind_bits_gained(x_this.hi - x_this.lo, xd.hi - xd.lo)
	      );
	    #endif

	    yp = f(&x_this, data);

	    #if DEBUG
	      fprintf(stderr, "  yp.a = "); ia_print(stderr, yp.a); fprintf(stderr, "\n");
	      fprintf(stderr, "  yp.b = "); ia_print(stderr, yp.b); fprintf(stderr, "\n");
	    #endif

	    if ((yp.a.lo > yp.a.hi) || (yp.b.lo > yp.b.hi))
	      { t_this = zf_type_undefined; 
                y_this.lo = One; y_this.hi = Zero;
              }
	    else 
	      { y_this.lo = FMIN(yp.a.lo, yp.b.lo); 
		y_this.hi = FMAX(yp.a.hi, yp.b.hi);
		if (y_this.lo > Zero)
		  { t_this = zf_type_positive; }
		else if (y_this.hi < Zero)
		  { t_this = zf_type_negative; }
		else if ((y_this.lo == Zero) && (y_this.hi == Zero))
		  { t_this = zf_type_root; }
		else
		  { /* Decide if it is small enough: */
                    Float xm;
		    ROUND_NEAR;
		    xm = Half * x_this.lo + Half * x_this.hi;
                    if ((x_this.hi <= max_root_hi) || (xm <= x_this.lo) || (xm >= x_this.hi))
                      { t_this = zf_type_root; }
                    else
                      { t_this = zf_type_complex; }
                  }
	      }
          }

        if (t_this != zf_type_complex)
          { /* Output this interval */
            if (x_prev.lo <= x_prev.hi)
              { if (t_prev == t_this)
                  { /* Merge with previous interval: */
                    x_prev.hi = x_this.hi;
                    if (t_prev != zf_type_undefined) 
                      { if (y_this.lo < y_prev.lo) y_prev.lo = y_this.lo;
                        if (y_this.hi > y_prev.hi) y_prev.hi = y_this.hi;
                      }
                  }
                else
                  { /* Report previous interval: */
                    int stop = report (&x_prev, &y_prev, t_prev, data);
                    if (stop) return (t_this);
	            x_prev = x_this;
                    y_prev = y_this;
                    t_prev = t_this;
                  }
              }
            else
              { /* Save this interval: */
                x_prev = x_this;
                y_prev = y_this;
                t_prev = t_this;
              }
          }
        else 
          { /* Split x_this into smaller pieces, and stack them: */
            zerofind_split_interval(&x_this, &y_this, l_this, &yp.a, &yp.b, &top, &nst); 
          }
        
	/* Unstack the next interval: */
	if (nst == 0)
	  { /* All done */
            x_this.lo = One; x_this.hi = Zero; 
          }
	else
	  { x_this = top->x;
	    y_this = top->y;
	    l_this = top->level;
	    t_this = top->type;
            nst--; top--;
	  }
      }

    /* Flush last sub-interval: */
    if (x_prev.lo <= x_prev.hi) 
      { (void) report(&x_prev, &y_prev, t_prev, data); }
    return(zf_type_undefined);
  }
  
void zerofind_split_interval(
    Interval *xv,
    Interval *yv,
    int lv,
    Interval *ya,
    Interval *yb,
    zf_stack_entry **topp,
    int *nstp
  )
  {
    Interval xm;
    int l_new = lv+1;
    Float xsplit;

    zerofind_refine_interval(xv, ya, yb, &xm);

    #if DEBUG
      fprintf(stderr, "  xm  = "); ia_print(stderr, xm); fprintf(stderr, "\n");
    #endif

    /* Stack high-end piece: */
    if (xm.hi < xv->hi) 
      { Interval xh, yh;
        zf_type th;
        xh.lo = xm.hi;
        xh.hi = xv->hi;
        if (yb->lo> Zero)
          { th = zf_type_positive;
            yh.lo = Zero;
            yh.hi = yb->hi;
          }
        else if (yb->hi < Zero)
          { th = zf_type_negative;
            yh.lo = yb->lo;
            yh.hi = Zero;
          }
        else
          { error("zerofind_split_interval: bad split"); }

        zerofind_stack_interval(&xh, &yh, th, l_new, topp, nstp);
      }

    /* Stack or subdivide middle piece: */
    { Float max_mid_width, xv_width;
      ROUND_UP;
      max_mid_width = (MINSHRINK * xv->hi) - (MINSHRINK * xv->lo);
      xv_width = xv->hi - xv->lo;
      if (xv_width > max_mid_width)
        { ROUND_NEAR; xsplit = Half * xm.lo + Half * xm.hi; }
      else
        { xsplit = xm.lo; }
    }
    if ((xsplit > xm.lo) && (xsplit < xm.hi))
      { /* Split middle piece: */
        Interval xt;
        xt.lo = xsplit; xt.hi = xm.hi;
        zerofind_stack_interval(&xt, yv, zf_type_complex, l_new, topp, nstp);
        xt.lo = xm.lo; xt.hi = xsplit;
        zerofind_stack_interval(&xt, yv, zf_type_complex, l_new, topp, nstp);
      }
    else
      { /* Middle piece shrunk enough, stack it whole: */
        zerofind_stack_interval(&xm, yv, zf_type_complex, l_new, topp, nstp);
      }
    
    /* Stack low-end piece: */
    if (xv->lo < xm.lo) 
      { Interval xl, yl;
        zf_type tl;
        xl.lo = xv->lo;
        xl.hi = xm.lo;
        if (ya->lo > Zero)
          { tl = zf_type_positive;
            yl.lo = Zero;
            yl.hi = ya->hi;
          }
        else if (ya->hi < Zero)
          { tl = zf_type_negative;
            yl.lo = ya->lo;
            yl.hi = Zero;
          }
        else
          { error("zerofind_split_interval: bad split"); }
        zerofind_stack_interval(&xl, &yl, tl, l_new, topp, nstp);
      }
  }

void zerofind_refine_interval(Interval *xv, Interval *ya, Interval *yb, Interval *xm)
  {
    /* Compute xm = sub-inderval of xv where y may be zero: */
    if (ya->hi < Zero)
      { assert(yb->hi >= Zero, "zerofind_refine_interval: all negative");
        xm->lo = zerofind_mix_low(xv->lo, yb->hi, xv->hi, -ya->hi);
      }
    else if (ya->lo > Zero)
      { assert(yb->lo <= Zero, "zerofind_refine_interval: all positive");
        xm->lo = zerofind_mix_low(xv->lo, -yb->lo, xv->hi, ya->lo);
      }
    else
      { xm->lo = xv->lo; }

    if (yb->hi < Zero)
      { assert(ya->hi >= Zero, "zerofind_refine_interval: all negative");
        xm->hi = zerofind_mix_high(xv->lo, -yb->hi, xv->hi, ya->hi);
      }
    else if (yb->lo > Zero)
      { assert(ya->lo <= Zero, "zerofind_refine_interval: all positive");
        xm->hi = zerofind_mix_high(xv->lo, yb->lo, xv->hi, -ya->lo);
      }
    else
      { xm->hi = xv->hi; }
  }

void zerofind_stack_interval(
    Interval *xv,
    Interval *yv,
    zf_type tv,
    int lv,
    zf_stack_entry **topp,
    int *nstp
  )
  { zf_stack_entry *top = *topp;
    int nst = *nstp;
    if (nst > 0)
      { if ((tv == top->type) && (tv != zf_type_complex))
	  { /* Merge $xv$ with top-of-stack interval */
	    top->x.lo = xv->lo;
	    if (tv != zf_type_undefined) 
              { if (yv->lo < top->y.lo) top->y.lo = yv->lo;
                if (yv->hi > top->y.hi) top->y.hi = yv->hi;
              }
            top->level = lv;
            return;
	  }
      }
    /* Couldn't merge, really push: */
    assert( nst < STACKSIZE, "zerofind: stack overflow");
    top++; nst++;
    top->x = *xv;
    top->y = *yv;
    top->level = lv;
    top->type = tv;
    (*nstp) = nst;
    (*topp) = top;
  }
  
Float compute_max_root_hi (Float lo, Float epsilon, Float delta)
  { Float hi;
  
    ROUND_DOWN; hi = lo + delta;

    if (epsilon != Zero)
      { Float hi_eps;
        if (lo >= Zero)
	  { Float c;
	    ROUND_UP; c = One - epsilon;
	    ROUND_DOWN; hi_eps = lo / c;
	  }
	else
  	  { ROUND_DOWN; hi_eps = lo + (-lo) * epsilon; }
        hi = FMAX(hi, hi_eps);
      }

    #if DEBUG
      fprintf(stderr, "compute_max_root_hi: lo = %e  hi = %e \n", lo, hi);
    #endif
    
    assert(hi > lo, "compute_max_root_hi: epsilon and delta are too small");

    return(hi);
  }

Float zerofind_mix_low (Float x1, Float w1, Float x2, Float w2)
  {
    Float sw, hdx, xr;
    assert(x1 < x2, "zerofind_mix_low: x out of order");
    assert((w1 >= Zero) && (w2 >= Zero), "zerofind_mix_low: negative weights");
    ROUND_UP; sw = w1 + w2;
    if (sw == Zero) return (x1);
    ROUND_DOWN;
    hdx = ((Half*x2) + (Half*(-x1)))*(w2/sw);
    xr = (x1 + hdx) + hdx;
    assert ((xr >= x1) && (xr <= x2), "zerofind_mix_low: bug");
    return (xr);
  }
    
Float zerofind_mix_high (Float x1, Float w1, Float x2, Float w2)
  {
    Float sw, hdx, xr;
    assert(x1 < x2, "zerofind_mix_high: x out of order");
    assert((w1 >= Zero) && (w2 >= Zero), "zerofind_mix_high: negative weights");
    ROUND_UP; sw = w1 + w2;
    if (sw == Zero) return (x2);
    ROUND_DOWN;
    hdx = ((Half*x2) + (Half*(-x1))) * (w1/sw);
    ROUND_UP;
    xr = (x2 - hdx) - hdx;
    assert ((xr >= x1) && (xr <= x2), "zerofind_mix_high: bug");
    return (xr);
  }

float zerofind_bits_gained(Float curw, Float orgw)
  {
    float curbits, orgbits;
    assert (curw <= orgw, "zerofind_bits_gained: bad arguments");
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

