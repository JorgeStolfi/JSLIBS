/* See {zf.h} */
/* Last edited on 2023-02-17 18:36:39 by stolfi */

#include <zf.h>

#include <flt.h>
#include <ia.h>
#include <affirm.h>
#include <bool.h>

#include <math.h>
#include <assert.h>
#include <stdio.h>

typedef struct zf_stack_entry_t {
    Interval xr;
    Interval yr;
    zf_kind_t kind;
    int level;
  } zf_stack_entry_t;

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

Float zf_compute_max_root_hi (Float xlo, double epsilon, double delta);
  /* Returns the maximum {hi} endpoint for a root-size interval whose 
    lower endpoint is {xlo}.  
   
    (By definition, an interval {xr} is root-size if and only if 
    {xr.hi <= xr.lo + max(delta, epsilon*max(|xr.lo|,|xr.hi|))}. */

bool_t zf_output_interval
  ( Interval *x_this, 
    Interval *y_this,
    zf_kind_t *k_this,
    Interval *x_prev, 
    Interval *y_prev,
    zf_kind_t *k_prev,
    zf_report_func_t *report
  );
  /* Feeds a new interval {x_this} into the interval filter. Namely,
    if the kind {k_prev} of the buffered interval {x_prev} matches the
    kind {k_this} of the new interval {x_this}, the procedure simply
    merges {x_this} into {x_prev}, joins {y_this} into {y_prev}, and
    reports nothing. Otherwise, the procedure reports {x_prev},
    {y_prev}, and {k_prev} (if {k_prev} is not {zf_kind_undefined});
    then saves {x_this}, {y_this}, and {k_this} into {x_prev},
    {y_prev}, {k_prev}.
    
    Returns TRUE if {report} was called and returned TRUE. Note that
    this condition refers to the {x_prev} interval (if it was reported),
    not to {x_this} (which was merged into {x_prev}). */
    
void zf_split_and_stack_trapezoid
  ( ia_trapez_t *tr,
    Interval yr,
    int level,
    zf_stack_entry_t **topp,  /* Top-of-stack entry */
    int *nstp               /* Num of entries on stack */
  );
  /* Pushes back a trapezoid onto the stack, whole or in pieces.
    
    Assumes that the graph of {F} is contained both in the trapezoid
    {tr}, and in the rectangle {tr->x} by {yr}. If these data imply
    that the function cannot have a zero in the interval {tr.x},
    simply pushes the box {tr.x,yr} onto the stack, with its
    classification (either {zf_kind_positive}, {zf_kind_negative}, or
    {zf_kind_undefined}).
    
    Othwerwise splits the (mixed) interval {xr} into at least two and at
    most four pieces, and pushes them all onto the stack.
    
    Ensures that any mixed intervals pushed onto the stack are
    either undivisible or no more than MINSHRINK times as wide as {xr}
    itself. */
    
void zf_stack_interval
  ( Interval *xr,
    Interval *yr,
    zf_kind_t kind,
    int level,
    zf_stack_entry_t **topp,
    int *nstp
  );
  /* Pushes the interval {xr} and its properties ({yr}, {kind}, {level})
    onto the stack, merging it with the top entry if possible. */
    
float zf_bits_gained(Float curw, Float orgw);
  /* Number of information bits gained when reducing an interval
    of width {orgw} to one of width {curw}. */
  
double log2 (double x);  

/*** IMPLEMENTATIONS ***/

#define STACKSIZE (300)
#define MINSHRINK (0.50f)
#define MINMERGE (0.25f)
#define DEBUG (0)

zf_kind_t zf_enum_zeros
  ( zf_eval_func_t *eval,
    zf_report_func_t *report,
    Interval xd,
    double epsilon,
    double delta
  )
  { int nst = 0;                            /* Count of stacked intervals */
    zf_stack_entry_t stack[STACKSIZE];      /* Stacked intervals */
    zf_stack_entry_t *top = &(stack[0])-1;  /* The top-of-stack entry */

    Interval x_this = xd;         /* The current interval */
    Interval y_this = ia_full();  /* The range of defined values of F in {x_this} */

    zf_kind_t k_this = zf_kind_mixed;  /* The kind of {x_this}. */
    int l_this = 0;                    /* Generation number of {x_this} */
    
    Interval x_prev = (Interval) { One, Zero };   /* Previous interval */
    Interval y_prev;   /* Its range of values */
    zf_kind_t k_prev;  /* Its kind */
    
    while (x_this.lo <= x_this.hi)
      {
        /* Current interval invariant:
           /\ ((x \in x_this) /\ (F(x) != undefined)) =>
                  /\ F(x) \in y_this
                  /\ (k_this == zf_kind_positive) => F(x) > 0
                  /\ (k_this == zf_kind_negative) => F(x) < 0
                  /\ k_this != zf_kind_undefined
                  /\ k_this != zf_kind_root
        */

        /* Stack invariant: 
           /\ top == &(stack[nst-1])
           /\_{i=1}^{nst-1} stack[i].xr.hi == stack[i-1].xr.lo
           /\ (nst > 0 => xv.hi == stack[nst-1].xr.lo)
           /\_{i=1}^{nst-1} stack[i].level >= stack[i-1].level
           /\ (nst > 0 => l_this >= stack[nst-1].level)
        */
        
        /* filter invariant:
            /\ (x_prev.lo <= x_prev.hi) => 
                  /\ (x_prev.hi == x_this.lo)
                  /\ (k_prev != zf_kind_undefined) => (y_prev.lo <= y_prev.hi)
            /\ (x_prev.lo > x_prev.hi) =>
                  /\ (xd.lo == x_this.lo))
        */

        if (k_this != zf_kind_mixed)
          { /* Output this interval */
            bool_t stop = zf_output_interval
              ( &x_this, &y_this, &k_this,
                &x_prev, &y_prev, &k_prev,
                report
              );
            if (stop) { return k_this; }
          }
        else
          { /* Try to resolve the kind of {x_this} by evaluating {F} on it: */

            /* Avoid evaluating {F} on excessively small intervals; */
            /*   if {x_this} is too small, grow it by eating up the stack: */
            Float mrhi_this = zf_compute_max_root_hi (x_this.lo, epsilon, delta);
            while ((x_this.hi < mrhi_this) && (nst > 0))
              { /* Merge {x_this} with the next stack entry, or with a piece thereof: */
                affirm (top->xr.lo == x_this.hi, "zerofind: stack not contiguous");
                y_this = ia_join(y_this, top->yr);
                if (mrhi_this >= top->xr.hi)
                  { /* Eat the whole entry: */
                    x_this.hi = top->xr.hi; 
                    nst--; top--;
                  }
                else
                  { /* Carve out a piece of the next entry: */
                    x_this.hi = mrhi_this;
                    top->xr.lo = x_this.hi;
                    /* Update the top entry's level to account for the carving: */
                    if (top->level < l_this) { top->level = l_this; }
                  }
              }

#           if DEBUG
              fprintf(stderr, "\n");
              fprintf(stderr, "  x_this = "); ia_print(stderr, x_this); fprintf(stderr, "\n");
              fprintf(stderr, "  y_this = "); ia_print(stderr, y_this); fprintf(stderr, "\n");
              fprintf(stderr, "  l_this = %3d", l_this);
              fprintf(stderr, "  bits = %6.2f\n", 
                zf_bits_gained(x_this.hi - x_this.lo, xd.hi - xd.lo)
              );
#           endif
            
            /* Evaluate {f} on {x_this}: */
            Interval yt;
            ia_butfly_t bt;
            eval(&x_this, &yt, &bt);
            affirm(bt.tp[0].x.lo == x_this.lo, "butterfly does not cover lo end of x range");
            affirm(bt.tp[0].x.hi == bt.tp[1].x.lo, "non-contiguous butterfly");
            affirm(bt.tp[1].x.hi == x_this.hi, "butterfly does not cover hi end of x range");

#           if DEBUG
              fprintf(stderr, "  f(x_this) = "); ia_print(stderr, yt); 
              fprintf(stderr, "\n");
              fprintf(stderr, "  bt = ");
              ia_butfly_print(stderr, &bt, "\n       ");
              fprintf(stderr, "\n");
#           endif

            /* Update {y_this} and its classification {k_this}: */
            y_this = ia_meet(y_this, yt);
            k_this = zf_classify_interval(&x_this, &y_this, mrhi_this);
            if (k_this != zf_kind_mixed)
              { /* This interval cannot or need not be subdivided, ignore the butterfly. */
                /* We push it back on the stack, rather than calling output. */
                /* The reason is that it may get merged with other intervals. */
                zf_stack_interval(&x_this, &y_this, k_this, l_this, &top, &nst);
              }
            else 
              { /* Split {x_this} into smaller pieces, if possible, and stack them. */
                Float xm = bt.tp[0].x.hi; /* Midpoint {xm} of butterfly: */
                if (xm == x_this.lo)
                  { /* Butterfly is just the hi trapezoid, split it: */
                    assert(bt.tp[1].x.lo == x_this.lo);
                    assert(bt.tp[1].x.hi == x_this.hi);
                    zf_split_and_stack_trapezoid(&(bt.tp[1]), y_this, l_this, &top, &nst);
                  }
                else if (xm == x_this.hi)
                  { /* Butterfly is just the lo trapezoid, split it: */
                    assert(bt.tp[0].x.lo == x_this.lo);
                    assert(bt.tp[0].x.hi == x_this.hi);
                    zf_split_and_stack_trapezoid(&(bt.tp[0]), y_this, l_this, &top, &nst);
                  }
                else
                  { /* Separate the butterfly into two trapezoids, stack them: */
                    affirm(xm > x_this.lo, "{xm} ouside {x_this}");
                    affirm(xm < x_this.hi, "{xm} ouside {x_this}");

                    /* Stack hi trapezoid: */
                    Interval xb = (Interval){xm, x_this.hi};
                    Float mrhi_b = zf_compute_max_root_hi(xb.lo, epsilon, delta);
                    Interval ybt1 = ia_join(bt.tp[1].yxlo, bt.tp[1].yxhi);
                    Interval yb = ia_meet(y_this, ybt1);
                    zf_kind_t kb = zf_classify_interval(&xb, &yb, mrhi_b);
                    zf_stack_interval(&xb, &yb, kb, l_this+1, &top, &nst);

                    /* Stack lo trapezoid: */
                    Interval xa = (Interval){x_this.lo, xm};
                    Float mrhi_a = zf_compute_max_root_hi (xa.lo, epsilon, delta);
                    Interval ybt0 = ia_join(bt.tp[0].yxlo, bt.tp[0].yxhi);
                    Interval ya = ia_meet(y_this, ybt0);
                    zf_kind_t ka = zf_classify_interval(&xa, &ya, mrhi_a);
                    zf_stack_interval(&xa, &ya, ka, l_this+1, &top, &nst);
                  }
              }
          }
        
        /* Unstack the next interval: */
        if (nst == 0)
          { /* All done */
            x_this.lo = One; x_this.hi = Zero; 
          }
        else
          { x_this = top->xr;
            y_this = top->yr;

            l_this = top->level;
            k_this = top->kind;
            nst--; top--;
          }
      }

    /* Flush last sub-interval: */
    if (x_prev.lo <= x_prev.hi) 
      { (void) report(&x_prev, &y_prev, k_prev); }
    return(zf_kind_undefined);
  }
  
bool_t zf_output_interval
  ( Interval *x_this, 
    Interval *y_this,
    zf_kind_t *k_this,
    Interval *x_prev, 
    Interval *y_prev,
    zf_kind_t *k_prev,
    zf_report_func_t *report
  )
  {
    if (x_prev->lo <= x_prev->hi)
      { /* Previous interval was not empty: */
        if ((*k_prev) == (*k_this))
          { /* Merge with previous interval of same kind: */
            affirm(x_prev->hi == x_this->lo, "filter is not contiguous");
            x_prev->hi = x_this->hi;
            if ((*k_prev) != zf_kind_undefined) 
              { if (y_this->lo < y_prev->lo) y_prev->lo = y_this->lo;
                if (y_this->hi > y_prev->hi) y_prev->hi = y_this->hi;
              }
            return FALSE;
          }
        else
          { /* Report previous interval: */
            bool_t stop = report (x_prev, y_prev, (*k_prev));
            if (stop) return TRUE;
          }
      }
    /* Save this interval: */
    (*x_prev) = (*x_this);
    (*y_prev) = (*y_this);
    (*k_prev) = (*k_this);
    return FALSE;
  }

zf_kind_t zf_classify_interval
  ( Interval *xr,
    Interval *yr,
    Float max_root_hi
  )
  { if (yr->lo > yr->hi)
      { yr->lo = One; yr->hi = Zero;
        return(zf_kind_undefined); 
      }
    else 
      { if (yr->lo > Zero)
          { return(zf_kind_positive); }
        else if (yr->hi < Zero)
          { return(zf_kind_negative); }
        else if ((yr->lo == Zero) && (yr->hi == Zero))
          { return(zf_kind_root); }
        else
          { /* Decide if it is small enough: */
            Float xm;
            ROUND_NEAR;
            xm = Half * xr->lo + Half * xr->hi;
            if 
              ( (xr->hi <= max_root_hi) ||
                (xm <= xr->lo) || (xm >= xr->hi)
              )
              { return(zf_kind_root); }
            else
              { return(zf_kind_mixed); }
          }
      }
  }
  
void zf_split_and_stack_trapezoid
  ( ia_trapez_t *tr,
    Interval yr,
    int level,
    zf_stack_entry_t **topp,
    int *nstp
  )
  { /* Shrink {yr} ito tight-fit {tr}: */
    yr = ia_meet(yr, ia_join(tr->yxlo, tr->yxhi));
    
    /* Reclassify the interval: */
    zf_kind_t kt = zf_kind_mixed; /* For now. */
    if (yr.lo > yr.hi)
      { kt = zf_kind_undefined; }
    else if (yr.lo > 0)
      { kt = zf_kind_positive; }
    else if (yr.hi < 0)
      { kt = zf_kind_negative; }
    
    if (kt != zf_kind_mixed)
      { /* There are no zeros in this interval, so we can stack it whole. */
        zf_stack_interval(&(tr->x), &yr, kt, level, topp, nstp);
      }
    else
      { /* At this point the trapezoid {tr} includes zero, so it must be split. */
        affirm ((yr.lo <= Zero) && (yr.hi >= Zero), "zf_split_trapezoid: bad interval");
        Interval xm, ym;
        int l_new = level+1;
        Float xsplit;
        
        zf_refine_interval(tr, &xm, &ym);

#     if DEBUG
          fprintf(stderr, "  tr = "); ia_trapez_print(stderr, tr); fprintf(stderr, "\n");
          fprintf(stderr, "  xm  = "); ia_print(stderr, xm); 
          fprintf(stderr, "  ym  = "); ia_print(stderr, ym); 
          fprintf(stderr, "\n");
#     endif

        ym = ia_meet(ym, yr);

        /* Stack high-end piece: */
        if (xm.hi < tr->x.hi) 
          { Interval xh, yh;
            zf_kind_t kh = zf_kind_undefined;
            xh.lo = xm.hi;
            xh.hi = tr->x.hi;
            if (tr->yxhi.lo > Zero)
              { kh = zf_kind_positive;
                yh.lo = Zero;
                yh.hi = yr.hi;
              }
            else if (tr->yxhi.hi < Zero)
              { kh = zf_kind_negative;
                yh.lo = yr.lo;
                yh.hi = Zero;
              }
            else
              { fatalerror("zf_split_trapezoid: bad split"); }

            zf_stack_interval(&xh, &yh, kh, l_new, topp, nstp);
          }

        /* Stack or subdivide middle piece: */
        { Float max_mid_width, x_width;
          ROUND_UP;
          max_mid_width = (MINSHRINK * tr->x.hi) - (MINSHRINK * tr->x.lo);
          x_width = tr->x.hi - tr->x.lo;
          if (x_width > max_mid_width)
            { ROUND_NEAR; xsplit = Half * xm.lo + Half * xm.hi; }
          else
            { xsplit = xm.lo; }
        }
        if ((xsplit > xm.lo) && (xsplit < xm.hi))
          { /* Split middle piece: */
            Interval xt;
            xt.lo = xsplit; xt.hi = xm.hi;
            zf_stack_interval(&xt, &ym, zf_kind_mixed, l_new, topp, nstp);
            xt.lo = xm.lo; xt.hi = xsplit;
            zf_stack_interval(&xt, &ym, zf_kind_mixed, l_new, topp, nstp);
          }
        else
          { /* Middle piece shrunk enough, stack it whole: */
            zf_stack_interval(&xm, &ym, zf_kind_mixed, l_new, topp, nstp);
          }

        /* Stack low-end piece: */
        if (tr->x.lo < xm.lo) 
          { Interval xl, yl;
            zf_kind_t kl = zf_kind_undefined;
            xl.lo = tr->x.lo;
            xl.hi = xm.lo;
            if (tr->yxlo.lo > Zero)
              { kl = zf_kind_positive;
                yl.lo = Zero;
                yl.hi = yr.hi;
              }
            else if (tr->yxlo.hi < Zero)
              { kl = zf_kind_negative;
                yl.lo = yr. lo;
                yl.hi = Zero;
              }
            else
              { fatalerror("zf_split_trapezoid: bad split"); }
            zf_stack_interval(&xl, &yl, kl, l_new, topp, nstp);
          }
      }
  }

void zf_refine_interval(ia_trapez_t *tr, Interval *xm, Interval *ym)
  {
    /* Compute xm = sub-inderval of xv where yr may be zero: */

    /* Compute low end {xm->lo} of interval: */
    if (tr->yxlo.hi < Zero)
      { affirm(tr->yxhi.hi >= Zero, "zf_refine_interval: all negative");
        xm->lo = flt_interp_lo(tr->yxlo.hi, tr->x.lo, tr->yxhi.hi, tr->x.hi, 0);
      }
    else if (tr->yxlo.lo > Zero)
      { affirm(tr->yxhi.lo <= Zero, "zf_refine_interval: all positive");
        xm->lo = flt_interp_lo(tr->yxhi.lo, tr->x.hi, tr->yxlo.lo, tr->x.lo, 0);
      }
    else
      { xm->lo = tr->x.lo; }

    /* Compute high end {xm->hi} of interval: */
    if (tr->yxhi.hi < Zero)
      { affirm(tr->yxlo.hi >= Zero, "zf_refine_interval: all negative");
        xm->hi = flt_interp_hi(tr->yxhi.hi, tr->x.hi, tr->yxlo.hi, tr->x.lo, 0);
      }
    else if (tr->yxhi.lo > Zero)
      { affirm(tr->yxlo.lo <= Zero, "zf_refine_interval: all positive");
        xm->hi = flt_interp_hi(tr->yxlo.lo, tr->x.lo, tr->yxhi.lo, tr->x.hi, 0);
      }
    else
      { xm->hi = tr->x.hi; }
      
    /* Compute vertical range of trapezoid clipped to {xm}: */
    if  (tr->yxlo.lo == tr->yxhi.lo)
      { ym->lo = tr->yxlo.lo; }
    else
      { Float xt = (tr->yxlo.lo < tr->yxhi.lo ? xm->lo : xm->hi);
        ym->lo = flt_interp_lo(tr->x.lo, tr->yxlo.lo, tr->x.hi, tr->yxhi.lo, xt); 
      }

    if  (tr->yxlo.hi == tr->yxhi.hi)
      { ym->hi = tr->yxlo.hi; }
    else
      { Float xt = (tr->yxlo.hi > tr->yxhi.hi ? xm->lo : xm->hi);
        ym->hi = flt_interp_hi(tr->x.lo, tr->yxlo.hi, tr->x.hi, tr->yxhi.hi, xt); 
      }
  }
  
void zf_stack_interval
  ( Interval *xr,
    Interval *yr,
    zf_kind_t kind,
    int level,
    zf_stack_entry_t **topp,
    int *nstp
  )
  { zf_stack_entry_t *top = *topp;
    int nst = *nstp;
    if (nst > 0)
      { if ((kind == top->kind) && (kind != zf_kind_mixed))
          { /* Merge {xr} with top-of-stack interval */
            top->xr.lo = xr->lo;
            if (kind != zf_kind_undefined) 
              { if (yr->lo < top->yr.lo) top->yr.lo = yr->lo;
                if (yr->hi > top->yr.hi) top->yr.hi = yr->hi;
              }
            affirm(top->level <= level, "zf_stack_interval: decreasing levels");
            return;
          }
      }
    /* Couldn't merge, really push: */
    affirm( nst < STACKSIZE, "zerofind: stack overflow");
    top++; nst++;
    top->xr = *xr;
    top->yr = *yr;
    top->level = level;
    top->kind = kind;
    (*nstp) = nst;
    (*topp) = top;
  }
  
Float zf_compute_max_root_hi (Float lo, double epsilon, double delta)
  { Float hi;
  
    ROUND_DOWN; hi = (Float)(lo + delta);

    if (epsilon != Zero)
      { Float hi_eps;
        if (lo >= Zero)
          { Float c;
            ROUND_UP; c = One - (Float)epsilon;
            ROUND_DOWN; hi_eps = lo / c;
          }
        else
          { ROUND_DOWN; hi_eps = lo + (-lo) * (Float)epsilon; }
        hi = FMAX(hi, hi_eps);
      }

#   if DEBUG
      fprintf(stderr, "zf_compute_max_root_hi: lo = %e  hi = %e \n", lo, hi);
#   endif
    
    affirm(hi > lo, "zf_compute_max_root_hi: epsilon and delta are too small");

    return(hi);
  }

float zf_bits_gained(Float curw, Float orgw)
  { float curbits, orgbits;
    affirm (curw <= orgw, "zf_bits_gained: bad arguments");
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

double log2 (double x)
  { return (double)(1.4426950408889634074L * (long double)log(x)); }
