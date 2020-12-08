// See minu_brent.h
// Last edited on 2009-02-09 11:54:53 by stolfi

#define _GNU_SOURCE
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>

#include <minu_gen.h>
#include <minu_brent.h>

#define Eps (1.0e-6) /* Square root of relative machine precision. */
#define Phi    (1.61803398874989484821) /* Golden ratio */
#define Ihp    (0.61803398874989484821) /* 1/Phi */
#define IhpSqr (0.3819660112501051518)  /* 1/Phi^2 */
 
// INTERNAL PROTOTYPES

void mub_message (char *msg);
void mub_newline(void);
void mub_print_abscissa(char *msg, double u);
void mub_print_point(char *msg, double u, double fu);

// IMPLEMENTATIONS

void minu_brent_minimize
  (
    void *parms,
    EvalProc eval,
    double *xP,
    double *fxP,
    double *dfxP,
    double tol,
    double dist,
    double *aP,
    double *bP,
    CheckProc check,
    bool_t debug
  )
  {
    #define PrintStatus() \
      {  \
        mub_newline(); \
        mub_print_abscissa("a", a); \
        mub_print_point("v", v, fv); \
        mub_print_point("x", x, fx); \
        mub_print_point("w", w, fw); \
        mub_print_abscissa("b", b); \
        mub_newline(); \
      }

    /* The current state of the search: */
    double a = (*aP);     /* Lower endpoint of current interval */
    double b = (*bP);     /* Upper endpoint of current interval. */
    /* The three probe points: */
    double x = (*xP);     /* Abscissa of best probe point. */
    double w = x;         /* Abscissa of second-smallest probe point, or {x}. */
    double v = x;         /* Abscissa of third-smallest probe point, or {x}. */
    double fx = (*fxP);   /* Value of function at {x}. */
    double fw = fx;       /* Value of function at {w}. */
    double fv = fx;       /* Value of function at {v}. */
    double dfx = (*dfxP); /* Derivative of function at {x} (not used). */ 
    /* The steps: */
    double d = 0.0;     /* Last step size. */
    double e = 0.0;     /* Next-to-last step size. */

    demand(a <= x, "initial guess {x} is less than {a}");
    demand(x <= b, "initial guess {x} is greater than {b}");

    /* Tolerance for the stopping criterion: */
    double tol3 = tol/3.0; 

    bool_t stop = FALSE;

    /* Initial check: */
    if (check != NULL) { stop = check(parms, a, b, x, fx, 0.0,0.0, 0.0,0.0); }

    while (! stop)
    {
      /* Check the loop invariant: */
      assert(fx <= fw);
      assert(fw <= fv);
      assert(a <= x);
      assert(x <= b);

      bool_t aInf = (a == -INFINITY);
      bool_t bInf = (b == +INFINITY);

      /* If {a} or {b} is infinite, it means we are still looking for the range: */
      if (aInf){ assert((x <= w) && (w <= v)); }
      if (bInf){ assert((x >= w) && (w >= v)); }

      if (debug){ PrintStatus(); }

      /* Error tolerances: */
      double tol1 = Eps*fabs(x) + tol3;
      double tol2 = 2*tol1;

      /* Check stopping criterion: */
      if ((! aInf) && (! bInf)) { stop = (fmax(fabs(x-a), fabs(b-x)) <= tol2); }
      if (stop) { break; }

      double r;            /* Step size. */
      double Dfx,Dq;       /* Estimated derivative at `x' is `Dfx/Dq'. */
      double DDfx,DDq;     /* Estimated second derivative at `x' is `DDfx/DDq'. */

      if (fabs(e) > tol1)
        { /* Fit parabola through {v}, {x}, {w}; */
          ComputeDerivatives(v, fv, w, fw, x, fx, &Dfx,&Dq, &DDfx,&DDq);
          r = e; e = d;
        }
      else
        { Dfx = 0.0; Dq = 0; DDfx = 0.0; DDq = 0.0;
          r = 0.0;
        }

      /* Compute the displacement {d} from {x} to the next probe: */
      if (DDfx < 0.0){ Dfx = -Dfx; DDfx = -DDfx; }
      if (
         (fabs(Dfx)*DDq < 0.5*fabs(DDfx*r)*Dq) && 
         (aInf || (Dfx*DDq < DDfx*(x - a)*Dq)) && 
         (bInf || (Dfx*DDq > DDfx*(x - b)*Dq))
      ) 
        { /* A parabolic-interpolation step: */
          if (debug){ mub_message("parabolic step"); }
          d = -(Dfx/DDfx)*(DDq/Dq);
        }
      else if (bInf || aInf)
        { /* A golden-section extrapolation: */
          if (debug){ mub_message("golden extension"); }
          e = (x != w ? x - w : (x != v ? x - v : 0.0));
          d = (e == 0.0 ? dist/Phi : e*Phi);
        }
      else
        { /* A golden-section interpolation: */
          if (debug){ mub_message("golden section"); }
          e = (b - x > x - a ? b - x : a - x);
          d = IhpSqr * e;
        }
      if (debug)
        { fprintf(stderr, "  x = %+13.10f  d = %+13.10f  x+d = %+13.10f\n",x,d,x+d); }
      
      /* Make sure that {d} is not zero: */
      if (fabs(d) < tol1) { d = (d < 0.0 ? -tol1 : tol1); }
      if (debug)
        { fprintf(stderr, "  x = %+13.10f  d = %+13.10f  x+d = %+13.10f\n",x,d,x+d); }
      
      assert(a <= x);
      assert(x <= b);
      assert(((b - x) > tol2) || ((x - a) > tol2));
      assert(fabs(d) >= tol1);

      if (debug) 
        { PrintStatus();
          fprintf(stderr, "  x = %20.15f  d = %20.15g  x+d = %20.15f\n", x, d, x+d);
        }

      /* Make sure that the new probe is not too close to {x}: */
      if (fabs(d) < tol1) { d = (d > 0 ? +tol1 : -tol1); }

      /* Make sure that the new probe is not too close to {a} or {b}: */
      double u = x + d; /* Tentative splitting point. */
      if (d < 0) 
        { /* Try leaving {x+d} in {[a __ x]}, else put in {[x __ b]}: */
          if (! aInf)
            { if ((x - a) < tol2)
                { /* No space in {[a__x]}, squirt to {[x__b]}: */ 
                  assert((b - x) >= tol2);
                  d = +tol1;
                }
              else if ((u - a) < tol1)
                { /* Too close to {a}, push it away: */ d = (a - x) + tol1; }
            }
        }
      else
        { /* Try leaving {x+d} in {[x __ b]}, else put in {[x __ b]}: */
          if (! bInf)
            { if ((b - x) < tol2)
                { /* No space in {[x__b]}, squirt to {[a__x]}: */ 
                  assert((x - a) >= tol2);
                  d = -tol1;
                }
              else if ((b - u) < tol1)
                { /* Too close to {b}, push it away: */ d = (b - x) - tol1; }
            }
        }
      if (debug)
        { fprintf(stderr, "  x = %+13.10f  d = %+13.10f  x+d = %+13.10f\n",x,d,x+d); }
        
      /* Final splitting point: */
      u = x + d;
      
      /* Paranoia checks: */
      assert((u - a) >= 0.99*tol1);
      assert((b - u) >= 0.99*tol1);
      assert(fabs(u - x) >= 0.99*tol1);

      /* Do the probe: */
      double fu;  /* Function value at {u}. */
      double dfu; /* Derivative at {u} (returned but ignored). */
      if (eval(parms, u, &fu, &dfu)){ break; }

      /* update  {a}, {b}: */
      if (fu >= fx)
        { if (u < x) { a = u; } else { b = u; } }
      if (fu <= fx)
        { if (u < x) { b = x; } else { a = x; } }
      assert(a <= u);
      assert(u <= b);
      
      /* Update {x}, {v}, {w} to keep the three best points, distinct if possible: */
      if (fu <= fx)
        { /* {u,fu} replaces {x,fx}: */
          if (w != x) { v = w; fv = fw; }
          w = x; fw = fx;
          x = u; fx = fu; dfx = dfu;
          if (check != NULL) { stop = check(parms, a, b, x, fx, Dfx,Dq, DDfx,DDq); }
          if (stop) { break; }
        }
      else if (fu <= fw)
        { /* {u,fu} replaces {w,fw}: */
          assert(w != x);
          v = w; fv = fw;
          w = u; fw = fu;
        }
      else if (fu <= fv)
        { /* {u,fu} replaces {v,fv} or {w,fw}: */
          assert(v != x);
          if (w == x)
            { w = u; fw = fu; }
          else
            { v = u; fv = fu; }
        }
      else if (v == w)
        { /* {u,fu} enters the set as the third best point {v,fv}: */
          v = u; fv = fu;
        }
    }

    /* Return the output values: */
    (*aP) = a;
    (*bP) = b;
    (*xP) = x;
    (*fxP) = fx;
    (*dfxP) = dfx;
  }
  
void mub_message (char *msg)
  {
    fprintf(stderr, "%s\n", msg);
  }

void mub_newline()
  {
    fprintf(stderr, "\n");
  }

void mub_print_abscissa(char *msg, double u)
  { 
    fprintf(stderr, "    %s = %20.15g\n", msg, u);
  }

void mub_print_point(char *msg, double u, double fu)
  { 
    fprintf(stderr, "    %s = %20.15g  F(%s) = %20.15g\n", msg, u, msg, fu);
  }
