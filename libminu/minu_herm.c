// See minu_herm.h
// Last edited on 2024-12-05 10:35:16 by stolfi

#include <minu_gen.h>
#include <minu_herm.h>

#include <affirm.h>

#include <values.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Eps (1.0e-6) /* Square root of relative machine precision. */
#define Phi  (1.61803398874989484821) /* Golden ratio */
#define Ihp (0.61803398874989484821) /* 1/Phi */

#define ALPHA Ihp

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

// INTERNAL PROTOTYPES

bool_t muhm_eval_and_check
  ( EvalProc eval, CheckProc check, void *parms,
    double t, double *ft, double *dft,
    double *x, double *fx, double *dfx,
    double a, double b
  );
/*
  Calls `eval' on `t' with result returned in `*ft'.
  (unless `t==*x', in which case simply copies `*fx').
  If eval returns `TRUE', returns `TRUE'; else, if `*ft<=*fx',
  updates `*x,*fx,*dfx', calls `check' (with bogus second derivative),
  and returns whathever it returns. */
  
double muhm_ensure_progress
  ( double s,
    double a, double u, double v, double b,
    double delta,
    double reqRatio,
    bool_t debug
  );
/* 
  Adjusts the estimated minimum `s' so as to ensure that
  the worst-case shrinking ratio at the next step,
  namely `t = MAX(s-u,v-s)/(v-u)', is no greter
  than `reqRatio', and that `s' is at least 
  `delta' away from the bracketing points `u,v'. */

double muhm_estimate_minimum
  ( double u, double fu, double dfu, 
    double v, double fv, double dfv
  );
/*
  Returns the estimated position of the minumum of `F(x)' in `[u__v]',
  given its values and derivatives at `u' and `v'. */

void muhm_quadratic_roots(double A, double B, double C, double *x1, double *x2);
/*
  Stores in `*x1' and `*x2' the real part of the roots of `A*x^2 + B*x + C'. */

void muhm_message(char *msg);
void muhm_newline(void);
void muhm_print_abscissa(char *msg, double u);
void muhm_print_point(char *msg, double u, double fu, double dfu, double ddfu);

// IMPLEMENTATIONS

void minu_herm_minimize
  (
    void *parms,
    EvalProc eval,
    double *x,
    double *fx,
    double *dfx,
    double tol,
    double dist,
    double *a,
    double *b,
    CheckProc check,
    bool_t debug
  )
{
  double u, fu, dfu;    /* Lower bracketing point */
  double v, fv, dfv;    /* Upper bracketing point */
    
  #define PrintStatus() \
    {  \
      muhm_newline(); \
      muhm_print_abscissa("  a", *a); \
      muhm_print_point((*x == u ? ">>u" : "  u"), u, fu, dfu, 1.0); \
      muhm_print_point((*x == s ? ">>s" : "  s"), s, fs, dfs, 1.0); \
      muhm_print_point((*x == v ? ">>v" : "  v"), v, fv, dfv, 1.0); \
      muhm_print_abscissa("  b", *b); \
      muhm_newline(); \
    }

  double s, fs, dfs;       /* Splitting point. */
  double xOld, aOld, bOld; /* Previous value of `x', `a', `b'. */
  double uOld, vOld;       /* Previous bracketing points `u', `v'. */
  double delta;            /* `tol' plus a machine-precision error */
  double reqRatio;         /* Maximum shrink ratio required at next shrinking step. */
  double ratio;         
  double step;
  affirm(*a < *x, "");
  affirm(*b > *x, "");

  /* Report initial interval and abort if good enough: */
  if (check!=NULL)
    { if (check(parms, *a, *b, *x, *fx, *dfx,1.0, 0.0,0.0)){ return; } }

  delta = tol + Eps*ABS(*x);

  affirm(delta > 0.0, "");

  /* Choose initial bracketing points `u', `v': */
  step = ABS(dist);
  u = MAX(*a, *x - step - delta);
  v = MIN(*b, *x + step + delta); 
  if (v - u < delta + delta){ return; }

  /* Evaluate `f' at `u,v', checking for termination: */
  if (muhm_eval_and_check(eval,check,parms, u,&fu,&dfu, x,fx,dfx, *a,*b))
    { return; }
  if (dfu <= 0.0){ *a = u; }

  if (muhm_eval_and_check(eval,check,parms, v,&fv,&dfv, x,fx,dfx, *a,*b))
    { return; }
  if (dfv >= 0.0){ *b = v; }
  
  xOld = *x; aOld = *a; bOld = *b;
  reqRatio = 1.0;

  /* Save original middle point: */
  s = *x; fs = *fx; dfs = *dfx;
  
  while (1)
    {
      /* Loop invariants:
          1. Probe points are definitely ordered:
            `v - u' >= delta'
          2. We are inside the current interval (which only shrinks):
            `*a <= u < v <= b'
          5. `s' is between `u' and `v': `u <= s <= v'
          6. `x' is the best point seen so far, or 
             the most recent in case of ties.
          7. `x' either `u', `v', or `s'.
      */

      affirm(*a <= u, "");
      affirm(u <= s, "");
      affirm(s <= v, "");
      affirm(v <= *b, "");
      affirm(u < v, "");
      affirm((*x == u) || (*x == s) || (*x == v), "");

      if (debug){ PrintStatus(); }

      delta = tol + Eps*ABS(*x);
      /* Check for termination criterion: */
      if (v - u < delta){ return; }

      uOld = u; vOld = v;
      if ((s != u) && (s != v))
        { /* We have a usable middle point `s'. */
          /* Reduce the interval `[u__v]' to `[u__s]' or `[s__v]': */
          bool_t shrink_u = FALSE; bool_t shrink_v = FALSE;
          if (debug){ muhm_message("three points"); }
          if (*x == u)
            { shrink_v = TRUE; }
          else if (*x == v)
            { shrink_u = TRUE; }
          else 
            { if (dfs >= 0) { shrink_v = TRUE; }
              if (dfs <= 0) { shrink_u = TRUE; }
            }
          if (shrink_u) 
            { if (debug){ muhm_message("increasing u"); }
              u = s; fu = fs; dfu = dfs; *a = s;
            }
          if (shrink_v) 
            { if (debug){ muhm_message("decreasing v"); }
              v = s; fv = fs; dfv = dfs; *b = s;
            }
          ratio = (v - u)/(vOld - uOld);
          if (debug){ muhm_print_abscissa("shrk", ratio); }
          reqRatio = ALPHA/MAX(ALPHA, ratio);
        }
      else if((dfu > 0) && (dfv < 0))
        { /* We don't have a useful middle point `s'. */
          /* There are local minima in `[a__u]' and `[v__b]'. */
          double h = MAX(v - u, delta);
          double mu = fu - dfu*MIN(h, u - *a);
          double mv = fv + dfv*MIN(h, *b - v);
          double t;
          if (debug){ muhm_message("three points, convex"); }
          if ((mu < mv) || ((mu == mv) && (*x == v)))
            { /* Extrapolates `u', and drags `v' unless it is the current best: */
              if (debug){ muhm_message("extrapolating down"); }
              if (*x == v)
                { u = MAX(*a, u - h);
                  if (u == uOld) { *a = *x; return; }
                }
              else
                { t = MAX(*a, u - 2*h); v = u; fv = fu; dfv = dfu; *b = u;  u = t;
                  if (u == uOld) { *b = *x; return; }
                }
              if (eval(parms, u, &fu, &dfu)) { return; }
              if (fu <= *fx) { *x = u; *fx = fu; *dfx = dfu; }
            }
          else if ((mu > mv) || ((mu == mv) && (*x == u)))
            { /* Extrapolates `v', and drags `u' unless it is the current best: */
              if (debug){ muhm_message("extrapolating up"); }
              if (*x == u)
                { v = MIN(*b, v + h);
                  if (v == vOld) { *b = *x; return; }
                }
              else
                { t = MIN(*b, v + 2*h); u = v; fu = fv; dfu = dfv; *a = v; v = t;
                  if (v == vOld) { *a = *x; return; }
                }
              if (eval(parms, v, &fv, &dfv)) { return; }
              if (fv <= *fx) { *x = v; *fx = fv; *dfx = dfv; }
            }
          else
            { affirm(FALSE, ""); }
          reqRatio = 1.0;
        }
      else
        { /* We don't have a usable middle point. */
          /* There is a local minimum in `[u__v]'. */
          if (debug){ muhm_message("two points, non-convex"); }
          /* Estimate its position `s'. */
          
          s = muhm_estimate_minimum(u, fu, dfu, v, fv, dfv);

          /* Adjust `s' so as to ensure sufficient progress: */
          if (debug){ muhm_print_abscissa("sraw", s); }
          if (debug){ muhm_print_abscissa("reqr", reqRatio); }
          s = muhm_ensure_progress(s, *a, u, v, *b, delta, reqRatio, debug);
          if (debug){ muhm_print_abscissa("sref", s); }

          /* Evaluate new point: */
          if (eval(parms, s, &fs, &dfs)){ return; }
          if (fs <= *fx){ *x = s; *fx = fs; *dfx = dfs; }
        }

      if ((*x != xOld) || (*a != aOld) || (*b != bOld))
        { /* Report current interval and point, and abort if good enough: */
          if (check != NULL)
            { if (check(parms, *a, *b, *x, *fx, *dfx,1.0, 0.0,0.0)){ return; } }
          xOld = *x; aOld = *a; bOld = *b;
        }
    }
  #undef PrintStatus
}

bool_t muhm_eval_and_check
  ( EvalProc eval, CheckProc check, void *parms,
    double t, double *ft, double *dft, 
    double *x, double *fx, double *dfx, 
    double a, double b
  )
{ if (t == (*x))
    { *ft = *fx; *dft = *dfx; }
  else
    { if (eval(parms, t, ft, dft))
        { return TRUE; }
      if (*ft < *fx)
        { *x = t; *fx = *ft; *dfx = *dft;
          if (check != NULL)
            { if (check(parms, a, b, *x, *fx, *dfx,1.0, 0.0,0.0)){ return TRUE; } }
        }
    }
  return FALSE;
}

double muhm_estimate_minimum
  ( double u, double fu, double dfu, 
    double v, double fv, double dfv
  )
{
  double w = v - u;
  /* Map interval `[u__v]' to `[-1__+1]': */
  double wdfu = w*dfu;
  double wdfv = w*dfv;
  /* Find Hermite interpolator  `A*x^3 +B*x^2 + C*x + D': */
  double B = (wdfv - wdfu)/2;
  double D = (fv + fu)/2 - B;
  double R = (fv - fu)/2;     /* A+C */
  double S = (wdfv + wdfu)/2; /* 3A + C */
  double A = (S-R)/2;
  double C = R - A;
  double x1, x2, x, s;
  /* Compute roots of quadratic `3*A*x^2 + 2*B*x + C': */
  muhm_quadratic_roots(A,B,C, &x1, &x2);
  /* Check which is the correct root: */
  x = (A > 0.0 ? x2 : x1);
  /* Clamp to interval: */
  if (x < -1) { x = -1; }
  if (x > +1) { x = +1; }
  /* Compare interpolated value to bracket point values: */
  if (((A*x + B)*x + C)*x + D > MIN(fu,fv))
    { x = (fu < fv ? -1 : +1); }
  /* Map back to `[u__v]': */
  s = u + (x + 1)/2*w;
  if (s < u) { s = u; }
  if (s > v) { s = v; }
  return s;
}

void muhm_quadratic_roots(double A, double B, double C, double *x1, double *x2)
{
  if (A != 0)
    { double delta;
      delta = B*B - 4*A*C;
      if (delta > 0)
        { double sd = sqrt(delta);
          double q = -0.5*(B + (B > 0 ? sd : -sd));
          *x1 = q/A;
          *x2 = C/q;
        }
      else
        { *x1 = -B/2/A;
          *x2 = *x1; 
        }
    }
  else if (B != 0)
    { if (B > 0)
        { *x1 = -INFINITY; *x2 = -C/B; }
      else
        { *x1 = -C/B; *x2 = INFINITY; }
    }
  else
    { *x1 = 0; *x2 = 0; }
}
  
double muhm_ensure_progress
  ( double s,
    double a, double u, double v, double b,
    double delta,
    double reqRatio,
    bool_t debug
  )
{
  double reqSep = (1 - reqRatio)*(v - u);
  double minSep = MAX(delta, reqSep);
  double sLo = u + minSep;
  double sHi = v + minSep;
  if (sLo > sHi)
    { return (0.5*u + 0.5*v); }
  else 
    { return (s < sLo ? sLo : (s > sHi ? sHi : s)); }
}

void muhm_message(char *msg)
{
  fprintf(stderr, "%s\n", msg);
}

void muhm_newline()
{
  fprintf(stderr, "\n");
}

void muhm_print_abscissa(char *msg, double u)
{ 
  fprintf(stderr, "    %s = %16g \n", msg, u);
}

void muhm_print_point(char *msg, double u, double fu, double dfu, double ddfu)
{ 
  fprintf(stderr, "    %s = %16g  F(%s) = %16g ", msg, u, msg, fu);
  if ((dfu != 0.0) || (ddfu != 0.0))
    { fprintf(stderr, " dF(%s) = %16g  ddF(%s) = %16g ", msg, dfu, msg, ddfu); }
  fprintf(stderr, "\n");
}

