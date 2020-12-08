// See minu_js.h
// Last edited on 2008-05-25 02:04:40 by stolfi

#define _GNU_SOURCE
#include <minu_gen.h>
#include <minu_js.h>

#include <affirm.h>

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define Eps (1.0e-6) /* Square root of relative machine precision. */
#define Phi  (1.61803398874989484821) /* Golden ratio */
#define Ihp (0.61803398874989484821) /* 1/Phi */

#define ALPHA Ihp

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

// INTERNAL PROTOTYPES

bool_t mujs_eval_and_check
  ( EvalProc eval, CheckProc check, void *parms,
    double t, double *ft, double *x, double *fx, double *dfx,
    double a, double b
  );
/*
  Calls `eval' on `t' with result returned in `*ft'.
  (unless `t==*x', in which case simply copies `*fx').
  If eval returns `TRUE', returns `TRUE'; else, if `*ft<=*fx',
  updates `*x,*fx,*dfx', calls `check' (with bogus second derivative),
  and returns whathever it returns. */
  
double mujs_ensure_progress
  ( double s,
    double a, double u, double m, double v, double b,
    double delta,
    double reqProgress,
    bool_t debug
  );
/* 
  Adjust the estimated minimum `s' to ensure ``sufficient progress''.
  Asking `reqProgress == 1' forces a golden search.
  Asking `reqProgress == 0' imposes only a minimum distance `delta'
  from the other probes. */

void mujs_message(char *msg);
void mujs_newline(void);
void mujs_print_abscissa(char *msg, double u);
void mujs_print_point(char *msg, double u, double fu, double dfu, double ddfu);

// IMPLEMENTATIONS

void minu_js_minimize
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
  double u, fu;    /* Lower bracketing point */
  double m, fm;    /* Middle point */
  double v, fv;    /* Upper bracketing point */
    
  #define PrintStatus() \
    {  \
      mujs_newline(); \
      mujs_print_abscissa("a", *a); \
      mujs_print_point("u", u, fu, 0, 0); \
      mujs_print_point("m", m, fm, 0, 0); \
      mujs_print_point("v", v, fv, 0, 0); \
      mujs_print_abscissa("b", *b); \
      mujs_newline(); \
    }

  double s, fs, dfs;   /* Splitting point. */
  double Dfx,Dq;       /* Estimated derivative at `x' is `Dfx/Dq'. */
  double DDfx,DDq;     /* Estimated second derivative at `x' is `DDfx/DDq'. */
  double xOld;         /* Previous value of `x' */
  double delta;        /* `tol' plus a machine-precision error */
  double reqProgress;  /* Required progress: 1 == strict golden, 0 == any. */
  double progress;     /* Actual progress: 1 == golden, 1/Phi == none, oo == wow. */  
  double step, back;
  affirm(*a < *x, "");
  affirm(*b > *x, "");

  /* Report initial interval and abort if good enough: */
  if (check!=NULL)
    { if (check(parms, *a, *b, *x, *fx, 0.0,0.0, 0.0,0.0)){ return; } }

  delta = tol + Eps*ABS(*x);

  affirm(delta > 0.0, "");

  /* Choose initial bracketing points `u', `v': */
  step = ABS(dist);
  back = step/Phi;
  if (dist > 0.0)
    { v = MIN(*b, *x + step + delta); 
      u = MAX(*a, *x - back - delta);
    }
  else
    { u = MAX(*a, *x - step - delta);
      v = MIN(*b, *x + back + delta); 
    }
  if (v - u < delta + delta){ return; }

  /* Evaluate `f' at `u,v', checking for termination: */
  if (mujs_eval_and_check(eval, check, parms, u, &fu, x, fx, dfx, *a, *b)){ return; }
  if (mujs_eval_and_check(eval, check, parms, v, &fv, x, fx, dfx, *a, *b)){ return; }

  /* Pick middle point: */
  if ((*x - u < delta) || (v - *x < delta))
    { m = 0.5*(u + v);
      affirm(m!=*x, "");
      if (eval(parms, m, &fm, &dfs)){ return; }
    }
  else
    { m = *x; fm = *fx; }

  xOld = *x;
  reqProgress = delta/(v - u);

  while (1)
    {
      affirm(*a <= u, "");
      affirm(u < m, "");
      affirm(m < v, "");
      affirm(v <= *b, "");

      if (fu >= fm){ *a = u; }
      if (fv >= fm){ *b = v; }

      delta = tol + Eps*ABS(*x);

      /* Loop invariants:
          1. Probe points are definitely ordered:
            `u + delta <= m <= v - delta'
          2. We are inside the current interval (which only shrinks):
            `u >= a, v <= b'
          4. `x,fx,dfx' is the the best point found so far,
            (the most recent, if there were ties). 
          5. `x' is between `u' and `v':
             `u <= x <= v'
      */

      if (debug){ PrintStatus(); }

      /* Check for termination criterion: */
      if ((*x - delta <= *a) && (*x + delta >= *b)){ return; }

      /* Estimate derivatives at `x': */
      if (ABS(*x - m) < delta)
        { ComputeDerivatives (u, fu, v, fv, *x, *fx, &Dfx,&Dq, &DDfx,&DDq); }
      else if ((ABS((*x) - v) < delta) || (((*x) > m) && (ABS((*x) - u) >= delta)))
        { ComputeDerivatives (u, fu, m, fm, *x, *fx, &Dfx,&Dq, &DDfx,&DDq); }
      else
        { ComputeDerivatives (v, fv, m, fm, *x, *fx, &Dfx,&Dq, &DDfx,&DDq); }
      affirm(Dq != 0.0, "");
      affirm(DDq != 0.0, "");

      if (*x != xOld)
        { /* Report current interval and point, and abort if good enough: */
          if (check != NULL)
            { if (check(parms, *a, *b, *x, *fx, Dfx,Dq, DDfx,DDq)){ return; } }
          xOld = *x;
        }

      /* Estimate minimum position by parabolic fit: */
      if (DDfx <= 0.0)
        {
          if (debug){ mujs_message("convex triple"); }
          if (Dfx < 0.0)
            { s = INFINITY; }
          else if (Dfx > 0.0)
            { s = -INFINITY; }
          else
            { s = *x; }
        }
      else
        { s = *x - (Dfx/DDfx)*(DDq/Dq); }

      if (debug){ mujs_print_abscissa("sraw == ", s); }
      if (debug){ mujs_print_abscissa("rqPr == ", reqProgress); }
      s = mujs_ensure_progress(s, *a, u, m, v, *b, delta, reqProgress, debug);
      if (debug){ mujs_print_abscissa("sref == ", s); }

      /* Evaluate new point: */
      if (eval(parms, s, &fs, &dfs)){ return; }
      if (fs < *fx){ *x = s; *fx = fs; *dfx = dfs; }

      /* Discard one point and rename others: */
      if (s < u)
        { double t = (m-s)/(u-v)/Phi;
          progress = MAX(t, 1.0/t);
          if (debug){ mujs_message("extrapolating down"); }
          v = m; fv = fm;
          m = u; fm = fu;
          u = s; fu = fs;
        }
      else if (s > v)
        { double t = (s-m)/(u-v)/Phi;
          progress = MAX(t, 1.0/t);
          if (debug){ mujs_message("extrapolating up"); }
          u = m; fu = fm;
          m = v; fm = fv;
          v = s; fv = fs;
        }
      else
        { double uvOld = v - u;
          if (fs < fm)
            {
              if (s <= m)
                { if (debug){ mujs_message("new minimum in [u__m]"); }
                  v = m; fv = fm;
                }
              else
                { if (debug){ mujs_message("new minimum in [m__v]"); }
                  u = m; fu = fm;
                }
              m = s; fm = fs;
            }
          else
            { if (s <= m)
                { if (debug){ mujs_message("high point in [u__m]"); }
                  u = s; fu = fs;
                }
              else
                { if (debug){ mujs_message("high point in [m__v]"); }
                  v = s; fv = fs;
                }
            }
          progress = uvOld/(v-u)/Phi;
        }

      /* Compute min progress required in next step: */
      if (debug){ mujs_print_abscissa("prpg == ", progress); }
      reqProgress = MIN(1.0, MAX(delta/(v - u), ALPHA * reqProgress/progress));
    }
  #undef PrintStatus
}

bool_t mujs_eval_and_check
  ( EvalProc eval, CheckProc check, void *parms,
    double t, double *ft, double *x, double *fx, double *dfx, 
    double a, double b
  )
{ double dft;
  if (t == (*x))
    { *ft = *fx; }
  else
    { if (eval(parms, t, ft, &dft))
        { return TRUE; }
      if (*ft < *fx)
        { *x = t; *fx = *ft; *dfx = dft;
          if (check != NULL)
            { if (check(parms, a, b, *x, *fx, 0.0,0.0, 0.0,0.0)){ return TRUE; } }
        }
    }
  return FALSE;
}
  
double mujs_ensure_progress
  ( double s,
    double a, double u, double m, double v, double b,
    double delta,
    double reqProgress,
    bool_t debug
  )
{
  /* 
    We first compute intervals 
    | 
    |   `[auLo _ auHi]' contained in `[a _ u]'
    |   `[umLo _ umHi]' contained in `[u _ m]'
    |   `[mvLo _ mvHi]' contained in `[m _ v]'
    |   `[vbLo _ vbHi]' contained in `[v _ b]'
    | 
    such that placing the next probe `s' in one of those intervals 
    will ensure ``sufficient progress'' in the next step, whatever
    the value of `f(s)'.

    Currently, ``sufficient progress'' means that the probe `s' is
    at least `delta' away from the current probes, and the next
    interval `[u _ v]' will be of ``adequate size''.

    For extrapolation steps, with `s < u' or `s > v', ``adequate
    size'' means that the extrapolation ratio `t' (`(u-s)/(m-u)' or
    `(s-v)/(v-m)', respectively) is between `Phi*reqProgress' and
    `Phi/reqProgress'.  The `progress' made in such a step is by
    definition the ratio `MIN(t/Phi,Phi/t)'.

    For interpolation steps, with `s' in `(u _ m)' or or in `(m _
    v)', ``adequate size'' means the worst-case shinking factor `t'
    (`MAX(m-u, v-s)/(v-u)' or `MAX(v-m, s-u)/(v-u)', respectively)
    is at most `1/(reqProgress*Phi)'.  The `progress' made in such a
    step is by definition `1/(Phi*t')', where `t'' is the shrinking
    ratio actually attained in the step. 

    Note that `progress' may be less than `reqProgress' if the 
    distances `m-u' and `v-m' are too dissimilar, or 
    we are too close to the endpoints of `[a _ b]'.
    But this condition will hopefully be corrected in the 
    following steps.

    In either case, the next value of `reqProgress' 
    should be  `alpha * reqProgress/progress', where `alpha' is some
    ``venture coefficient'' between 0 and 1. */

  double phiDnStep = Phi*(m - u);
  double phiNuWidth = (v - u)/Phi;
  double phiUpStep = Phi*(v - m);
  double maxDnStep = phiDnStep/reqProgress;
  double minDnStep = MAX(delta, phiDnStep*reqProgress);
  double maxNuWidth = phiNuWidth/reqProgress;
  double minUpStep = MAX(delta, phiUpStep*reqProgress);
  double maxUpStep = phiUpStep/reqProgress;

  double auLo = MAX(a, u - maxDnStep);
  double auHi = u - minDnStep;
  double umLo = MAX(u + delta, v - maxNuWidth);
  double umHi = MIN(m - delta, u + maxNuWidth);
  double mvLo = MAX(m + delta, v - maxNuWidth);
  double mvHi = MIN(v - delta, u + maxNuWidth);
  double vbLo = v + minUpStep;
  double vbHi = MIN(b, v + maxUpStep);

  if ((umLo > umHi) && (mvLo > mvHi))
    { mujs_message("no space to split!");
      mujs_print_abscissa("umLo == ", umLo);
      mujs_print_abscissa("umHi == ", umHi);
      mujs_print_abscissa("mvLo == ", mvLo);
      mujs_print_abscissa("mvHi == ", mvHi);
    }
  affirm((umLo <= umHi) || (mvLo <= mvHi), "");
  /* Adjust estimate to preserve loop invariants: */
  if (s < u)
    {
      /* Try to probe before `u' */
      if (auLo <= auHi)
        { s = MAX(auLo, MIN(auHi, s)); }
      else
        { if (debug){ mujs_message("can't split [a_u], splitting [u_v]"); }
          s = (umLo <= umHi ? umLo : mvLo);
        }
    }
  else if (s > v)
    {
      /* Try to probe after `v' */
      if (vbLo <= vbHi)
        { s = MIN(vbHi, MAX(vbLo, s)); }
      else
        { if (debug){ mujs_message("can't split [v_b], splitting [u_v]"); }
          s = (mvLo <= mvHi ? mvHi : umHi); 
        }
    }
  else if (s <= m)
    { if (umLo <= umHi)
        { s = MAX(umLo, MIN(umHi, s)); }
      else
        { if (debug){ mujs_message("can't split [u_m], splitting [m_v]"); }
          s = mvLo;
        }
    }
  else
    { if (mvLo <= mvHi)
        { s = MAX(mvLo, MIN(mvHi, s)); }
      else
        { if (debug){ mujs_message("can't split [m_v], splitting [u_m]"); }
          s = umHi;
        }
    }
  return s;
}

void mujs_message(char *msg)
{
  fprintf(stderr, "%s\n", msg);
}

void mujs_newline()
{
  fprintf(stderr, "\n");
}

void mujs_print_abscissa(char *msg, double u)
{ 
  fprintf(stderr, "    %s = %20g\n", msg, u);
}

void mujs_print_point(char *msg, double u, double fu, double dfu, double ddfu)
{ 
  fprintf(stderr, "    %s = %20g F(%s) = %20g", msg, u, msg, fu);
  if ((dfu != 0.0) || (ddfu != 0.0))
    { fprintf(stderr, " dF(%s) = %20g ddF(%s) = %20g", msg, dfu, msg, ddfu); }
  fprintf(stderr, "\n");
}

