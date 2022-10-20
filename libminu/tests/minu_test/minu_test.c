// See minu_test.h
// Last edited on 2022-10-20 06:28:27 by stolfi

#include <minu_gen.h>
#include <stdint.h>
#include <minu_test.h>

#include <affirm.h>
#include <pswr.h>
#include <jsstring.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <pswr.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct TestParms {
  Problem *prb;      /* Problem parameters */
  Minimizer *opt;    /* Optimizer parameters */
  PSStream *ps;         /* Postscript file */
  unsigned nCalls;   /* Number of function calls */
} TestParms;

// PROTOTYPES

void minu_test_init_plot(PSStream *ps, int32_t page, Problem *prb);
void minu_test_end_plot
  ( PSStream *ps, unsigned nCalls, double error, Minimizer *opt, Problem *prb );
void minu_test_draw_segment
  ( PSStream *ps, Problem *prb, double xa, double ya, double xb, double yb );
void minu_test_plot_dot(PSStream *ps, Problem *prb, double x, double y, bool_t black);
void minu_test_print_status 
  ( unsigned nCalls, double err, 
    double x, double fx, double dfx, double dq, double ddfx, double ddq
  );
    
bool_t minu_test_single_eval
  ( void *parms, /* Actually a TestParms */
    double x, double *fx, double *dfx
  );

bool_t minu_test_single_check
  ( void *parms, /* Actually a TestParms */ 
    double a, double b, 
    double x, double fx,
    double dfx, double dq, double ddfx, double ddq
  );

bool_t minu_test_multiple_eval
  ( void *parms, /* Actually a TestParms */ 
    double x, double *fx, double *dfx
  );

bool_t minu_test_multiple_check 
  ( void *parms, /* Actually a TestParms */ 
    double a, double b,
    double x, double fx,
    double dfx, double dq, double ddfx, double ddq
  );

// IMPLEMENTATIONS

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

Performance minu_test_single
  ( PSStream *ps,         /* Postscript file */\
    int32_t page,          /* Page number in document */
    Minimizer *opt,    /* The minimization tool */
    Problem *prb,      /* Function and parameters */
    bool_t debug         /* Passed to the minimizer */
  )
  {
    double x, fx, dfx;
    double a = prb->xMin;
    double b = prb->xMax;

    double error;
    TestParms test_parms = (TestParms){prb, opt, ps, 0};
    TestParms *tp = &test_parms;

    fprintf(stderr, "\n");
    fprintf(stderr, "=============================================================================\n");
    fprintf(stderr, "===== minu_test_single(%s, %s)\n\n", opt->name, prb->name);

    minu_test_init_plot(ps, page, prb);
    x = prb->xStart;
    prb->eval(prb, x, &fx, &dfx);
    tp->nCalls = 1;
    opt->minimize(
      /*parms*/ tp,
      /*eval*/  minu_test_single_eval,
      /*x*/     &x,
      /*fx*/    &fx,
      /*dfx*/   &dfx,
      /*tol*/   prb->tol,
      /*dist*/  prb->dist,
      /*a*/     &a,
      /*b*/     &b,
      /*check*/ minu_test_single_check,
      /*debug*/ debug
    );
    error = prb->error(prb, x, fx);
    fprintf(stderr, "\n");
    fprintf(stderr, "===== Final Result:\n");
    minu_test_print_status(tp->nCalls, error, x, fx, dfx,1.0, 0.0,0.0);

    fprintf(stderr, "=============================================================================\n");
    fprintf(stderr, "\n");
    minu_test_end_plot(ps, tp->nCalls, error, opt, prb);
    return (Performance){
        /*nTests*/    1,
        /*avgCalls*/  (double)tp->nCalls,
        /*devCalls*/  0.0,
        /*maxCalls*/  tp->nCalls,
        /*avgError*/  MAX(0.0, error),
        /*devError*/  0.0,
        /*maxError*/  MAX(0.0, error),
        /*nFailures*/ (error > 0.0 ? 1 : 0)
     };
  }

void minu_test_init_plot(PSStream *ps, int32_t page, Problem *prb)
  {
    double xa, ya, xb, yb, dfxb;
    int32_t NSteps = 400;
    double dx = (prb->xMax - prb->xMin) / ((double)NSteps);
    int32_t i;

    pswr_new_canvas(ps, NULL);
    pswr_set_window(ps, 
      0.0,   1.0,   0.0,   1.0,
      95.0,  520.0, 275.0, 700.0
    );
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    xb = prb->xMin;
    prb->eval(prb, xb, &yb, &dfxb);
    for (i = 1; i <= NSteps; i++)
      { xa = xb; ya = yb;
        xb = xb + dx; prb->eval(prb, xb, &yb, &dfxb);
        minu_test_draw_segment(ps, prb, xa, ya, xb, yb);
      }
  }

void minu_test_draw_segment
  ( PSStream *ps,
    Problem *prb,
    double xa,
    double ya,
    double xb,
    double yb
  )
  {
    double xf = 1.0/(prb->xMax - prb->xMin);
    double yf = 1.0/(prb->yMax - prb->yMin);
    double xpa = xf*(xa - prb->xMin);
    double ypa = yf*(ya - prb->yMin);
    double xpb = xf*(xb - prb->xMin);
    double ypb = yf*(yb - prb->yMin);
    pswr_segment(ps, xpa, ypa, xpb, ypb);
  }

void minu_test_plot_dot
  ( PSStream *ps,
    Problem *prb,
    double x,
    double y,
    bool_t black
  )
  {
    double xf = 1.0/(prb->xMax - prb->xMin);
    double yf = 1.0/(prb->yMax - prb->yMin);
    double xp = xf*(x - prb->xMin);
    double yp = yf*(y - prb->yMin);
    if (black)
      { pswr_set_fill_color(ps, 0.0,0.0,0.0); }
    else
      { pswr_set_fill_color(ps, 1.0,1.0,1.0); }
    pswr_dot(ps, xp, yp, 0.4, TRUE, TRUE);
  }

void minu_test_end_plot
  ( PSStream *ps, 
    unsigned nCalls,
    double error,
    Minimizer *opt,
    Problem *prb
  )
  {
    char buf[120];
    pswr_frame(ps);
    pswr_add_caption(ps, txtcat("Minimizer: ", opt->name), 0.0);
    pswr_add_caption(ps, txtcat("Problem: ", prb->name), 0.0);
    sprintf(buf, "  tol = %f", prb->tol);    pswr_add_caption(ps, buf, 0.0);
    sprintf(buf, "  dist = %f", prb->dist);  pswr_add_caption(ps, buf, 0.0);
    sprintf(buf, "  xRange = [%f __ %f]", prb->xMin, prb->xMax); pswr_add_caption(ps, buf, 0.0);
    sprintf(buf, "  yRange = [%f __ %f]", prb->yMin, prb->yMax); pswr_add_caption(ps, buf, 0.0);
    sprintf(buf, "  nCalls = %d", nCalls);  pswr_add_caption(ps, buf, 0.0);
    sprintf(buf, "  error = %f", error);  pswr_add_caption(ps, buf, 0.0);
  }

void minu_test_print_status 
  ( unsigned nCalls,
    double err, 
    double x, double fx, double dfx,
    double dq, double ddfx, double ddq
  )
  {
    fprintf(stderr, "*");
    fprintf(stderr, " calls = %6d error = %9.2f\n", nCalls, err);
    fprintf(stderr, " ");
    fprintf(stderr, " x = %12f fx = %12f", x, fx);
    if (dq != 0.0)
      {
        fprintf(stderr, " dfx = %12f", dfx);
        if (dq != 1.0)
          { fprintf(stderr, "/%-12f", dq); }
      }
    if (ddq != 0.0)
      {
        fprintf(stderr, " ddfx = %12f", ddfx);
        if (ddq != 1.0)
          { fprintf(stderr, "/%-12f", ddq); }
      }
    fprintf(stderr, "\n");
  }

bool_t minu_test_single_eval(void *parms, double x, double *fx, double *dfx)
  {
    TestParms *tp = (TestParms *)parms;
    Problem *prb = tp->prb;
    prb->eval(prb, x, fx, dfx);
    minu_test_plot_dot(tp->ps, prb,
      MIN(prb->xMax, MAX(prb->xMin, x)), 
      MIN(prb->yMax, MAX(prb->yMin, *fx)),
      FALSE
    );
    tp->nCalls++;
    return (tp->nCalls >= prb->maxCalls);
  }

bool_t minu_test_single_check
  ( void *parms, 
    double a, double b, 
    double x, double fx,
    double dfx, double dq, double ddfx, double ddq
  )
  {
    TestParms *tp = (TestParms *)parms;
    Problem *prb = tp->prb;
    double err = prb->error(prb, x, fx);
    minu_test_print_status(tp->nCalls, err, x, fx, dfx,dq, ddfx,ddq);
    minu_test_plot_dot(tp->ps, prb, x, fx, TRUE);
    return ((err <= 0.0) || (tp->nCalls >= prb->maxCalls));
  }

Performance minu_test_multiple 
  ( Minimizer *opt,    /* The minimization tool */
    Problem *prb,      /* Function and parameters */
    unsigned nTests    /* Number of tests to perform. */
  )
  {
    double totCalls, totCallsSq = 0.0;
    double totError, totErrorSq = 0.0;
    unsigned maxCalls = 0;
    double maxError = 0.0;
    unsigned nFailures = 0;
    TestParms test_parms = (TestParms){prb, opt, NULL, 0};
    TestParms *tp = &test_parms;
    int32_t i;

    fprintf(stderr, "\n");
    fprintf(stderr, "=============================================================================\n");
    fprintf(stderr, "===== minu_test_multiple (%s,%s)\n\n", opt->name, prb->name);
    affirm(nTests >= 2, "");

    for (i = 1; i <= nTests; i++)
      {
        double x, fx, dfx;
        double a = prb->xMin;
        double b = prb->xMax;

        x = prb->xStart + prb->dist * (2*drandom() - 1);
        prb->eval(prb, x, &fx, &dfx);
        tp->nCalls = 1;
        opt->minimize(
          /*parms*/ tp,
          /*eval*/  minu_test_multiple_eval,
          /*x*/     &x,
          /*fx*/    &fx,
          /*dfx*/   &dfx,
          /*tol*/   prb->tol,
          /*dist*/  prb->dist,
          /*a*/     &a,
          /*b*/     &b,
          /*check*/ minu_test_multiple_check,
          /*debug*/ FALSE
        );
        { double nc = (double)tp->nCalls;
          totCalls = totCalls + nc;
          totCallsSq = totCallsSq + nc * nc;
          maxCalls = MAX(maxCalls, tp->nCalls);
        }
        { double ee = prb->error(prb, x, fx);
          if (ee > 0.0)
            { 
              totError = totError + ee;
              totErrorSq = totErrorSq + ee * ee;
              maxError = MAX(maxError, ee);
              nFailures++;
            }
        }
      }
    {
      double nt = (double)nTests;
      double avgCalls = totCalls / nt;
      double varCalls = (totCallsSq - avgCalls * totCalls) / (nt - 1.0);
      double devCalls = sqrt(varCalls);

      double avgError = totError / nt;
      double varError = (totErrorSq - avgError * totError) / (nt - 1.0);
      double devError = sqrt(varError);

      fprintf(stderr, "\n");
      fprintf(stderr, "=== Calls:  avg = %9.2f  dev = %9.2f  max = %9d\n", 
        avgCalls, devCalls, maxCalls);
      fprintf(stderr, "=== Error:  avg = %9.2f  dev = %9.2f  max =  %9.2f  failures = %6d\n",
        avgError, devError, maxError, nFailures);
      fprintf(stderr, "=============================================================================\n");
      fprintf(stderr, "\n");
      return 
        (Performance){
          /*nTests*/    nTests,
          /*avgCalls*/  avgCalls,
          /*devCalls*/  devCalls,
          /*maxCalls*/  maxCalls,
          /*avgError*/  avgError,
          /*devError*/  devError,
          /*maxError*/  maxError,
          /*nFailures*/ nFailures
       };
    }
  }

bool_t minu_test_multiple_eval(void *parms, double x, double *fx, double *dfx)
  {
    TestParms *tp = (TestParms *)parms;
    Problem *prb = tp->prb;
    prb->eval(prb, x, fx, dfx);
    tp->nCalls++;
    return (tp->nCalls >= prb->maxCalls);
  }

bool_t minu_test_multiple_check 
  ( void *parms, 
    double a, double b,
    double x, double fx,
    double dfx, double dq, double ddfx, double ddq
  )
  { 
    TestParms *tp = (TestParms *)parms;
    Problem *prb = tp->prb;
    double err = prb->error(prb, x, fx);
    return ((err <= 0.0) || (tp->nCalls >= prb->maxCalls));
  }
