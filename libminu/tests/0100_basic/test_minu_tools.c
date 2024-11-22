// See test_minu_tools.h
// Last edited on 2024-11-08 16:54:59 by stolfi

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <jsfile.h>
#include <jsrandom.h>
#include <jsstring.h>
#include <epswr.h>
#include <affirm.h>

#include <minu_gen.h>

#include <test_minu_tools.h>

typedef struct TestParms {
  Problem *prb;         /* Problem parameters */
  Minimizer *opt;       /* Optimizer parameters */
  epswr_figure_t *eps;  /* Encapsulated Postscript file */
  unsigned nCalls;      /* Number of function calls */
} TestParms;

// PROTOTYPES

epswr_figure_t *test_minu_tools_new_figure
  ( int32_t i_opt, 
    Minimizer *opt, 
    int32_t i_prb, 
    Problem *prb, 
    char *tag
  );
  /* Creates a new Encapsulated Postscript stream {eps} for a test of
    minimizer {opt} on problem {prb}. The figure will be written to the
    file "out/aumt_opt{OO}_prb{PP}_{tag}.eps" where {OO} is {i_opt} and
    {PP} is {i_prb}, both formatted as "%02d". Also plots the graph of
    the problem's function, in the range specified in {prb}. */

void test_minu_tools_end_figure
  ( epswr_figure_t *eps, 
    unsigned nCalls, 
    double error, 
    Minimizer *opt, 
    Problem *prb
  );
  /* Writes captions to figure {eps} and closes the file. */
  
void test_minu_tools_draw_segment
  ( epswr_figure_t *eps, 
    Problem *prb, 
    double xa, 
    double ya, 
    double xb, 
    double yb
  );

void test_minu_tools_plot_dot
  ( epswr_figure_t *eps, 
    Problem *prb, 
    double x, 
    double y, 
    bool_t black
  );

void test_minu_tools_print_status 
  ( unsigned nCalls,
    double err, 
    double x, 
    double fx,
    double dfx, 
    double dq, 
    double ddfx, 
    double ddq
  );
    
bool_t test_minu_tools_single_eval
  ( void *parms, /* Actually a TestParms */
    double x, 
    double *fx, 
    double *dfx
  );

bool_t test_minu_tools_single_check
  ( void *parms, /* Actually a TestParms */ 
    double a, double b, 
    double x, double fx,
    double dfx, double dq, double ddfx, double ddq
  );

bool_t test_minu_tools_multiple_eval
  ( void *parms, /* Actually a TestParms */ 
    double x, double *fx, double *dfx
  );

bool_t test_minu_tools_multiple_check 
  ( void *parms, /* Actually a TestParms */ 
    double a, double b,
    double x, double fx,
    double dfx, double dq, double ddfx, double ddq
  );

// IMPLEMENTATIONS

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

Performance test_minu_tools_single
  ( int32_t i_opt,   /* Index of minimizer, for file name. */
    Minimizer *opt,  /* The minimization tool */
    int32_t i_prb,   /* Index of problem, for file name. */
    Problem *prb,    /* Function and parameters */
    bool_t debug     /* Passed to the minimizer */
  )
  { epswr_figure_t *eps = test_minu_tools_new_figure(i_opt, opt, i_prb, prb, "sngl");
    TestParms tp = (TestParms){prb, opt, eps, 0};
    fprintf(stderr, "\n");
    fprintf(stderr, "=============================================================================\n");
    fprintf(stderr, "===== test_minu_tools_single(%s, %s)\n\n", opt->name, prb->name);

    double x = prb->xStart;
    double fx, dfx;
    prb->eval(prb, x, &fx, &dfx);
    tp.nCalls = 1;
    double a = prb->xMin;
    double b = prb->xMax;
    opt->minimize(
      /*parms*/ &tp,
      /*eval*/  test_minu_tools_single_eval,
      /*x*/     &x,
      /*fx*/    &fx,
      /*dfx*/   &dfx,
      /*tol*/   prb->tol,
      /*dist*/  prb->dist,
      /*a*/     &a,
      /*b*/     &b,
      /*check*/ test_minu_tools_single_check,
      /*debug*/ debug
    );
    double error = prb->error(prb, x, fx);
    
    fprintf(stderr, "\n");
    fprintf(stderr, "===== Final Result:\n");
    test_minu_tools_print_status(tp.nCalls, error, x, fx, dfx,1.0, 0.0,0.0);

    fprintf(stderr, "=============================================================================\n");
    fprintf(stderr, "\n");
    test_minu_tools_end_figure(eps, tp.nCalls, error, opt, prb);
    
    return (Performance){
        /*nTests*/    1,
        /*avgCalls*/  (double)tp.nCalls,
        /*devCalls*/  0.0,
        /*maxCalls*/  tp.nCalls,
        /*avgError*/  MAX(0.0, error),
        /*devError*/  0.0,
        /*maxError*/  MAX(0.0, error),
        /*nFailures*/ (error > 0.0 ? 1 : 0)
     };
  }

Performance test_minu_tools_multiple 
  ( int32_t i_opt,     /* Index of minimizer, for file name. */
    Minimizer *opt,    /* The minimization tool */
    int32_t i_prb,     /* Index of problem, for file name. */
    Problem *prb,      /* Function and parameters */
    unsigned nTests    /* Number of tests to perform. */
  )
  {
    double totCalls, totCallsSq = 0.0;
    double totError, totErrorSq = 0.0;
    unsigned maxCalls = 0;
    double maxError = 0.0;
    unsigned nFailures = 0;
    TestParms tp = (TestParms){prb, opt, NULL, 0};
          
    fprintf(stderr, "\n");
    fprintf(stderr, "=============================================================================\n");
    fprintf(stderr, "===== test_minu_tools_multiple (%s,%s)\n\n", opt->name, prb->name);
    affirm(nTests >= 2, "");
 
    for (int32_t it = 1; it <= nTests; it++)
      { double x = prb->xStart + prb->dist * (2*drandom() - 1);
        double fx, dfx;
        prb->eval(prb, x, &fx, &dfx);
        tp.nCalls = 1;
        double a = prb->xMin;
        double b = prb->xMax;
        opt->minimize(
          /*parms*/ &tp,
          /*eval*/  test_minu_tools_multiple_eval,
          /*x*/     &x,
          /*fx*/    &fx,
          /*dfx*/   &dfx,
          /*tol*/   prb->tol,
          /*dist*/  prb->dist,
          /*a*/     &a,
          /*b*/     &b,
          /*check*/ test_minu_tools_multiple_check,
          /*debug*/ FALSE
        );
        double nc = (double)tp.nCalls;
        totCalls = totCalls + nc;
        totCallsSq = totCallsSq + nc * nc;
        maxCalls = MAX(maxCalls, tp.nCalls);
        
        double ee = prb->error(prb, x, fx);
        if (ee > 0.0)
          { 
            totError = totError + ee;
            totErrorSq = totErrorSq + ee * ee;
            maxError = MAX(maxError, ee);
            nFailures++;
          }
      }

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

bool_t test_minu_tools_single_eval(void *parms, double x, double *fx, double *dfx)
  {
    TestParms *tp = (TestParms *)parms;
    Problem *prb = tp->prb;
    prb->eval(prb, x, fx, dfx);
    test_minu_tools_plot_dot(tp->eps, prb,
      MIN(prb->xMax, MAX(prb->xMin, x)), 
      MIN(prb->yMax, MAX(prb->yMin, *fx)),
      FALSE
    );
    tp->nCalls++;
    return (tp->nCalls >= prb->maxCalls);
  }

bool_t test_minu_tools_single_check
  ( void *parms, 
    double a, double b, 
    double x, double fx,
    double dfx, double dq, double ddfx, double ddq
  )
  {
    TestParms *tp = (TestParms *)parms;
    Problem *prb = tp->prb;
    double err = prb->error(prb, x, fx);
    test_minu_tools_print_status(tp->nCalls, err, x, fx, dfx,dq, ddfx,ddq);
    test_minu_tools_plot_dot(tp->eps, prb, x, fx, TRUE);
    return ((err <= 0.0) || (tp->nCalls >= prb->maxCalls));
  }

bool_t test_minu_tools_multiple_eval(void *parms, double x, double *fx, double *dfx)
  {
    TestParms *tp = (TestParms *)parms;
    Problem *prb = tp->prb;
    prb->eval(prb, x, fx, dfx);
    tp->nCalls++;
    return (tp->nCalls >= prb->maxCalls);
  }

bool_t test_minu_tools_multiple_check 
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

epswr_figure_t *test_minu_tools_new_figure
  ( int32_t i_opt, 
    Minimizer *opt, 
    int32_t i_prb, 
    Problem *prb, 
    char *tag
  )
  { char *fname = NULL;
    char *fname = jsprintf("aumt_opt%02d_prb%02d_%s", i_opt, i_prb, tag);
    double fontHeight = 10;
    double hSize = 425; /* Plot width (pt). */
    double vSize = 425; /* Plot height (pt). */
    double mrg = 4.0;   /* Figure margin (pt). */
    int32_t nCap = 8;
    double mrg_bot = mrg + nCap*fontHeight + mrg;
    bool_t eps_verbose = FALSE;
    epswr_figure_t *eps = epswr_new_named_figure
        ( "out", fname, NULL, -1, NULL,
          hSize, vSize, mrg, mrg, mrg_bot, mrg, 
          eps_verbose
        );
    epswr_set_text_geometry(eps, FALSE, 0,hSize, mrg-mrg_bot, -mrg, 0.0);
    epswr_set_text_font(eps, "Courier", fontHeight);
        
    epswr_set_client_window(eps, 0.0,1.0, 0.0,1.0);

    int32_t NSteps = 1000;
    double dx = (prb->xMax - prb->xMin) / ((double)NSteps);
    double xb = prb->xMin;
    double yb, dfxb;
    prb->eval(prb, xb, &yb, &dfxb);
    epswr_set_pen(eps, 1.0,0.0,0.5, 0.10, 0.0, 0.0);
    for (int32_t i = 1; i <= NSteps; i++)
      { double xa = xb; 
        double ya = yb;
        xb = xb + dx; 
        prb->eval(prb, xb, &yb, &dfxb);
        test_minu_tools_draw_segment(eps, prb, xa, ya, xb, yb);
      }
    free(fname);
    return eps;
  }

void test_minu_tools_draw_segment
  ( epswr_figure_t *eps,
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
    epswr_segment(eps, xpa, ypa, xpb, ypb);
  }

void test_minu_tools_plot_dot
  ( epswr_figure_t *eps,
    Problem *prb,
    double x,
    double y,
    bool_t black
  )
  { double xf = 1.0/(prb->xMax - prb->xMin);
    double yf = 1.0/(prb->yMax - prb->yMin);
    double xp = xf*(x - prb->xMin);
    double yp = yf*(y - prb->yMin);
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.10, 0.0, 0.0);
    if (black)
      { epswr_set_fill_color(eps, 0.0,0.0,0.0); }
    else
      { epswr_set_fill_color(eps, 1.0,1.0,1.0); }
    epswr_dot(eps, xp, yp, 0.6, TRUE, TRUE);
  }

void test_minu_tools_end_figure
  ( epswr_figure_t *eps, 
    unsigned nCalls,
    double error,
    Minimizer *opt,
    Problem *prb
  )
  {
    char buf[1200];
    
    epswr_set_fill_color(eps, 0,0,0);
    
    epswr_text(eps, txtcat("Minimizer: ", opt->name), FALSE, 0.0, TRUE, FALSE);
    
    epswr_text(eps, txtcat("Problem: ", prb->name), FALSE, 0.0, TRUE, FALSE);
    
    sprintf(buf, "  tol = %f", prb->tol);   
    epswr_text(eps, buf, FALSE, 0.0, TRUE, FALSE);
    
    sprintf(buf, "  dist = %f", prb->dist);  
    epswr_text(eps, buf, FALSE, 0.0, TRUE, FALSE);
    
    sprintf(buf, "  xRange = [%f __ %f]", prb->xMin, prb->xMax); 
    epswr_text(eps, buf, FALSE, 0.0, TRUE, FALSE);
    
    sprintf(buf, "  yRange = [%f __ %f]", prb->yMin, prb->yMax); 
    epswr_text(eps, buf, FALSE, 0.0, TRUE, FALSE);
    
    sprintf(buf, "  nCalls = %d", nCalls);  
    epswr_text(eps, buf, FALSE, 0.0, TRUE, FALSE);
    
    sprintf(buf, "  error = %f", error);  
    epswr_text(eps, buf, FALSE, 0.0, TRUE, FALSE);
    
    epswr_end_figure(eps);
  }


void test_minu_tools_print_status 
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
