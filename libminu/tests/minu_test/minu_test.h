// Tools for testing univariate minimizers
// Last edited on 2023-02-03 22:15:40 by stolfi

#ifndef minu_test_H
#define minu_test_H

#include <minu_gen.h>
#include <stdint.h>

typedef void (*ProblemFunc) (void *prb, double x, double *fx, double *dfx);

typedef double (*ProblemError) (void *prb, double x, double fx);

typedef double (*Function) (double x);

typedef bool_t (*GoodProc) (double x, double fx);

typedef struct Problem 
{
  char *name;                /* Problem's name. */
  ProblemFunc eval;          /* Function to optimize. */
  ProblemError error;        /* Error evaluator. */
  double xMin; double xMax;  /* Nominal range of "x". */
  double yMin; double yMax;  /* Range of "f(x)", for plotting. */
  double xStart;             /* Starting "x" value */
  double tol;                /* Required precision */
  double dist;               /* Estimated distance from "xStart" to minimum */
  unsigned maxCalls;         /* Maximum number of calls to "f" */
} Problem;

typedef struct Minimizer 
{
  char *name;                /* Name of method */
  MinimizeProc minimize;     /* Minimizer procedure */
  bool_t usesDerivative;       /* TRUE if method uses derivative info. */
} Minimizer;
    
typedef struct Performance 
{
  unsigned nTests;     /* Number of tests run */
    /* The standard deviations are significative only if "nTests >= 2". */
  double avgCalls;     /* Average number of function calls per test */
  double devCalls;     /* Standard deviation of function calls per test */
  unsigned maxCalls;   /* Max number of function calls per test. */
     /* Negative error is interpreted as zero error: */    
  double avgError;     /* Average solution error per test */
  double devError;     /* Standard deviation of solution error per test */
  float maxError;      /* Max solution error among all tests */ 
  unsigned nFailures;  /* Number of times that error was positive. */
} Performance;

// PROCEDURES 

Performance minu_test_single
  ( PSStream *ps,         /* Postscript file */
    int32_t page,          /* Page number in document */
    Minimizer *opt,    /* The minimization tool */
    Problem *prb,      /* Function and parameters */
    bool_t debug         /* Passed to the minimizer */
  );
/* 
  Runs the `minimizer' on the `problem' and prints the results.
  Also appends to `psf' a new page with a plot of the function,
  showing the `minimizer''s probes. */
  
Performance minu_test_multiple
  ( Minimizer *opt,    /* The minimization tool */
    Problem *prb,      /* Function and parameters */
    unsigned nTests    /* Number of tests to perform. */
  );
/* 
  Runs `nTests' tests with random starting points. 
  Prints averages etc. */

#endif
