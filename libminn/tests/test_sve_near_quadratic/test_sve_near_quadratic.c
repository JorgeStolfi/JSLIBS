/* test_sve_near_quadratic --- test of {sve_minn.h} for a nearly quadratic func */
/* Last edited on 2017-03-13 21:02:11 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#include <fftw3.h>

#include <bool.h>
#include <sign.h>
#include <argparser.h>
#include <jsrandom.h>
#include <jsstring.h>
#include <affirm.h>
#include <rmxn.h>
#include <rmxn_extra.h>
#include <rn.h>
#include <jsfile.h>
#include <vec.h>

#include <sve_minn.h>

/* GENERAL PARAMETERS */

typedef struct options_t
  { } options_t;

/* INTERNAL PROTOTYPES */

int main (int argc, char **argv);

options_t *get_options(int argc, char **argv);
  /* Parses the command-line options. */

void test_minimizer(int id, int n);
  /* Tests the minimizer on problem with index {id} and dimension {n}. */ 

void write_solution(char *prefix, char *tag, int n, sve_goal_t *F, double x[], double Fx);
  /* Writes the solution to problem with index {id} and dimension {n},
    also shows it to {stderr}, together with the function value {Fx}. */

void write_vector(char *prefix, char *tag, int n, double x[]);
  /* Writes the solution to problem with index {id} and dimension {n}. */

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { /* options_t *o = get_options(argc, argv); */
    int id;
    for (id = 0; id < 10; id++) 
      { int n = (id % 3) + 1;
        test_minimizer(id, n);
      }
    fclose(stderr);
    fclose(stdout);
    return (0);
  }
  
void test_minimizer(int id, int n)
  { 
    fprintf(stderr, "=== %s ===\n", __FUNCTION__);
    fprintf(stderr, "test id = %d\n", id);
    fprintf(stderr, "dimension = %d\n", n);
    
    /* Output file names: */
    char *prefix = NULL; /* Prefix for output file names. */
    asprintf(&prefix, "out/%03d-%02d", id, n);
    
    /* Shake the dice: */
    srand(4615 + id);  srandom(4615 + id);
    
    /* Working storage for the goal function: */
    double A[n*n]; /* A random linear matrix. */
    double M[n*n]; /* Matrix of quadric. */
    double dx[n];
    
    /* Choose the true minimum: */
    double x_tru[n];
    int k;
    for (k = 0; k < n; k++) { x_tru[k] = 2*drandom() - 1; }
    
    /* Fill the matrices: */
    int i, j;
    for (i = 0; i < n; i++)
      { for (j = 0; j < n; j++)
          { A[i*n + j] = 2*drandom() - 1; }
      }
    rmxn_mul_tr(n, n, n, A, A, M);
    
    auto double F(int n, double x[]); 
      /* The goal function for optimization. */
      
    double F(int n, double x[])
      { assert(n == n);
        /* Subtract the true minimum: */
        rn_sub(n, x, x_tru, dx);
        /* Apply a slight nonlinear deformation: */
        double r2 = rn_norm(n, dx);
        double scale = 1.0 + 0.0001*r2;
        rn_scale(n, scale, dx, dx);
        /* Evaluate the quadric at {dx}: */
        double S = 0;
        int i0, i1;
        for (i0 = 0; i0 < n; i0++)
          { for (i1 = 0; i1 < n; i1++)
              { S += dx[i0]*dx[i1]*M[i0*n + i1]; }
          }
        return S;
      }
      
    int nok = 0;      /* Counts iterations (actually, calls to {OK}). */
    
    auto bool_t OK(int n, double x[], double Fx); 
      /* Acceptance criterion function. */
      
    bool_t OK(int n, double x[], double Fx)
      { assert(n == n);
        fprintf(stderr, "iteration %d:\n", nok);
        write_solution(NULL, "tmp", n, &F, x, Fx);
        nok++;
        fprintf(stderr, "\n");
        return 0;
      }

    /* Output the true solution: */
    fprintf(stderr, "true minimum:\n");
    double Fx_tru = F(n, x_tru);
    write_solution(prefix, "tru", n, &F, x_tru, Fx_tru);

    double x[n];     /* Initial guess and final solution. */
    for (k = 0; k < n; k++) { x[k] = 2*drandom() - 1; }
    
    /* Evaluate and write the initial solution: */
    fprintf(stderr, "initial guess:\n");
    double Fx = F(n, x);
    write_solution(prefix, "ini", n, &F, x, Fx);
    
    /* Optimize iteratively: */
    double dMax = +INFINITY;
    bool_t dBox = FALSE;
    double rMin = 0.0000001;
    double rMax = 10.0;
    double rIni = 0.5;
    double stop = 0.01*rMin;
    sign_t dir = -1;
    int maxIters = 200;
    bool_t debug = TRUE;
    
    sve_minn_iterate(n, &F, &OK, x, &Fx, dir, dMax, dBox, rIni, rMin, rMax, stop, maxIters, debug);
    
    /* Evaluate and write the final solution: */
    fprintf(stderr, "final solution:\n");
    write_solution(prefix, "fin", n, &F, x, Fx);
  }

void write_solution(char *prefix, char *tag, int n, sve_goal_t *F, double x[], double Fx)
  { /* Print and write the pulse: */
    fprintf(stderr, "  point =\n");
    write_vector(NULL, tag, n, x); 

    double FxN = F(n, x);
    fprintf(stderr, "  function value = %+24.16e %+24.16e\n", Fx, FxN);
    demand(Fx == FxN, "inconsistent function value");
    /* if (prefix != NULL) { write_vector(prefix, tag, n, x); }a */
  }

void write_vector(char *prefix, char *tag, int n, double x[])
  { char *fname = NULL;
    FILE *wr;
    if (prefix != NULL) 
      { asprintf(&fname, "%s-%s-a.dat", prefix, tag);
        wr = open_write(fname, TRUE);
      }
    else
      { wr = stderr; }
    int k;
    for (k = 0; k < n; k++)
      { fprintf(wr, "%5d %12.8f\n", k, x[k]); }
    if ((wr != stderr) && (wr != stdout)) { fclose(wr); }
    if (fname != NULL) { free(fname); }
  }
  
options_t *get_options(int argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    argparser_finish(pp);
    return o;
  }
