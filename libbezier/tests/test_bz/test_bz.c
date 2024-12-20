/* Last edited on 2024-12-05 10:22:13 by stolfi */
/* Test of bz_basic.h and bz_patch.h routines. */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>
#include <interval.h>

#include <bz_basic.h>

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);

void test_bezier(bz_degree_t g, double a, double b, double m);
  /* Tests of basic routines for polynomials of degree {g},
    over the interval {[a _ b]}, splitting the domain at value {m}. */
    
void test_bezier_aux(bz_degree_t g, double C[], double a, double b, double m, bz_index_t ix);
  /* Tests of basic routines for a polynomial {f} of degree {g} with Bézier
    coeffs {c[0..g]} over the interval {[a _ b]}, splitting the domain at value
    {m}. If {i} is non-negative, compares {f} to the Bernstein-Bézier polynomial
    with degree {g} and index {ix}. */
 
void evaluate
  ( char *name,      /* Polynoial name. */
    bz_degree_t g,   /* Degree of curve. */
    double c[],      /* Bezier coeffs. */
    double u,        /* Start of interval. */
    double v,        /* End of interval */
    double x,        /* Argument value. */
    bool_t debug,    /* Print debugging info. */
    bz_degree_t ord, /* Max derivative order desired. */
    double f[]       /* OUT: Image vector (size {n}). */
  );
  /* Calls {bz_eval} with the given parameters; then, if {debug} is true,
    prints the arguments (incluindg the coeffs {c[0..g]} and the results
    {f[0..ord]}. */

void print_coeffs(bz_degree_t g, double c[]);
  /* Prints Bézier coeffs {c[0..g]} to {stderr}. */
  
void print_results(bz_degree_t ord, double f[]);
   /* Prints results {f[0..ord]} to {stderr}. */
  
void check_val(double va, double vb, double tol, char *msg, int32_t kk);
  /* Compare{va} and {vb}, bombs out with error {msg,kk} if
    they are not practically equal. */
 
/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { 
    double a = 2.0;
    double b = 5.0;
    bz_degree_t g;
    bz_index_t kx;
    for (g = 0; g <= 4; g++)
      { for (kx = -1; kx <= +1; kx++)
          { double x = 0.5*(a+b) + 0.25*(b-a)*kx;
            test_bezier(g, a, b, x);
          }
      }
    fprintf(stderr, "=== DONE ===\n");
    return(0);
  }

void test_bezier(bz_degree_t g, double a, double b, double m)
  {
    double C[g+1];
    bz_index_t i, j;
    /* Bernstein-Bézier polynomials: */
    for (i = 0; i <= g; i++)
      { for (j = 0; j <= g; j++) { C[j] = 0; }
        C[i] = 1;
        test_bezier_aux(g, C, a, b, m, i);
      }
      
    /* A random polynomial: */
    for (j = 0; j <= g; j++) { C[j] = drandom(); }
    test_bezier_aux(g, C, a, b, m, -1);
        
    fprintf(stderr, "\n");
  }
  
void test_bezier_aux(bz_degree_t g, double C[], double a, double b, double m, bz_index_t ix)
  { 
    bool_t debug = FALSE;
    
    bz_degree_t gD = (bz_degree_t)(g >= 1 ? g-1 : 0);  /* Degree of 1st derivative, clipped at 0. */
    bz_degree_t gDD = (bz_degree_t)(g >= 2 ? g-2 : 0); /* Degree of 2nd derivative, clipped at 0. */
    
    double CDD[gDD+1];  /* Bézier coeffs of 2nd derivative. */ 
    double CD[gD+1];    /* Bézier coeffs of 1st derivative. */
    double CDI[g+1];    /* Bézier coeffs of integral of derivative. */
    double C0[g+1];     /* Bézier coeffs of first piece. */
    double C1[g+1];     /* Bézier coeffs of second piece. */
    
    int32_t NS = 35; /* Sampling intervals. */
    int32_t i, j, k;
     
    fprintf(stderr, "--- testing [ %7.4f _ %7.4f _ %7.4f ] g = %d ix = %d ---\n", a, m, b, g, ix);
        
    if (debug)
      { fprintf(stderr, "coeffs C = \n");
        print_coeffs(g,C);
        fprintf(stderr, "\n");
      }
    /*  */
    bz_split(g, C, m - a, C0, b - m, C1);
    if (debug)
      { fprintf(stderr, "coeffs C0 = \n");
        print_coeffs(g,C0);
        fprintf(stderr, "coeffs C1 = \n");
        print_coeffs(g,C1);
      }
    /* Some basic checks: */
    affirm(C0[0] == C[0],  "ini coeff bug");
    affirm(C1[g] == C[g],  "fin coeff bug");
    affirm(C0[g] == C1[0], "join coeff bug");
    
    if (g >= 0)
      { /* fprintf(stderr, "--- testing {bz_diff} ---"); */
        bz_diff(g, C, b - a, CD);
        if (debug)
          { fprintf(stderr, "coeffs CD = \n");
            print_coeffs(gD,CD);
          }
      }
    
    if (g >= 1)
      { /* fprintf(stderr, "--- testing {bz_diff} again ---"); */
        bz_diff(gD, CD, b - a, CDD);
        if (debug)
          { fprintf(stderr, "coeffs CDD = \n");
            print_coeffs(gDD,CDD);
          }
      }
    
    /* fprintf(stderr, "--- testing {bz_integ} ---"); */
    bz_integ(gD, CD, b - a, CDI);
    for(j = 0; j <= g; j++) { CDI[j] += C[0]; }
    if (debug)
      { fprintf(stderr, "coeffs CDI = \n");
        print_coeffs(g,CDI);
      }
    
    /* fprintf(stderr, "--- testing {bz_eval} ---"); */
    
    bz_degree_t ord = (bz_degree_t)(g+1);     /* Max derivative order for function. */
    bz_degree_t ordD = (bz_degree_t)(gD+1);   /* Max derivative order for 1st derivative. */
    bz_degree_t ordDD = (bz_degree_t)(gDD+1); /* Max derivative order for 2nd derivative. */
    
    double f[ord+1], fDI[ord+1], f0[ord+1], f1[ord+1];
    double fD[ordD+1], fDD[ordDD+1];
    
    double tol = 1.0e-12; /* Accuracy for zero-order inerpolation. */
    
    for (i = 0; i <= NS; i++)
      { double x = a + (b - a)*((double)i)/((double)NS);
        
        evaluate("C", g, C, a, b, x, debug, ord, f);
        if (i == 0)  { check_val(f[0], C[0], tol, "eval - initial value", -1); }
        if (i == NS) { check_val(f[0], C[g], tol, "eval - final value",   -1); }
          
        evaluate("C0",  g, C0,  a, m, x, debug, ord, f0);
        evaluate("C1",  g, C1,  m, b, x, debug, ord, f1);
        evaluate("CDI", g, CDI, a, b, x, debug, ord, fDI);
        if (g >= 1) { evaluate("CD",  gD,  CD,  a, b, x, debug, ordD,  fD); }
        if (g >= 2) { evaluate("CDD", gDD, CDD, a, b, x, debug, ordDD, fDD); }
        
        double tolf = tol;
        for (k = 0; k <= ord; k++)
          { tolf = (k+1)*tolf;
            check_val(f[k], f0[k],  tolf, "split/eval - first half - order %d", k);
            check_val(f[k], f1[k],  tolf, "split/eval - second half - order %d", k);
            check_val(f[k], fDI[k], tolf, "diff/integ/eval - order %d", k);
            if ((g >= 1) && (k >= 1)) { check_val(f[k], fD[k-1],  tolf, "diff/eval - order %d", k); }
            if ((g >= 2) && (k >= 2)) { check_val(f[k], fDD[k-2], tolf, "diff/diff/eval - order %d", k); }
          }
        
        if (ix >= 0)
          { double fbb = bz_bernstein(g, ix, (x - a)/(b - a));
            check_val(fbb, f[0], tol, "eval/bernstein", -1);
          }
      }
  }

void evaluate
  ( char *name,      /* Polynoial name. */
    bz_degree_t g,   /* Degree of curve. */
    double c[],      /* Bezier coeffs. */
    double u,        /* Start of interval. */
    double v,        /* End of interval */
    double x,        /* Argument value. */
    bool_t debug,    /* Print debugging info. */
    bz_degree_t ord, /* Max derivative order desired. */
    double f[]       /* OUT: Image vector (size {n}). */
  )
  {
    bz_eval(g, c,u,v,x, ord,f);
    if (debug)
      { fprintf(stderr, "evaluating %s\n", name);
        fprintf(stderr, "domain = [ %18.15f _ %18.15f ] arg = %18.15f coeffs = \n", u, v, x);
        print_coeffs(g, c);
        fprintf(stderr, "results = \n");
        print_results(ord, f);
        fprintf(stderr, "\n");
      }
  }

void print_coeffs(bz_degree_t g, double c[])
  { int32_t i;
    for (i = 0; i <= g; i++) { fprintf(stderr, "%3d %18.15f\n", i, c[i]); }
  }

void print_results(bz_degree_t ord, double f[])
  { int32_t i;
    for (i = 0; i <= ord; i++) { fprintf(stderr, "%3d %18.15f\n", i, f[i]); }
  }


void check_val(double va, double vb, double tol, char *msg, int32_t kk)
  { 
    if (fabs(va - vb) > tol)
      { char *fmsg = NULL;
        if (kk >= 0)
          { char *fmsg = jsprintf(msg, kk); }
        else
          { char *fmsg = jsprintf(msg); }
        fprintf(stderr, "** %24.16e != %24.16e %s\n", va, vb, fmsg);
        assert(FALSE);
        free(fmsg);
      }
  }
