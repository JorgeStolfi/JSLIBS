#define PROG_NAME "testjsqroots"
#define PROG_DESC "test of {jsqroots.h}"
#define PROG_VERS "1.0"

/* Last edited on 2008-07-19 20:11:57 by stolfi */ 
/* Created on 2007-01-02 by J. Stolfi, UNICAMP */

#define testjsmath_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <values.h>

#include <jsqroots.h>
#include <jsrandom.h>
#include <affirm.h>
#include <bool.h>

#define bug(FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

double zrandom(void); /* A random double with random magnitude. */

void test_roots_quadratic(int nt);

#define N_double_nice 10

static double double_nice[N_double_nice] = 
  { 
    4.45014771701440277e-308, /* 8/DBL_MAX. */                   
    7.45834073120020674e-155, /* 1/sqrt(DBL_MAX). */
    5.00000000000000000e-001,
    9.99999999999999889e-001,
    1.00000000000000000e+000, 
    1.00000000000000022e+000,
    2.00000000000000000e+000,
    1.34078079299425971e+154, /* sqrt(DBL_MAX). */
    2.24711641857789464e+307, /* DBL_MAX/8. */                   
    0
  };

int main (int argn, char **argv)
  { test_roots_quadratic(200);  
    return 0;
  }

bool_t check_roots_quadratic
  ( double A, double B, double C,                /* Coefficients. */
    double r1T, double r2T, double imT, int sdT, /* "True" solution. */
    double r1C, double r2C, double imC, int sdC, /* Computed solution. */
    bool_t verbose
  );

void test_roots_quadratic(int nt)
  { fprintf(stderr, "Checking {roots_quadratic}...\n");
    
    /* Some simple cases: */
    double A, B, C;
    double r1C, r2C, imC;
    int sdC;
    
    A = 00.0; B = 00.0; C = 00.0;
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)check_roots_quadratic(A, B, C, NAN, NAN, NAN, 00, r1C, r2C, imC, sdC, TRUE);
    
    A = +1.0; B = 00.0; C = 00.0;
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)check_roots_quadratic(A, B, C, 00.000, 00.000, 00.000, 00, r1C, r2C, imC, sdC, TRUE);
    
    A = +1.0; B = 00.0; C = -1.0;
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)check_roots_quadratic(A, B, C, +1.000, -1.000, 00.000, +1, r1C, r2C, imC, sdC, TRUE);
    
    A = +1.0; B = 00.0; C = -16.0;
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)check_roots_quadratic(A, B, C, +4.000, -4.000, 00.000, +1, r1C, r2C, imC, sdC, TRUE);
    
    A = +1.0; B = 00.0; C = +1.0;
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)check_roots_quadratic(A, B, C, 00.000, 00.000, +1.000, -1, r1C, r2C, imC, sdC, TRUE);
    
    A = +1.0; B = -2.0; C = +1.0;
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)check_roots_quadratic(A, B, C, +1.000, +1.000, 00.000, 00, r1C, r2C, imC, sdC, TRUE);
    
    A = +1.0; B = +2.0; C = +1.0;
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)check_roots_quadratic(A, B, C, -1.000, -1.000, 00.000, 00, r1C, r2C, imC, sdC, TRUE);
    
    A = +1.0; B = -1.0; C = +1.0;
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)check_roots_quadratic(A, B, C, +0.500, +0.500, 00.000, sqrt(3)/2, r1C, r2C, imC, sdC, TRUE);
    
    int i,j,k,sd,sx,sy;
    int nbug = 0; /* Number of buggy cases. */
    int ngud = 0; /* Number of checked-ok cases. */
    for (i = 0; i < N_double_nice + nt; i++)
      { /* Choose the absolute midpoint {M} of the two roots: */
        double M = (i < nt ? zrandom() : double_nice[i-nt]);
        for (j = 0; j < N_double_nice + nt; j++)
          { /* Choose the squared half-spacing of the roots: */
            double V = (j < nt ? zrandom() : double_nice[j-nt]);
            /* Define coeffs {A,B,C} and roots {r1,r2,im}: */
            double RR, RI; /* Real and imag half-spacing between the roots. */
            /* Consider all three signs of the determinant: */
            for (sd = -1; sd <= +1; sd++)
              { if (sd > 0)
                  { RR = sqrt(V); RI = 0.0; }
                else if (sd < 0)
                  { RR = 0.0; RI = sqrt(V); }
                else
                  { RR = 0.0; RI = 0.0; }
                for (k = 0; k < N_double_nice + nt; k++)
                  { /* Choose the abs value {S} of the {A} coefficient: */
                    double S = (k < nt ? zrandom() : double_nice[k-nt]);
                    /* Variously flip the axes: */
                    for (sx = -1; sx <= +1; sx += 2)
                      { for (sy = -1; sy <= +1; sy += 2)
                          { /* Define the equation and the "true" solutions: */
                            double r1 = (M - RR) * sx;
                            double r2 = (M + RR) * sx;
                            double im = RI;
                            double A = S * sy;
                            double B = -2*M * S * sy * sx;
                            double C = (r1*r2 + im*im) * S * sy;
                            /* Try solving the quadratic: */
                            double r1C = NAN, r2C = NAN, imC = NAN;
                            int sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
                            /* Check against the expected solution: */
                            bool_t vb = (ngud < 10);
                            bool_t ok = check_roots_quadratic(A, B, C, r1, r2, im, sd, r1C, r2C, imC, sdC, vb);
                            if (! ok) { nbug++; } else { ngud++; }
                            if (nbug > 50) { fprintf(stderr, "aborted!\n"); exit(1); }
                          }
                      }
                  }
              }
          }
      }
  }
  
#define roots_quadratic_REL_TOLERANCE (1.0e-10)
#define roots_quadratic_BIAS (1.0e-153/roots_quadratic_REL_TOLERANCE)

bool_t check_roots_quadratic
  ( double A, double B, double C,                    /* Coefficients. */
    double r1T, double r2T, double imT, int sdT, /* "True" solution. */
    double r1C, double r2C, double imC, int sdC, /* Computed solution. */
    bool_t verbose
  )
  {
    if (verbose)
      { fprintf(stderr, "check_roots_quadratic\n");
        fprintf(stderr, "  A  = %g  B = %g  C = %g\n", A, B, C);
      }

    /* Sort the real parts to simplify the comparison: */
    if (r1T > r2T) { double t = r1T; r1T = r2T; r2T = t; }
    if (r1C > r2C) { double t = r1C; r1C = r2C; r2C = t; }
  
    double r1E = r1C - r1T;
    double r2E = r2C - r2T;
    double imE = imC - imT;
  
    double r1R = 2 * r1E/(fabs(r1T) + fabs(r1C) + roots_quadratic_BIAS);
    double r2R = 2 * r2E/(fabs(r2T) + fabs(r2C) + roots_quadratic_BIAS);
    double imR = 2 * imE/(fabs(imT) + fabs(imC) + roots_quadratic_BIAS);
    
    double maxR = fmax(fmax(fabs(r1R), fabs(r2R)), fabs(imR));
    if ((maxR > roots_quadratic_REL_TOLERANCE) || (sdT*sdC < 0))
      { fprintf(stderr, "  ** maxR = %12.4e  sd = %+02d --> %+02d\n", maxR, sdT, sdC);
        fprintf(stderr, "  A  = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C);
        fprintf(stderr, "  r1 = %24.16e --> %24.16e  err = %24.16e  rel = %24.16e\n", r1T, r1C, r1E, r1R);
        fprintf(stderr, "  r2 = %24.16e --> %24.16e  err = %24.16e  rel = %24.16e\n", r2T, r2C, r2E, r2R);
        fprintf(stderr, "  im = %24.16e --> %24.16e  err = %24.16e  rel = %24.16e\n", imT, imC, imE, imR);
        return FALSE;
      }
    else
      { return TRUE; }
  }
  
#define LOG_DBL_MAX (709.78271289338397)

double zrandom(void)
  { 
    return drandom()*exp(LOG_DBL_MAX*(2*drandom() - 1));
  }
