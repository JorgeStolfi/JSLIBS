#define PROG_NAME "testjsqroots"
#define PROG_DESC "test of {jsqroots.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-15 19:15:34 by stolfi */ 
/* Created on 2018-06-30 by J. Stolfi, UNICAMP */

#define testjsqroots_COPYRIGHT \
  "Copyright © 2018  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <jsqroots.h>
#include <jsrandom.h>
#include <affirm.h>
#include <bool.h>

double tjsr_random_double(void); /* A random double with random magnitude. */

void tjsr_one_test_from_roots(double M, double r1T, double r2T, double imT, int32_t sdT, bool_t verbose);
  /* Tests {roots_quadratic} on the equation {A*x^2 + B*x + C = 0} with expected 
    roots {r1T,r2T,imT} and discriminant sign {sdT}.
    The coefficient {A} of {x^2} is set to {M}, and {B,C}
    are scaled accordingly. */

bool_t tjsr_check_roots
  ( double A, double B, double C,                /* Coefficients. */
    double r1T, double r2T, double imT, int32_t sdT, /* "True" solution. */
    double r1C, double r2C, double imC, int32_t sdC, /* Computed solution. */
    bool_t verbose
  );
  /* Compares the computed roots {r1C,r2C,imC} and computed sign 
    of discriminant {sdC} of the equation {A*x^2 + B*x + C = 0}
    against the "true" (expected) values {r1T,r2T,imT,sdT}. */

void tjsr_compare(double reC, double imC, double reT, double imT, double *xE, double *xR);
  /* Compares {reC + imC*I} to {reT + imT*I} where {I = sqrt(-1)},
    and returns the absolute and relative Euclidean half-errors in {*xE,*xR}. 
    Takes {NAN}s into account. */

void tjsr_do_multiple_tests(int32_t nt);

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

int32_t main (int32_t argn, char **argv)
  { tjsr_do_multiple_tests(200);  
    return 0;
  }

#define tjsr_MAX_COEF (DBL_MAX/8.0)
  /* Max abs value of coefficients {A,B,C} to use in mass tests. */
  
#define tjsr_MIN_COEF (DBL_MIN*8.0)
  /* Min abs value of nonzero coefficients {A,B,C} to use in mass tests. */

/* 
  TOLERANCES 
  
  The root finders should be improved so that we can make these
  tolerances more stringent. */

#define tjsr_ROOTS_REL_TOL (1.5e-8)
  /* Max relative half-error allowed in roots. */
  
#define tjsr_ROOTS_TINY (1.0e-69)
  /* Max root that is allowed to be approximated to zero... */
  
#define tjsr_REL_ERR_BIAS (tjsr_ROOTS_TINY/tjsr_ROOTS_REL_TOL)
  /* Max relative error allowed in roots. */

void tjsr_compare(double reC, double imC, double reT, double imT, double *xE, double *xR)
  { 
    if ((isnan(reC) != isnan(reT)) || (isnan(imC) != isnan(imT)))
      { /* Disagree about {NAN}s in real or imaginary part: */
        (*xE) = (*xR) = 1.0;
      }
    else 
      { /* Absolute half-errors in each part and whole: */
        double reE = (isnan(reC) ? 0.0 : fabs(reC/2 - reT/2));
        double imE = (isnan(imC) ? 0.0 : fabs(imC/2 - imT/2));
        /* Half-magnitudes: */
        double magC = hypot((isnan(reC) ? 0.0 : reC/2), (isnan(imC) ? 0.0 : imC/2));
        double magT = hypot((isnan(reT) ? 0.0 : reT/2), (isnan(imT) ? 0.0 : imT/2));
        (*xE) = hypot(reE, imE);
        (*xR) = (*xE)/(magC + magT + tjsr_REL_ERR_BIAS);
      }
  }  

void tjsr_do_multiple_tests(int32_t nt)
  { fprintf(stderr, "Checking {roots_quadratic}...\n");
    
    /* Some simple cases: */
    double A, B, C;
    double r1T, r2T, imT;
    double r1C, r2C, imC;
    int32_t sdC, sdT;
    
    /* Check some degenerate cases: */
    
    A = 00.0; B = INF; C = 00.0;
    r1T = NAN; r2T = NAN; imT = NAN; sdT = 00; 
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)tjsr_check_roots(A, B, C,  r1T, r2T, imT, sdT,  r1C, r2C, imC, sdC, TRUE);
    
    A = 00.0; B = 00.0; C = 00.0;
    r1T = NAN; r2T = NAN; imT = 0.00; sdT = 00; 
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)tjsr_check_roots(A, B, C,  r1T, r2T, imT, sdT,  r1C, r2C, imC, sdC, TRUE);
    
    A = 00.0; B = 00.0; C = 10.0;
    r1T = NAN; r2T = NAN; imT = 0.00; sdT = 00; 
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)tjsr_check_roots(A, B, C,  r1T, r2T, imT, sdT,  r1C, r2C, imC, sdC, TRUE);
    
    A = 00.0; B = 16.0; C = 25.0;
    r1T = -25.0/16.0; r2T = NAN; imT = 0.00; sdT = +1; 
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)tjsr_check_roots(A, B, C,  r1T, r2T, imT, sdT,  r1C, r2C, imC, sdC, TRUE);
    
    A = 00.0; B = 23.0; C = 0.0;
    r1T = 0.0; r2T = NAN; imT = 0.00; sdT = +1; 
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)tjsr_check_roots(A, B, C,  r1T, r2T, imT, sdT,  r1C, r2C, imC, sdC, TRUE);
    
    A = -1.4806137525181331e+125; B = -8.8356704178207576e-29; C = -1.3181876704774955e-182;
    r1T = -2.9837864205953832e-154; r2T = -2.9837864205953832e-154; imT = 0.00; sdT = 00; 
    r1C = r2C = imC = 4615.17031947;
    sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)tjsr_check_roots(A, B, C,  r1T, r2T, imT, sdT,  r1C, r2C, imC, sdC, TRUE);
    
    for (int32_t sgnM = -1; sgnM <= +1; sgnM += 2)
      { double M = 8*sgnM;
    
        /* Check some cases with real distinct roots: */

        tjsr_one_test_from_roots(M,  -1.000, +1.000, 00.000, +1, TRUE); /* {A = +1.0; B = 00.0; C = -1.0;} */
        tjsr_one_test_from_roots(M,  -4.000, +4.000, 00.000, +1, TRUE); /* {A = +1.0; B = 00.0; C = -16.0;} */

        /* Check some cases with real double-roots: */

        tjsr_one_test_from_roots(M,  00.000, 00.000, 00.000, 00, TRUE); /* {A = +1.0; B = 00.0; C = 00.0;} */
        tjsr_one_test_from_roots(M,  00.000, 00.000, 00.000, 00, TRUE); /* {A = +1.0; B = 00.0; C = 00.0;} */
        tjsr_one_test_from_roots(M,  +1.000, +1.000, 00.000, 00, TRUE); /* {A = +1.0; B = -2.0; C = +1.0;} */
        tjsr_one_test_from_roots(M,  -1.000, -1.000, 00.000, 00, TRUE); /* {A = +1.0; B = +2.0; C = +1.0;} */

        /* Check some cases with complex conjugate roots: */

        tjsr_one_test_from_roots(M,  00.000, 00.000, +1.000, -1, TRUE); /* {A = +1.0; B = 00.0; C = +1.0;} */
        tjsr_one_test_from_roots(M,  +2.000, +2.000, +1.000, -1, TRUE); /* {A = +1.0; B = -4.0; C = +5.0;} */
      }

    /* Checks random non-degenerate cases: */
    int32_t i,j,k,sd,sx,sy;
    int32_t nbug = 0; /* Number of buggy cases. */
    int32_t ngud = 0; /* Number of checked-ok cases. */
    for (i = 0; i < N_double_nice + nt; i++)
      { /* Choose the absolute midpoint {M} of the two roots: */
        double M = (i < nt ? tjsr_random_double() : double_nice[i-nt]);
        for (j = 0; j < N_double_nice + nt; j++)
          { /* Choose the squared half-spacing of the roots (non-zero): */
            double V = (j < nt ? tjsr_random_double() : double_nice[j-nt]);
            if (fabs(V) < DBL_MIN) { V = (V >= 0 ? DBL_MIN : -DBL_MIN); }
            for (sd = -1; sd <= +1; sd++)
              { /* Choose the imaginary part {RI} and the half-separation {RR} of the real parts: */
                double RR, RI; /* Real and imag half-spacing between the roots. */
                if (sd > 0)
                  { RR = sqrt(V); 
                    /* Make sure that {RR} is large enough compared to {M}: */
                    RR = fmax(RR, 1.0e-10 * fabs(M));
                    RI = 0.0;
                  }
                else if (sd < 0)
                  { RR = 0.0; RI = sqrt(V); }
                else
                  { RR = 0.0; RI = 0.0; }
                /* Try various scaling factors for all coeffs: */
                for (k = 0; k < N_double_nice + nt; k++)
                  { /* Choose a scaling factor {S}: */
                    double S = (k < nt ? tjsr_random_double() : double_nice[k-nt]);
                    /* Variously flip the axes: */
                    for (sx = -1; sx <= +1; sx += 2)
                      { for (sy = -1; sy <= +1; sy += 2)
                          { /* Define the true roots {r1T,r2T,imT}: */
                            double r1T = M * sx - RR;
                            double r2T = M * sx + RR;
                            double imT = RI;
                            /* Last-minute paranoia: */
                            if (sd > 0)  { assert(r1T < r2T); assert(RI == 0.0); }
                            if (sd == 0) { assert(r1T == r2T); assert(RI == 0.0); }
                            if (sd < 0)  { assert(r1T == r2T); assert(RI > 0.0); }
                            /* Define coeffs {A,B,C} from the roots and axis signs: */
                            double A = sy;
                            double B = -2*M * sx * sy;
                            double C = (r1T*r2T + imT*imT) * sy;
                            /* Find the max absolute coef value: */
                            double CoefMax = fmax(fabs(A), fmax(fabs(B), fabs(C)));
                            /* If the coefficients are way too big or way too small, skip this test: */
                            if (CoefMax <= tjsr_MAX_COEF)
                              { /* Find the min nonzero absolute coef value: */
                                double CoefMin = +INF;
                                if (fabs(A) != 0) { CoefMin = fmin(CoefMin, fabs(A)); }
                                if (fabs(B) != 0) { CoefMin = fmin(CoefMin, fabs(B)); }
                                if (fabs(C) != 0) { CoefMin = fmin(CoefMin, fabs(C)); }
                                if (CoefMin >= tjsr_MIN_COEF)
                                  { /* Adjust the scaling {S} so that it will not overflow or underflow: */
                                    double R = S;
                                    if (fabs(R)*CoefMax > tjsr_MAX_COEF) 
                                      { R = tjsr_MAX_COEF/CoefMax; }
                                    if ((CoefMin != +INF) && (fabs(R)*CoefMin < tjsr_MIN_COEF))
                                      { R = tjsr_MIN_COEF/CoefMin; }
                                    /* Scale the equation: */
                                    A = R*A; B = R*B; C = R*C;
                                    /* Try solving the quadratic: */
                                    double r1C = NAN, r2C = NAN, imC = NAN;
                                    int32_t sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
                                    /* Check against the expected solution: */
                                    bool_t vb = (ngud < 10);
                                    bool_t ok = tjsr_check_roots(A, B, C, r1T, r2T, imT, sd, r1C, r2C, imC, sdC, vb);
                                    if (! ok) { nbug++; } else { ngud++; }
                                    if (nbug > 50) { fprintf(stderr, "aborted!\n"); exit(1); }
                                  }
                              }
                          }
                      }
                  }
              }
          }
      }
  }
  
void tjsr_one_test_from_roots(double M, double r1T, double r2T, double imT, int32_t sdT, bool_t verbose)
  { double A, B, C;
    assert(! isnan(r1T));
    assert(! isnan(r2T));
    assert(! isnan(imT));
    assert((sdT >= -1) && (sdT <= +1));
    if (sdT > 0)
      { /* Two real roots: */
        assert(r1T < r2T);
        assert(imT == 0.0);
        A = M;
        B = -(r1T + r2T)*M;
        C = r1T*r2T*M;
      }
    else if (sdT < 0)
      { /* Two conjugate imaginary roots: */
        assert(r1T == r2T);
        assert(imT > 0.0);
        A = M;
        B = -2*r1T*M;
        C = (r1T*r1T + imT*imT)*M;
      }
    else
      { /* Two equal real roots: */
        assert(r1T == r2T);
        assert(imT == 0.0);
        A = M;
        B = -2*r1T*M;
        C = r1T*r1T*M;
      }
    double r1C, r2C, imC;
    r1C = r2C = imC = 4615.17031947; /* To check for undef returns. */
    int32_t sdC = roots_quadratic(A, B, C, &r1C, &r2C, &imC);
    (void)tjsr_check_roots(A, B, C, r1T, r2T, imT, sdT, r1C, r2C, imC, sdC, verbose);
  }

bool_t tjsr_check_roots
  ( double A, double B, double C,                /* Coefficients. */
    double r1T, double r2T, double imT, int32_t sdT, /* "True" solution. */
    double r1C, double r2C, double imC, int32_t sdC, /* Computed solution. */
    bool_t verbose
  )
  {
    if (verbose)
      { fprintf(stderr, "tjsr_check_roots  A = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C); }

    bool_t bug = FALSE; /* Set when we find something wrong. */
    bool_t warn = FALSE; /* Set when we find something questionable. */
    
    /* For convenience, sort the true roots: */
    if ((! isnan(r2T)) && (isnan(r1T) || (r1T > r2T)))
      { double t = r1T; r1T = r2T; r2T = t; }

    /* Check true root order and sign of imaginary part: */
    if (sdT > 0) 
      { assert(isnan(r2T) || (r1T < r2T)); assert(isnan(imT) || (imT == 0.0)); } 
    else if (sdT < 0)
      { assert(isnan(r2T) || (r1T == r2T)); assert(isnan(imT) || (imT > 0.0)); }
    else
      { assert(isnan(r2T) || (r1T == r2T)); assert(isnan(imT) || (imT == 0.0)); }
    
    /* It is common for {stC} and {sdT} to disagree, so we do not compare them. */

    /* Check consistency of roots and discriminant sign: */
    if (isnan(r1C))
      { if (! isnan(r2C))
          { fprintf(stderr, "** inconsistent use of NAN:  sdC = %+02d  r1C = %24.16e  r2C = %24.16e\n", sdC, r1C, r2C);  bug = TRUE; }
      }
    if (sdC > 0)
      { if ((!isnan(imC)) && (imC != 0.0))
          { fprintf(stderr, "** nonzero imaginary part:  sdC = %+02d  imC = %24.16e\n", sdC, imC);  bug = TRUE; }
        if ((! isnan(r2C)) && (r1C >= r2C))
          { fprintf(stderr, "** roots not sorted:  sdC = %+02d  r1C = %24.16e  r2C = %24.16e\n", sdC, r1C, r2C);  bug = TRUE; }
      }
    else if (sdC < 0.0)
      { if ((!isnan(imC)) && (imC <= 0.0))
          { fprintf(stderr, "** non-positive imaginary part:  sdC = %+02d  imC = %24.16e\n", sdC, imC);  bug = TRUE; }
        if ((! isnan(r2C)) && (r1C != r2C))
          { fprintf(stderr, "** real parts not equal:  sdC = %+02d  r1C = %24.16e  r2C = %24.16e\n", sdC, r1C, r2C);  bug = TRUE; }
      }
    else
      { if ((!isnan(imC)) && (imC != 0.0))
          { fprintf(stderr, "** nonzero imaginary part:  sdC = %+02d  imC = %24.16e\n", sdC, imC);  bug = TRUE; }
        if ((! isnan(r2C)) && (r1C != r2C))
          { fprintf(stderr, "** real parts not equal:  sdC = %+02d  r1C = %24.16e  r2C = %24.16e\n", sdC, r1C, r2C);  bug = TRUE; }
      }

    /* Compute absolute and relative errors in roots: */
    double r1E, r2E;  /* Absolute errors, or 1.0 if disagree about NAN. */
    double r1R, r2R;  /* Relative errors, or 1.0 if disagree about NAN. */
    tjsr_compare(r1C, +imC, r1T, +imT, &r1E, &r1R);
    tjsr_compare(r2C, -imC, r2T, -imT, &r2E, &r2R);

   /* Maximum relative error in all roots: */
    double maxR = fmax(fabs(r1R), fabs(r2R));
    
    double tolR = tjsr_ROOTS_REL_TOL;
    
    if (maxR > tolR)
      { fprintf(stderr, "** relative error too large:  maxR = %12.4e  tolR = %12.4e\n", maxR, tolR);  bug = TRUE; }
    
    /* Print out data if there was any bug: */
    if (bug || warn)
      { fprintf(stderr, "  ---------------------------------------------------------\n");
        fprintf(stderr, "  A =   %24.16e  B =  %24.16e   C =   %24.16e\n", A, B, C);
        fprintf(stderr, "  sdC = %+02d    sdT = %+02d\n", sdC, sdT);
        fprintf(stderr, "  imC = %24.16e  imT = %24.16e\n", imC, imT);
        fprintf(stderr, "  r1C = %24.16e  r1T = %24.16e  err = %24.16e  rel = %24.16e\n", r1C, r1T, r1E, r1R);
        fprintf(stderr, "  r2C = %24.16e  r2T = %24.16e  err = %24.16e  rel = %24.16e\n", r2C, r2T, r2E, r2R);
        fprintf(stderr, "  ---------------------------------------------------------\n");
        if (! bug)
          { /* Just a warning: */
            fprintf(stderr, "  !! just a warning\n");
          }
        else if ((isnan(r1C) && isnan(r2C)) || isnan(imC))
          { /* Bug in a degenerate case, do not consider. */
            bug = FALSE; warn = TRUE;
            fprintf(stderr, "  !! error ignored\n");
          }
        else
          { /* Bug in a non-degenerate case: */
            fprintf(stderr, "  ** check failed\n");
          }
      }
    else
      { /* Check succeeded: */
        if (verbose) { fprintf(stderr, "  check succeeded\n"); }
      } 
    return (! bug);
  }
  
#define LOG_DBL_MAX (709.78271289338397)

double tjsr_random_double(void)
  { 
    return drandom()*exp(LOG_DBL_MAX*(2*drandom() - 1));
  }
