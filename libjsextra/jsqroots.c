/* See jsqroots.h */
/* Last edited on 2008-07-16 03:02:23 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <bool.h>
#include <assert.h>

#include <jsmath.h>
#include <jsqroots.h>

int32_t roots_quadratic(double A, double B, double C, double *r1, double *r2, double *im)
  {
    /* Garbage in, garbage out: */
    if ((! isfinite(A)) || (! isfinite(B)) || (! isfinite(C)))
      { (*r1) = (*r2) = (*im) = NAN; return 00; }
      
    /* Check for degenerate cases: */
    if (A == 0)
      { /* Linear equation {B*x + C == 0}: */
        (*im) = 0.0;
        if (B == 0)
          { /* Non-equation {C == 0}: */
            (*r1) = (*r2) = NAN;
            return 00;
          }
        else
          { if (C != 0)
              { /* Non-degenerate linear equation {B*x + C == 0}: */
                (*r1) = -C/B; /* May overflow or underflow... */
                (*r2) = NAN;
              }
            else
              { /* Special linear equation {B*x == 0}: */
                (*r1) = 0.0; 
                (*r2) = NAN;
              }
            return +1;
          }
      }
    else if (C == 0)
      { /* Special quadratic equation {A*x^2 + B*x == 0}, with {A != 0}: */
        (*r1) = 0.0; 
        (*r2) = -B/A; /* May overflow or underflow... */
        (*im) = 0.0;
        return +1;
      }
    else
      { /* Non-degenerate quadratic equation {A*x^2 + B*x + C == 0}: */ 
        return roots_proper_quadratic(A, B, C, r1, r2, im);
      }
  }
  
#define rq_MAX_DIS_EXP (1024)
#define rq_MIN_DIS_EXP (-900)
  /* Max and min exponents in discriminant. */

int32_t roots_proper_quadratic(double A, double B, double C, double *r1, double *r2, double *im)
  { 
    assert(isfinite(A));
    assert(isfinite(B));
    assert(isfinite(C));
    
    assert(A != 0);
    assert(C != 0);

    /* Determine the binary exponents of {A,C}: */
    int32_t eA, eC;
    (void)frexp(A, &eA);
    (void)frexp(C, &eC);

    /* Rescale {A = A/2^eS}, {C = C*2^eS} so that they have about the same magnitude: */
    int32_t eS = eA - (eA + eC)/2;
    if (eS != 0)
      { fprintf(stderr, "scaling {A,C} by 2^{±%d}\n", eS);
        fprintf(stderr, "  A = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C);
        A = ldexp(A, -eS); C = ldexp(C, +eS);
        eA -= eS; eC += eS;
        fprintf(stderr, "  A = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C);
      }

    if (B == 0)
      { /* Special quadratic  equation {A*x^2 + C = 0}: */
        double x2 = -C/A; /* May overflow or underflow... */
        if (x2 > 0)
          { double x = ldexp(sqrt(x2), -eS);
            (*r1) = -x; (*r2) = +x; (*im) = 0.0;
            return +1;
          }
        else if (x2 < 0)
          { double x = ldexp(sqrt(-x2), -eS);
            (*r1) = (*r2) = 0.0; (*im) = x;
            return -1;
          }
        else
          { /* Should not happen: */ assert(FALSE); }
      }

    /* General case, with all coefs non-zero. */
    /* Make sure that the discriminant will not {over,under}flow: */
    int32_t eB;
    (void)frexp(B, &eB);
    int32_t eBB4 = 2*eB - 2, eAC = eA + eC;
    int32_t eU = (eBB4 > eAC ? eBB4 : eAC) + 1; /* Upper bound to exp of discriminant. */
    int32_t eT = (eU > rq_MAX_DIS_EXP ? rq_MAX_DIS_EXP : (eU < rq_MIN_DIS_EXP ? rq_MIN_DIS_EXP : eU));
    int32_t eN = (eT - eU)/2;
    if (eN != 0)
      { fprintf(stderr, "scaling {A,B,C} by 2^{%d}\n", eN);
        fprintf(stderr, "  A = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C);
        A = ldexp(A, eN); B = ldexp(B, eN); C = ldexp(C, eN);
        fprintf(stderr, "  A = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C);
      }

    /* Alea jacta est: */
    double dis = (B/4)*B - A*C; /* Discriminant/4. */

    if ((eN != 0) || (eS != 0))
      { fprintf(stderr, "  dis = %24.16e\n", dis); }

    assert(isfinite(dis));

    if (dis > 0)
      { /* Two distinct real roots: */
        double q = (B/2) + copysign(sqrt(dis),B);
        (*r1) = - ldexp(C / q, -eS);
        (*r2) = - ldexp(q / A, -eS);
        (*im) = 0.0;
        return +1;
      }
    else if (dis < 0) 
      { /* Two complex conjugate roots: */
        (*r1) = (*r2) = -ldexp(B/A, -eS -1);
        (*im) = ldexp(sqrt(dis), -eS + 1);
        return -1;
      }
    else
      { /* Two equal real roots, {-B/A}: */
        (*r1) = (*r2) = -ldexp(B/A, -eS - 1);
        (*im) = 0.0;
        return 00;
      }
  }
