/* See jsqroots.h */
/* Last edited on 2024-11-15 19:14:25 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <bool.h>
#include <assert.h>

#include <jsmath.h>
#include <jsqroots.h>

#define jsqroots_DEBUG FALSE

int32_t roots_quadratic(double A, double B, double C, double *r1, double *r2, double *im)
  {
    
    /* Check for degenerate cases otherwise {roots_proper_quadratic}. */
    /* Beware of damn minus zero! */
    
    if (jsqroots_DEBUG)
      { fprintf(stderr, "  -------------------------------------------------------\n");
        fprintf(stderr, "  A = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C);
      }
    
    int32_t discr; /* Discriminat sign to return. */
    if ((! isfinite(A)) || (! isfinite(B)) || (! isfinite(C)))
      { /* Garbage in, garbage out: */
        (*r1) = (*r2) = (*im) = NAN; discr = 00;
      }
    else if (fabs(A) == 0)
      { /* Linear equation {B*x + C == 0}: */
        (*im) = 0.0;
        if (fabs(B) == 0)
          { /* Non-equation {C == 0}: */
            (*r1) = (*r2) = NAN;
            discr = 00;
          }
        else
          { if (fabs(C) != 0)
              { /* Non-degenerate linear equation {B*x + C == 0}: */
                (*r1) = -C/B; /* May overflow or underflow; can't do anything... */
                (*r2) = NAN;
              }
            else
              { /* Special linear equation {B*x == 0}: */
                (*r1) = 0.0; 
                (*r2) = NAN;
              }
            discr = +1;
          }
      }
    else if (fabs(C) == 0)
      { /* Special quadratic equation {A*x^2 + B*x == 0}, with {A != 0}: */
        double x = -B/A; /* May overflow or underflow; can't do anything... */
        if (fabs(x) == 0)
          { /* Somehow we got a double root at 0; pretend the discr was 0: */
            (*r1) = (*r2) = 0;
            discr = 00;
          }
        else
          { if (x < 0.0)
              { (*r1) = x; (*r2) = 0; }
            else
              { (*r1) = 0.0; (*r2) = x; }
            discr = +1;
          }
        (*im) = 0.0;
      }
    else
      { /* Non-degenerate quadratic equation {A*x^2 + B*x + C == 0}: */ 
        discr = roots_proper_quadratic(A, B, C, r1, r2, im);
      }

    if (jsqroots_DEBUG)
      { fprintf(stderr, "  r1 = %24.16e  r2 = %24.16e  im = %24.16e  discr = %+02d\n", (*r1), (*r2), (*im), discr);
        fprintf(stderr, "  -------------------------------------------------------\n");
      }
    
    return discr;
  }
  
#define jsqroots_MAX_DISCR_EXP (1023)
#define jsqroots_MIN_DISCR_EXP (-900)
  /* Max and min exponents in discriminant. */

int32_t roots_proper_quadratic(double A, double B, double C, double *r1, double *r2, double *im)
  { 
    assert(isfinite(A));
    assert(isfinite(B));
    assert(isfinite(C));
    
    assert(A != 0);
    assert(C != 0);
    
    /* To simplify the code below: */
    (*r1) = (*r2) = (*im) = 0.0; 

    /* Determine the binary exponents of {A,C}: */
    int32_t eA; (void)frexp(A, &eA);
    int32_t eC; (void)frexp(C, &eC);

    /* Rescale {A = A/2^eS}, {C = C*2^eS} so that they have about the same magnitude: */
    int32_t eS = (eA - eC)/2;
    if (eS != 0)
      { if (jsqroots_DEBUG) { fprintf(stderr, "  scaling {A,C} by {2^%d,2^%d}\n", -eS, +eS); }
        A = ldexp(A, -eS); C = ldexp(C, +eS);
        eA -= eS; eC += eS;
        if (jsqroots_DEBUG) { fprintf(stderr, "  A = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C); }
        assert(isfinite(A) && (A != 0));
        assert(isfinite(C) && (C != 0));
      }

    /* Check for special case (beware of minus zero!): */
    if (fabs(B) == 0)
      { /* Special quadratic equation {A*x^2 + C = 0}. */
        assert(A != 0.0);
        assert(C != 0.0);
        
        /* Compute the absolute value {x} of the roots: */
        double x = sqrt(fabs(C/A)); /* Should not overflow and have small exponent. */
        if (jsqroots_DEBUG) { fprintf(stderr, "  |x| = %24.16e\n", x); }
        assert(x > 0.0);
        
        /* Compensate by the pre-scaling of {A} and {C}: */ 
        if (jsqroots_DEBUG) { fprintf(stderr, "  scaling roots by 2^%d\n", -eS); }
        x = ldexp(x, -eS); /* Could it overflow? */
        if (jsqroots_DEBUG) { fprintf(stderr, "  now |x| = %24.16e\n", x); }
        if (x == 0)
          { /* Somehow we got a double real root at zero: */
            return 00;
          }
        else if ((A < 0) == (C > 0))
          { /* Symmetric real roots: */
            if (jsqroots_DEBUG) { fprintf(stderr, "  symmetric real roots |x| = %24.16e\n", x); }
            (*r1) = -x; (*r2) = +x; 
            return +1;
          }
        else
          { /* Symmetric distinct imaginary roots: */
            if (jsqroots_DEBUG) { fprintf(stderr, "  pure imaginary  roots |x| = %24.16e\n", x); }
            (*im) = x;
            return -1;
          }
      }

    /* General case, with all coefs non-zero. */
    /* Make sure that the discriminant {B^2/4 - A*C} will not {over,under}flow: */
    int32_t eB; (void)frexp(B, &eB);
    int32_t eBB4 = 2*eB - 2, eAC = eA + eC;
    int32_t eU = (eBB4 > eAC ? eBB4 : eAC) + 1; /* Upper bound to exp of discriminant. */
    int32_t eTmax = jsqroots_MAX_DISCR_EXP;
    int32_t eTmin = jsqroots_MIN_DISCR_EXP;
    int32_t eT = (eU > eTmax ? eTmax : (eU < eTmin ? eTmin : eU));  /* Target dis exp after scaling. */
    int32_t eN = (eT - eU)/2; /* Scaling exp on coefs to get discr exp {eT}. */
    if (eN != 0)
      { if (jsqroots_DEBUG) { fprintf(stderr, "  scaling {A,B,C} by 2^{%d}\n", eN); }
        A = ldexp(A, eN); B = ldexp(B, eN); C = ldexp(C, eN);
        if (jsqroots_DEBUG) { fprintf(stderr, "  A = %24.16e  B = %24.16e  C = %24.16e\n", A, B, C); }
        assert(isfinite(A) && (A != 0));
        assert(isfinite(B) && (B != 0));
        assert(isfinite(C) && (C != 0));
      }

    /* Alea jacta est: */
    double discr = (B/4)*B - A*C; /* Discriminant/4. */
    if (jsqroots_DEBUG) { fprintf(stderr, "  discriminant = %24.16e\n", discr); }

    if (! isfinite(discr)) 
      { fprintf(stderr, "  ##  A = %24.16e  B = %24.16e  C = %24.16e  discr = %g\n", A, B, C, discr); } 
    assert(isfinite(discr));

    /* Solve the equation, then compensate for {A,C} prescaling: */
    if (discr > 0)
      { /* Two distinct real roots: */
        double q = (B/2) + copysign(sqrt(discr),B);
        double x1 = -ldexp(C/q, -eS); /* Could it overflow? */
        double x2 = -ldexp(q/A, -eS); /* Could it overflow? */
        if (x1 == x2)
          { /* Somehow we got two equal roots.  Pretend that the discr was zero: */
            (*r1) = (*r2) = x1;
            return 00;
          }
        else
          { /* Sort the roots: */
            if (x1 > x2) { double t = x1; x1 = x2; x2 = t; }
            (*r1) = x1; (*r2) = x2; 
            return +1;
          }
      }
    else if (discr < 0) 
      { /* Two complex conjugate roots: */
        double xre = ldexp(-B/A/2, -eS); /* Could it overflow? */
        double xim = ldexp(sqrt(-discr)/fabs(A), -eS); /* Could it overflow? */
        (*r1) = (*r2) = xre;
        if (xim == 0.0)
          { /* Somehow we got two equal roots.  Pretend that the discr was zero: */
            return 00;
          }
        else
          { (*im) = xim; 
            return -1;
          }
      }
    else
      { /* Two equal real roots, {-B/A/2}: */
        double x = ldexp(-B/A/2, -eS); /* Could it overflow? */
        (*r1) = (*r2) = x; 
        return 00;
      }
  }
