/* See tbfind.h */
/* Last edited on 2024-11-15 19:15:23 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <affirm.h>
#include <assert.h>
#include <bool.h>

#include <tbfind.h>

#define TB_DEBUG FALSE

int32_t tb_find(double f(int32_t i), int32_t iMin, int32_t iMax)
  {
    int32_t i;
    if (TB_DEBUG) fprintf(stderr, "\n[");
    double fMin = (iMin > iMax ? 0 : f(iMin));
    double fMax = (iMax <= iMin ? fMin : f(iMax));
    if (0 <= fMin )
      { i = iMin; }
    else if (fMax < 0)
      { i = iMax+1; }
    else 
      { /* Now {f(iMin) < 0 <= f(iMax)} hence {iMin < iMax}. */
        /* Aswer must be in {iMin+1..iMax}. */
        bool_t use_lin = TRUE; /* TRUE uses linear interpolation, FALSE uses bisection. */
        while (iMax - iMin >= 2)
          { /* Now {f(iMin) < 0 <= f(iMax)} and there are at least two intervals. */
            if (TB_DEBUG) fprintf(stderr, "(%d:%d)", iMin,iMax);
            
            int32_t nOld = iMax - iMin; /* Number of intervals remaining. */
            
            int32_t iTry;
            if (use_lin)
              { /* Get an index {iTry} in {iMin+1..iMax-1} by linear linterpolation: */
                double r = (0 - fMin)/(fMax - fMin);
                iTry = iMin + (int32_t)(rint(r*(iMax - iMin)));
                if (iTry <= iMin) { iTry = iMin+1; }
                if (iTry >= iMax) { iTry = iMax-1; }
              }
            else
              { /* Get the middle element: */
                iTry = (iMin + iMax)/2;
              }
            assert(iTry > iMin);
            assert(iTry < iMax);

            if (TB_DEBUG) fprintf(stderr, "?%d", iTry);
            
            /* Bisection step with pivot {iTry}: */
            double fTry = f(iTry);
            if (fTry < 0)
              { iMin = iTry; fMin = fTry; if (TB_DEBUG) fprintf(stderr, ">"); }
            else
              { iMax = iTry; fMax = fTry; if (TB_DEBUG) fprintf(stderr, "<"); }
            
            /* Decide what to do next: */
            use_lin = (! use_lin) | ((iMax - iMin) < 0.55*nOld);
          }
        /* Only one interval, and {f(iMin) < 0 <= f(iMax)}. */
        i = iMax;
      }

    if (TB_DEBUG) 
      { fprintf(stderr, "=%d=(", i); 
        if (iMin > iMax)
          { fprintf(stderr, "-oo_+oo"); }
        else
          { if (i <= iMin) { fprintf(stderr, "-oo"); } else { fprintf(stderr, "%3f", f(i-1)); } 
            fprintf(stderr, "_"); 
            if (i > iMax) { fprintf(stderr, "+oo"); } else { fprintf(stderr, "%3f", f(i)); } 
          }
        fprintf(stderr, ")"); 
        fprintf(stderr, "]\n"); 
      }
    return i;
  }
