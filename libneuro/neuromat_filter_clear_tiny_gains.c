/* See {neuromat_filter_clear_tiny_gains.h}. */
/* Last edited on 2024-01-05 18:10:42 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_clear_tiny_gains.h>
       
int32_t neuromat_filter_clear_tiny_gains(int32_t nf, double H[], double eps, double fsmp, bool_t verbose)
  { demand((eps >= 0) && (eps <= 0.5), "invalid {eps}");
    if (eps > 0)
      { 
        /* Clear all small gains: */
        double beta = 0.9; /* 0.5; */
        double logA = -beta/(1-eps);
        for (int32_t kfa = 0; kfa <= nf/2; kfa++) 
          { int32_t kfb = (nf - kfa) % nf;
            if ((H[kfa] != 0) || (H[kfb] != 0))
              { double Ha, Hb;  /* Adjusted {H[kfa],H[kfb]}. */
                double Gabs = M_SQRT1_2*hypot(H[kfa], H[kfb]);
                if (Gabs < eps + 1.0e-140)
                  { /* Correct to  zero: */
                    Ha = Hb = 0;
                  }
                else
                  { /* Correct smoothly down: */
                    double logGmult = -(logA + beta*Gabs/(Gabs-eps));
                    double Gmult = exp(logGmult);
                    if ((verbose) && ((Gabs < 10*eps) || (Gabs >= 0.9)))
                      { fprintf(stderr, "  Gabs = %24.16e logGmult = %24.16e Gmult = %24.16e\n", Gabs, logGmult, Gmult); }
                    Ha = H[kfa] * Gmult, Hb = H[kfb] * Gmult;
                  }
                if ((Ha == 0) && (Hb == 0))
                  { if (verbose)
                      { fprintf(stderr, "  cleared gains for %.8f Hz", ((double)kfa)*fsmp/nf);
                        fprintf(stderr, " H[%d] = %24.16e ", kfa, H[kfa]);
                        if (kfb != kfa) 
                          { fprintf(stderr, " H[%d] = %24.16e ", kfb, H[kfb]);
                            fprintf(stderr, " |H| = %24.16e", Gabs);
                          }
                        fprintf(stderr, "\n");
                      }
                  }
                H[kfa] = Ha; H[kfb] = Hb;
              }
          }
      }
    /* Get highest surviving frequency: */
    { int32_t kfa = nf/2, kfb = (nf - kfa) % nf;
      while ((kfa >= 0) && (H[kfa] == 0) && (H[kfb] == 0)) { kfa--; kfb = (kfb + 1) % nf; }
      if (kfa < 0)
        { if (verbose) { fprintf(stderr, "  filter gains are all zero\n"); } }
      else
        { if (verbose) 
            { double fa = ((double)kfa)*fsmp/nf;
              fprintf(stderr, "  max preserved frequency = %.8f Hz", fa);
              if (kfb == kfa) 
                { fprintf(stderr, " gain H[%d] = %+24.16e\n", kfa, H[kfa]); }
              else
                { double fb = ((double)kfb)*fsmp/nf;
                  fprintf(stderr, " = %.8f Hz", fb);
                  double Hmax = M_SQRT1_2*hypot(H[kfa], H[kfb]);
                  fprintf(stderr, " gain = %24.16e", Hmax);
                  fprintf(stderr, "  H[%d] = %+24.16e  H[%d] = %+24.16e\n", kfa, H[kfa], kfb, H[kfb]);
                }
            }
        }
      return kfa;
    }
  }
