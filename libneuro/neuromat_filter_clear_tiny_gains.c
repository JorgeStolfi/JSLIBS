/* See {neuromat_filter_clear_tiny_gains.h}. */
/* Last edited on 2023-12-16 14:06:36 by stolfi */

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
  { int32_t kfa = nf/2; /* Indices of conjugate Hartley freqs. */
    /* Clear all small gains: */
    for (int32_t kfa = 0; kfa <= nf/2; kfa++) 
      { int32_t kfb = (nf - kfa) % nf;
        if ((H[kfa] != 0) || (H[kfb] != 0))
          { double Gabs = M_SQRT1_2*hypot(H[kfa], H[kfb]);
            double Gmult = ((Gabs - eps) < 1.0e-140 ? 0.0 : exp(1/(1-eps) - 1/(Gabs - eps)));
            double Ga = H[kfa] * Gmult, Gb = H[kfb] * Gmult;
            if ((Ga == 0) && (Gb == 0))
              { if (verbose)
                  { fprintf(stderr, "  cleared gains for %.8f Hz", ((double)kfa)*fsmp/nf);
                    fprintf(stderr, " H[%d] = %24.16e ", kfa, H[kfa]);
                    if (kfb != kfa) 
                      { fprintf(stderr, " H[%d] = %24.16e ", kfb, H[kfb]);
                        fprintf(stderr, " Gm = %24.16e\n", Gabs);
                      }
                    fprintf(stderr, "\n");
                  }
              }
            H[kfa] = Ga; H[kfb] = Gb;
          }
      }
    /* Get highest surviving frequency: */
    kfa = nf/2;
    { int32_t kfb = (nf - kfa) % nf;
      while ((kfa >= 0) && (H[kfa] == 0) && (H[kfb] == 0)) 
        { kfa--; kfb = (kfb + 1) % nf; }
      if (kfa < 0)
        { if (verbose) { fprintf(stderr, "filter gains are all zero\n"); } }
      else
        { if (verbose) 
            { double fa = ((double)kfa)*fsmp/nf;
              fprintf(stderr, "max preserved frequency = %.8f Hz", fa);
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
      }
    return kfa;
  }
