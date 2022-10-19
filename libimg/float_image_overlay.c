/* See {float_image_overlay.h}. */
/* Last edited on 2021-08-28 03:49:07 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
 
#include <bool.h>
#include <jsmath.h>
#include <affirm.h>
#include <float_image.h>

#include <float_image_overlay.h>

void float_image_overlay(float_image_t *A, float_image_t *B, int32_t icop, int32_t xlo, int32_t ylo)
  {
    int32_t NCA = (int32_t)A->sz[0];
    int32_t NXA = (int32_t)A->sz[1];
    int32_t NYA = (int32_t)A->sz[2];
    
    int32_t NCB = (int32_t)B->sz[0];
    int32_t NXB = (int32_t)B->sz[1];
    int32_t NYB = (int32_t)B->sz[2];
    
    demand(NCA == NCB, "channel counts do not match");
    demand((icop >= 0) && (icop < NCA), "invalid opacity channel index {icop}");

    int32_t xini = (int32_t)imax(0, xlo);
    int32_t xfin = (int32_t)imin(NXA-1, xlo+NXB-1);
    
    int32_t yini = (int32_t)imax(0, ylo);
    int32_t yfin = (int32_t)imin(NYA-1, ylo+NYB-1);
    
    for (int32_t iy = yini; iy <= yfin; iy++)
      { for (int32_t ix = xini; ix <= xfin; ix++)
          { float *oAP = float_image_get_sample_address(A, icop, ix, iy);
            float oA = (*oAP);
            float oB = float_image_get_sample(B, icop, ix, iy);
            double rB = oB, rA = (1.0 - rB)*oA;
            double oC = rB + rA;
            if (oC != 0)
              { double sB = rB/oC, sA = rA/oC; 
                for (int32_t ic = 0; ic < NCA; ic++)
                  { if (ic != icop)
                      { float *vAP = float_image_get_sample_address(A, ic, ix, iy);
                        double vA = (*vAP);
                        double vB = float_image_get_sample(B, ic, ix - xlo, iy - ylo);
                        (*vAP) = (float)(sB*vB + sA*vA);
                      }
                    else
                      { *(oAP) = oC; }
                  }
              }
          }
      }
  }
