/* See {float_image_overlay.h}. */
/* Last edited on 2021-08-25 10:53:24 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
 
#include <bool.h>
#include <affirm.h>
#include <float_image.h>

#include <float_image_overlay.h>

void float_image_overlay(float_image_t *A, float_image_t *B, int32_t icop)
  {
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    
    demand((int32_t)B->sz[0] == NC, "channel count mismatch");
    demand((int32_t)B->sz[1] == NX, "column count mismatch");
    demand((int32_t)B->sz[2] == NY, "row count mismatch");
    
    demand((icop >= 0) && (icop < NC), "invalid opacity chanel index {icop}");
    
    for (int32_t iy = 0; iy < NY; iy++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { float *oAP = float_image_get_sample_address(A, icop, ix, iy);
            float *oBP = float_image_get_sample_address(B, icop, ix, iy);
            double rB = (*oBP), rA = (1.0 - rB)*(*oAP);
            double oC = rB + rA;
            if (oC != 0)
              { double sB = rB/oC, sA = rA/oC; 
                for (int32_t ic = 0; ic < NC; ic++)
                  { if (ic == icop)
                      { float *vAP = float_image_get_sample_address(A, ic, ix, iy);
                        double vA = (*vAP);
                        double vB = float_image_get_sample(B, ic, ix, iy);
                        (*vAP) = (float)(sB*vB + sA*vA);
                      }
                    else
                      { *(oAP) = oC; }
                  }
              }
          }
      }
  }
