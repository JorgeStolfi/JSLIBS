/* See pst_vertex_map_shrink.h */
/* Last edited on 2025-02-27 14:42:35 by stolfi */

#include <math.h>
#include <assert.h>
#include <stdint.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <pst_basic.h>

#include <pst_vertex_map_shrink.h>

float_image_t *pst_vertex_map_shrink
  ( float_image_t *A,
    int32_t wch,
    int32_t NXB,
    int32_t NYB,
    double scale
  )
  {
    bool_t debug = FALSE;
    
    int32_t NC, NXA, NYA;
    float_image_get_size(A, &NC, &NXA, &NYA);
    
    /* Create the output image {B}: */
    if (debug) { fprintf(stderr, "creating shrunk image %d×%d ... \n", NXB, NYB); }
    float_image_t *B = float_image_new(NC, NXB, NYB);
    
    /* Fill the pixels of {B}: */
    for (int32_t yB = 0; yB < NYB; yB++)
      { for (int32_t xB = 0; xB < NXB; xB++)
          { /* Accumulate the weighted pixel sum over the window.
              For {c!=wch}, the value of {sum_wv[c]} is the weighted 
              sum of the samples {A[c,X',Y']} where {x'=2*x+rx},
              {y'=2*y+ry}, and {rx,ry} range in {-1..1}. The weight
              is {wA*wH[rx]*wH[ry]} where {wA} is the {A[wch,X',Y']}
              and {wH[-1..+1]=(1/2,1,1/2)}. */
            double sum_wA_wH = 0;   /* Sum of {wA*wH[rx]*wH[ry]}. */
            double sum_wA_wH2 = 0;  /* Sum of {wA*(wH[rx]*wH[ry])^2}. */
            double sum_wA_wH_vA[NC];  /* Sum of {wA*wH[rx]*wH[ry]} times height. */
            for (int32_t c = 0; c < NC; c++) { sum_wA_wH_vA[c] = 0; }
            for (int32_t ry = -1; ry <= +1; ry++)
              { int32_t yA = 2*yB + ry;
                for (int32_t rx = -1; rx <= +1; rx++)
                  { int32_t xA = 2*xB + rx;
                    if ((xA >= 0) && (xA < NXA) && (yA >= 0) && (yA < NYA))
                      { float vA[NC];
                        float_image_get_pixel(A, xA, yA, vA);
                        /* Get the reliability weight, if any: */
                        double wA = ((wch >= 0) && (wch < NC) ? vA[wch] : 1.0); 
                        demand(isfinite(wA) && (wA >= 0), "invalid weight value");
                        if (wA > 0)
                          { /* Accumulate weighted sums: */
                            double wH = (rx == 0 ? 1 : 0.5)*(ry == 0 ? 1 : 0.5);
                            sum_wA_wH += wA*wH; 
                            sum_wA_wH2 += wA*wH*wH; 
                            for (int32_t c = 0; c < NC; c++)
                              { sum_wA_wH_vA[c] += wA*wH*vA[c]; }
                          }
                      }
                  }
              }
            /* Store the weighted average (maybe {NAN}) in the {B} pixel: */
            float vB[NC]; 
            bool_t bad = FALSE;
            for (int32_t c = 0; c < NC; c++)
              { if (c == wch)
                  { vB[c] = (float)(sum_wA_wH2 == 0 ? 0.0 : (sum_wA_wH*sum_wA_wH)/sum_wA_wH2);
                    if ((! isfinite(vB[c])) || (vB[c] == 0)) { bad = TRUE; }
                  }
                else
                  { vB[c] = (float)(sum_wA_wH == 0 ? NAN : scale*sum_wA_wH_vA[c]/sum_wA_wH); 
                    if (! isfinite(vB[c])) { bad = TRUE; }
                  }
              }
              
            pst_ensure_pixel_consistency(NC, wch, bad, vB);
            float_image_set_pixel(B, xB, yB, vB);
          }
      }
    return B;
  }
