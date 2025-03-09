/* See pst_cell_map_shrink.h */
/* Last edited on 2025-02-27 14:19:44 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <pst_basic.h>

#include <pst_cell_map_shrink.h>

float_image_t *pst_cell_map_shrink
  ( float_image_t *A,
    int32_t wch,
    int32_t NXB,
    int32_t NYB,
    double scale
  )
  { 
    int32_t NC, NXA, NYA;
    float_image_get_size(A, &NC, &NXA, &NYA);

    /* Allocate the reduced map: */
    float_image_t *B = float_image_new(NC, NXB, NYB);

    /* Fill the pixels of {B}: */
    for (int32_t yB = 0; yB < NYB; yB++)
      { for (int32_t xB = 0; xB < NXB; xB++)
          { double sum_wA = 0; 
            double sum_wA_vA[NC];
            for (int32_t c = 0; c < NC; c++) { sum_wA_vA[c] = 0; }
            bool_t bad = FALSE;
            for (int32_t ry = 0; ry <= 1; ry++)
              { int32_t yA = 2*yB + ry;
                for (int32_t rx = 0; rx <= 1; rx++)
                  { int32_t xA = 2*xB + rx;
                    if ((xA < NXA) && (yA < NYA))
                      { float vA[NC];
                        float_image_get_pixel(A, xA, yA, vA);
                        /* Get the reliability weight, if any: */
                        double wA = ((wch >= 0) && (wch < NC) ? vA[wch] : 1.0); 
                        demand(isfinite(wA) && (wA >= 0), "invalid wIeight value");
                        if ((! isfinite(vA[0])) || (! isfinite(vA[1]))) { wA = 0; }
                        if (wA == 0)
                          { bad = TRUE; }
                        else
                          { /* Accumulate weighted sums: */
                            sum_wA += wA;
                            for (int32_t c = 0; c < NC; c++)
                              { sum_wA_vA[c] += wA*vA[c];
                                if (! isfinite(vA[c])) { bad = TRUE; }
                              }
                          }
                      }
                  }
              }
            /* Store the weighted average (maybe {NAN}) in the {B} pixel: */
            float vB[NC];
            for (int32_t c = 0; c < NC; c++)
              { if (c == wch)
                  { vB[c] = (float)(bad || sum_wA == 0 ? 0.0 : sum_wA);
                    if ((! isfinite(vB[c])) || (vB[c] == 0)) { bad = TRUE; }
                  }
                else
                  { vB[c] = (float)(bad || sum_wA == 0 ? NAN : scale*sum_wA_vA[c]/sum_wA); 
                    if (! isfinite(vB[c])) { bad = TRUE; }
                  }
              }
              
            pst_ensure_pixel_consistency(NC, wch, bad, vB);
            float_image_set_pixel(B, xB, yB, vB);
          }
      }
    return B;
  }
