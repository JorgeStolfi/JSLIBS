/* See pst_vertex_map_expand.h */
/* Last edited on 2025-02-28 06:36:16 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <pst_basic.h>

#include <pst_vertex_map_expand.h>

float_image_t *pst_vertex_map_expand
  ( float_image_t *B,
    int32_t wch,
    int32_t NXA,
    int32_t NYA,
    double scale
  )
  { 
    bool_t debug = FALSE;
  
    int32_t NC, NXB, NYB;
    float_image_get_size(B, &NC, &NXB, &NYB);
    
    if (debug) { fprintf(stderr, "creating expanded image %d×%d ... \n", NXB, NYB); }
    float_image_t *A = float_image_new(NC, NXA, NYA);

    for (int32_t yA = 0; yA < NYA; yA++)
      { int32_t yB = yA/2;
        int32_t ny = yA % 2;
        for (int32_t xA = 0; xA < NXA; xA++)
          { int32_t xB = xA/2;
            int32_t nx = xA % 2;
            bool_t bad = FALSE;
            float vA[NC];
            if ((nx == 0) && (ny == 0))
              { /* Just copy the pixel, scaled: */
                assert((xB >= 0) && (xB < NXB) && (yB >= 0) && (yB < NYB));
                float_image_get_pixel(B, xB, yB, vA);
                for (int32_t c = 0; c < NC; c++)
                  { if (c != wch) { vA[c] = (float)(vA[c]*scale); } }
              }
            else
              { /* Take average of nearby pixels: */
                double sum_wB = 0; /* Sum of {wB*wH} */
                double sum_wB_vB[NC]; /* Sum of {wB*wH*vB[0]} */
                for (int32_t c = 0; c < NC; c++) { sum_wB_vB[c] = 0; }
                for (int32_t ry = 0; ry <= ny; ry++)
                  { int32_t yBr = yB + ry;
                    for (int32_t rx = 0; rx <= nx; rx++)
                      { int32_t xBr = xB + rx;
                        if ((xBr >= 0) && (xBr < NXB) && (yBr >= 0) && (yBr < NYB))
                          { float vB[NC];
                            float_image_get_pixel(B, xBr, yBr, vB);
                            double wB = ((wch >= 0) && (wch < NC) ? vB[wch] : 1.0);
                            demand(isfinite(wB) && (wB >= 0), "invalid weight value");
                            if (wB > 0)
                              { sum_wB += wB;
                                for (int32_t c = 0; c < NC; c++)
                                  { sum_wB_vB[c] += wB*vB[c]; }
                              }
                          }
                      }
                  }
                bad = bad | (sum_wB == 0);
                for (int32_t c = 0; c < NC; c++)
                  { if (c == wch)
                      { vA[c] = (float)(sum_wB == 0 ? 0.0 : sum_wB);
                        if ((! isfinite(vA[c])) || (vA[c] == 0)) { bad = TRUE; }
                      }
                    else
                      { vA[c] = (float)(sum_wB == 0 ? NAN : scale*sum_wB_vB[c]/sum_wB);
                        if (! isfinite(vA[c])) { bad = TRUE; }
                      }
                  }

              }
            pst_ensure_pixel_consistency(NC, wch, bad, vA);
            float_image_set_pixel(A, xA, yA, vA);
          }
      }
    return A;
  }
