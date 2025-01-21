/* See float_image_expand_by_one.h */
/* Last edited on 2025-01-18 17:23:16 by stolfi */ 

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <float_image.h>

#include <float_image_expand_by_one.h>

float_image_t *float_image_expand_by_one(float_image_t *A, int32_t cW)
  {
    int32_t NCA, NXA, NYA;
    float_image_get_size(A, &NCA, &NXA, &NYA);
    
    int32_t NXE = NXA + 1;
    int32_t NYE = NYA + 1;
    float_image_t *E = float_image_new(NCA, NXE, NYE);
    
    for (int32_t y = 0; y < NYE; y++)
      { for (int32_t x = 0; x < NXE; x++)
          { double sum_wa[NCA];  /* Weighted sum of existing {A} pixels. */
            double sum_w = 0;    /* Sum of weights of existing {A} pixels. */
            double sum_m = 0;    /* Sum of reciprocals of existing {A} pixels. */
            uint32_t num_w = 0;  /* Number of existing {A} pixels. */
            for (int32_t c = 0; c < NCA; c++) { sum_wa[c] = 0; }
            for (int32_t dx = 0; dx <= 1; dx++)
              { for (int32_t dy = 0; dy <= 1; dy++)
                  { int32_t x1 = x - dx;
                    int32_t y1 = y - dy;
                    if ((x1 >= 0) && (x1 < NXA) && (y1 >= 0) && (y1 < NYA))
                      { double w = ((cW >= 0) && (cW < NCA) ? float_image_get_sample(A, cW, x1, y1) : 1.0);
                        demand(isfinite(w) && (w >= 0), "invalid weight");
                        for (int32_t c = 0; c < NCA; c++)
                          { if (c != cW)
                              { double a = float_image_get_sample(A, c, x1, y1);
                                sum_wa[c] += w*a;
                              }
                          }
                        sum_w += w;
                        sum_m += 1.0/w;
                        num_w++;
                      }
                  }
               }
            for (int32_t c = 0; c < NCA; c++)
              { if (c != cW) 
                  { double aexp = (sum_w == 0 ? 0.0 : sum_wa[c]/sum_w);
                    float_image_set_sample(E, c, x, y, (float)aexp);
                  }
                else
                  { double wexp = (num_w == 0 ? 0.0 : num_w/sum_m); /* Harmonic mean. */
                    float_image_set_sample(E, c, x, y, (float)wexp);
                  }
              }
          }
      }
    return E;
  }
