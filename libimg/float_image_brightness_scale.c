/* See float_image_brightness_scale.h */
/* Last edited on 2025-02-02 02:21:26 by stolfi */ 

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>

#include <float_image_brightness_scale.h>

float_image_t *float_image_brightness_scale_add(float_image_t *A, bool_t loY, float vMin, float vMax)
  {
    /* Get the size of {A}: */
    int32_t NCA, NXA, NYA;
    float_image_get_size(A, &NCA, &NXA, &NYA);
    
    /* Choose the height of the scale: */
    int32_t NYS_min = 6; /* Min scale height. */
    int32_t NYS = NYA/10;  /* Actual scale height. */
    if (NYS < NYS_min) { NYS = NYS_min; }
    
    int32_t NXR = NXA, NYR = NYA + NYS; 
    float_image_t *R = float_image_new(NCA, NXR, NYR);
    
    int32_t ylo_scale = (loY ? 0 : NYA);     /* First row of scale in {R}. */
    int32_t yhi_scale = ylo_scale + NYS - 1; /* Last row of scale in {R}. */

    int32_t ylo_orig = (loY ? NYS : 0);      /* First row of {A} copy in {R}. */
    int32_t yhi_orig = ylo_orig + NYA - 1;   /* Last row of {A} copy in {R}. */
    float_image_assign_rectangle(R, 0,NXA-1, ylo_orig,yhi_orig, A, 0,0);
    float_image_brightness_scale_paint(R, 0,NXR-1, ylo_scale,yhi_scale, vMin,vMax);
    return R;
  }

void float_image_brightness_scale_paint
  ( float_image_t *A,
    int32_t xMin, int32_t xMax,  
    int32_t yMin, int32_t yMax,  
    float vMin, float vMax
  )
  {
    int32_t NCA, NXA, NYA;
    float_image_get_size(A, &NCA, &NXA, &NYA);

    /* Get size of scale: */
    int32_t NXS = xMax - xMin + 1;  /* Width of scale. */
    int32_t NYS = yMax - yMin + 1;  /* Heights of scale. */
    
    /* Choose the count, width and height of the chips: */
    /* Valid chip counts are 10+1, 5+1, 4+1, 2+1, 0+1. */
    int32_t NXP_min = 2; /* Min chip width. */
    int32_t NP = 11; /* Number of chips in scale. */
    int32_t NXP = NXS/NP;
    if (NXP < NXP_min) { NP = 6; NXP = NXS/NP; }
    if (NXP < NXP_min) { NP = 5; NXP = NXS/NP; }
    if (NXP < NXP_min) { NP = 3; NXP = NXS/NP; }
    if (NXP < NXP_min) { NP = 1; NXP = NXS/NP; }
    assert((NP == 1) || (NP >= 3)); /* A scale with {NP=2} is useless. */
    
    int32_t xMrg = (NP*NXP_min + (NP-1) < NXA ? 1 : 0); /* Gray margin between chips. */
    /* Adjust chip width for margin: */
    NXP = (NXS - (NP-1)*xMrg)/NP;
    assert(NXP >= NXP_min);
    int32_t xMrg0 = (NXS - NP*NXP - (NP-1)*xMrg)/2;  /* Margin before first chip. */

    int32_t yMrg = (NYS < 6 ? 0 : 1); /* Margin above and below chips. */

    float vMid = (vMin+vMax)/2;
    for (int32_t c = 0; c < NCA; c++)
      { /* Fill whole scale space with middling value: */
        float_image_fill_channel_rectangle(A, c, xMin, xMax, yMin, yMax, vMid);
        /* If there is a single chip, it is already painted with the proper value. So: */
        if (NP > 1)
          { for (int32_t ip = 0; ip < NP; ip++)
              { int32_t xpMin = xMrg0 + ip*(NXP + xMrg);
                int32_t xpMax = xpMin + NXP - 1;
                double rp = ((double)ip)/((double)NP-1), sp = 1-rp;
                float vp = (float)(sp*vMin + rp*vMax);
                float_image_fill_channel_rectangle(A, c, xpMin,xpMax, yMin+yMrg,yMax-yMrg, vp);
              }
          }
      }
  }   
