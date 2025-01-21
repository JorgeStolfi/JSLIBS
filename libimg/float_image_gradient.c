/* See {float_image_gradient.h}. */
/* Last edited on 2025-01-18 12:59:29 by stolfi */

#include <math.h>
#include <assert.h>
#include <values.h>
#include <string.h>

#include <bool.h>
#include <affirm.h>
#include <wt_table.h>
#include <wt_table_binomial.h>
#include <float_image.h>
#include <float_image_gradient.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

void float_image_gradient_sobel(float_image_t *A, int32_t cA, float_image_t *DX, int32_t cX, float_image_t *DY, int32_t cY)
  {
    /* Get the image dimensions: */
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    demand((cA >= 0) && (cA < A->sz[0]), "invalid {A} channel");
    
    bool_t has_DX = ((cX >= 0) && (DX != NULL));
    if (has_DX)
      { float_image_check_size(DX, -1, NX, NY, "bad {DX} image");
        demand((cX >= 0) && (cX < DX->sz[0]), "invalid {DX} channel");
      }
    
    bool_t has_DY = ((cY >= 0) && (DY != NULL));
    if (has_DY)
      { float_image_check_size(DY, -1, NX, NY, "bad {DY} image");
        demand((cY >= 0) && (cY < DY->sz[0]), "invalid {DY} channel");
      }
      
    if ((! has_DX) && (! has_DY)) { return; }
    
    /* Fill images: */
    int32_t x, y;
    for (x = 0; x < NX; x++)
      { for (y = 0; y < NY; y++)
          { /* Compute the horizontal and vertical derivatives in channel {c}: */
            double fx, fy;
            float_image_get_gradient_sobel(A, cA, x, y, &fx, &fy);
            /* Store the gradient: */
            if (has_DX) { float_image_set_sample(DX, cX, x, y, (float)fx); }
            if (has_DY) { float_image_set_sample(DY, cY, x, y, (float)fy); }
          }
      }
  }

void float_image_gradient_sqr_sobel(float_image_t *A, int32_t cA, float_image_t *G, int32_t cG)
  {
    /* Get the image dimensions: */
    int32_t NCA = (int32_t)A->sz[0];
    demand((cA >= 0) && (cA < NCA), "invalid {A} channel");
    int32_t NCG = (int32_t)G->sz[0];
    demand((cG >= 0) && (cG < NCG), "invalid {G} channel");
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    float_image_check_size(G, -1, NX, NY, "bad gradient image");
    
    /* Fill it: */
    int32_t x, y;
    for (x = 0; x < NX; x++)
      { for (y = 0; y < NY; y++)
          { /* Compute the horizontal and vertical derivatives in channel {c}: */
            double fx, fy;
            float_image_get_gradient_sobel(A, cA, x, y, &fx, &fy);
            /* Store the gradient squared: */
            double g2 = fx*fx + fy*fy;
            float_image_set_sample(G, cG, x, y, (float)g2);
          }
      }
  }

float_image_t *float_image_gradient_sqr_relative
  ( float_image_t *A, 
    double noise, 
    bool_t mono
  )
  {
    /* Get the image dimensions: */
    int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);
    
    /* Allocate the output image: */
    float_image_t *R = float_image_new((mono ? 1 : NC), NX, NY);
    
    /* Scratch images for gradient squared and variance: */
    float_image_t *G2 = float_image_new(1, NX, NY);
    float_image_t *V2 = float_image_new(1, NX, NY);
    
    /* Weight table: */
    int32_t hw = 2, nw = 2*hw+1;
    double wt[nw];
    wt_table_binomial_fill((uint32_t)nw, wt, NULL);
    wt_table_normalize_sum((uint32_t)nw, wt);
 
    /* Process channel by channel: */
    double noise2 = noise*noise;
    double scale = (mono ? 1.0 : 1.0/NC);
    int32_t c;
    for (c = 0; c < NC; c++)
      { /* Compute the gradient squared and variance for channel {c}: */
        float_image_gradient_sqr_sobel(A, c, G2, 0);
        float_image_local_avg_var(A, c, hw, wt, NULL, 0, V2, 0);
        /* Divive one by the other, store/add in {R}: */
        int32_t x, y;
        for (x = 0; x < NX; x++)
          { for (y = 0; y < NY; y++)
              { double g2 = float_image_get_sample(G2, 0, x, y);
                double v2 = float_image_get_sample(V2, 0, x, y);
                double q = g2/(v2 + noise2)*scale;
                float *rP = float_image_get_sample_address(R, (mono ? 0 : c), x, y);
                if ((c == 0) || (! mono))
                  { (*rP) = (float)q; }
                else
                  { (*rP) = (float)((*rP) + q); }
              }
          }
      }
        
    float_image_free(G2); G2 = NULL;
    float_image_free(V2); V2 = NULL;

    return R;
  }
  
