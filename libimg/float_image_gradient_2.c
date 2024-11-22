/* See {float_image_gradient_2.h}. */
/* Last edited on 2023-11-26 06:43:23 by stolfi */

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
#include <float_image_gradient_2.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

void float_image_gradient_sqr_relative_2
  ( float_image_t *A,
    int cA,
    double noise, 
    float_image_t *G,
    int cG    
  )
  {
    /* Get the image dimensions: */
    int NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);
    demand(cA < NC, "invalid A channel");
    
    /* Check size of output image: */
    demand((cG >= 0) && (cG < G->sz[0]), "invalid G channel");
    demand(G->sz[1] == NX, "image X size mismatch");
    demand(G->sz[2] == NY, "image Y size mismatch");
    
    /* Scratch images for gradient squared and variance: */
    float_image_t *G2 = float_image_new(1, NX, NY);
    float_image_t *V2 = float_image_new(1, NX, NY);
    
    /* Weight table: */
    int hw = 2, nw = 2*hw+1;
    double wt[nw];
    wt_table_binomial_fill(nw, wt, NULL);
    wt_table_normalize_sum(nw, wt);
 
    /* Process channel by channel: */
    int cMin = (cA >= 0 ? cA : 0);
    int cMax = (cA >= 0 ? cA : NC-1);
    double scale = 1.0/(cMax - cMin + 1.0);
    double noise2 = noise*noise;
    int c;
    for (c = cMin; c <= cMax; c++)
      { /* Compute the gradient squared and variance for channel {c} of image {A}: */
        float_image_gradient_sqr_sobel(A, c, G2, 0);
        float_image_local_avg_var(A, c, hw, wt, NULL, 0, V2, 0);
        /* Divide one by the other, store/add in {G}: */
        int x, y;
        for (x = 0; x < NX; x++)
          { for (y = 0; y < NY; y++)
              { double g2 = float_image_get_sample(G2, 0, x, y);
                double v2 = float_image_get_sample(V2, 0, x, y);
                double q = g2/(v2 + noise2)*scale;
                float *rP = float_image_get_sample_address(G, cG, x, y);
                if (c == cMin)
                  { (*rP) = (float)q; }
                else
                  { (*rP) += (float)q; }
              }
          }
      }
        
    float_image_free(G2); G2 = NULL;
    float_image_free(V2); V2 = NULL;
  }
  
