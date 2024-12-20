/* See {float_image_mmorph.h}. */
/* Last edited on 2024-12-04 23:30:02 by stolfi */

#include <math.h>
#include <assert.h>
#include <values.h>
#include <string.h>

#include <bool.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_mmorph.h>

/* INTERNAL PROTOTYPES */

/* IMPLEMENTATIONS */

float_image_t *float_image_mmorph_dilate(float_image_t *A, int32_t hw, double wt[])
  {
    /* Get the image dimensions: */
    int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);
    
    /* Allocate the dilated image: */
    float_image_t *G = float_image_new(NC, NX, NY);
    
    /* Fill it: */
    int32_t x, y;
    for (x = 0; x < NX; x++)
      { for (y = 0; y < NY; y++)
          { int32_t c;
            for (c = 0; c < NC; c++)
              { /* Compute the local gradient squared {g2} in channel {c}: */
                /* Compute the horizontal and vertical derivatives in channel {c}: */
                double fd = float_image_get_dilated(A, c, x, y, hw, wt);
                /* Store the dilation: */
                float_image_set_sample(G, c, x, y, (float)fd);
              }
          }
      }
    return G;
  }
