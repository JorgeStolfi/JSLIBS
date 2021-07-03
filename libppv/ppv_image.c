/* See ppv_image.h */
/* Last edited on 2021-07-03 14:58:07 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <uint16_image.h>

#include <ppv_types.h>
#include <ppv_array.h>

#include <ppv_image.h>

uint16_image_t *ppv_image_from_array(ppv_array_t *A, ppv_sample_t maxval)
  {
    demand(maxval <= uint16_image_MAX_SAMPLE, "invalid maxval");
    ppv_dim_t d = A->d;
    demand((d == 2) || (d == 3), "array must have 2 or 3 axes");
    
    demand(A->size[0] <= INT32_MAX, "too many columns"); 
    int32_t cols = (int32_t)A->size[0];
    
    demand(A->size[1] <= INT32_MAX, "too many rows"); 
    int32_t rows = (int32_t)A->size[1];
    
    int32_t chns ;
    if (A->d == 2)
      { chns = 1; }
    else
      { demand(A->size[2] <= INT32_MAX, "too many channels");
        chns = (int32_t)A->size[2];
      }
    
    uint16_image_t *J = uint16_image_new(cols, rows, chns);
    J->maxval = (uint16_t)maxval;
    if ((cols > 0) && (rows > 0) && (chns > 0)) 
      { /* Copy samples: */
        ppv_index_t ix[d];
        for (int32_t y = 0; y < rows; y++)
          { uint16_t *row = J->smp[y];
            ix[1] = y;
            for (int32_t x = 0; x < cols; x++)
              { ix[0] = x;
                for (int32_t c = 0; c < chns; c++)
                  { if (d == 3) { ix[2] = c; }
                    ppv_sample_t v = ppv_get_sample(A, ix);
                    demand(v <= maxval, "invalid sample in array");
                    row[x*chns + c] = (uint16_t)v;
                  }
              }
          }
      }
    return J;
  }

ppv_array_t *ppv_image_to_array(uint16_image_t *J)
  {
    ppv_dim_t d = (J -> chns == 1 ? 2 : 3);
    demand((d == 2) || (d == 3), "array must have 2 or 3 axes");
    
    int32_t cols = J->cols; demand(cols <= ppv_MAX_SIZE, "too many columns"); 
    int32_t rows = J->rows; demand(rows <= ppv_MAX_SIZE, "too many rows"); 
    int32_t chns = J->chns; demand(chns <= ppv_MAX_SIZE, "too many channels");
    uint16_t maxval = J->maxval; 
    assert(maxval <= ppv_MAX_SAMPLE_VAL);
    
    ppv_size_t sz[d];
    sz[0] = cols;
    sz[1] = rows;
    if (d == 3) { sz[2] = chns; }
    ppv_nbits_t bps = 0;
    while (maxval >= (1 << bps)) { bps++; }
    assert(bps <= ppv_MAX_BPS);
    ppv_nbits_t bpw = ppv_best_bpw(bps); 
    
    ppv_array_t *A = ppv_array_new(d, sz, bps, bpw);
    if ((cols > 0) && (rows > 0) && (chns > 0)) 
      { /* Copy samples: */
        ppv_index_t ix[d];
        for (int32_t y = 0; y < rows; y++)
          { uint16_t *row = J->smp[y];
            ix[1] = y;
            for (int32_t x = 0; x < cols; x++)
              { ix[0] = x;
                for (int32_t c = 0; c < chns; c++)
                  { if (d == 3) { ix[2] = c; }
                    uint16_t v = row[x*chns + c];
                    demand(v <= maxval, "invalid sample in image");
                    ppv_set_sample(A, ix, (ppv_sample_t)v);
                  }
              }
          }
      }
    return A;
  
  }
