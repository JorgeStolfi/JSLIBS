/* See ppv_image.h */
/* Last edited on 2021-07-09 01:13:29 by jstolfi */

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

uint16_image_t *ppv_image_from_array(ppv_array_t *A)
  {
    demand(A->maxsmp <= uint16_image_MAX_SAMPLE, "invalid maxsmp");
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
    J->maxval = (uint16_t)A->maxsmp;
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
                    ppv_sample_t smp = ppv_get_sample(A, ix);
                    row[x*chns + c] = (uint16_t)smp;
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
    uint16_t maxsmp = J->maxval; 
    assert(maxsmp <= ppv_MAX_SAMPLE_VAL);
    
    ppv_size_t sz[d];
    sz[0] = cols;
    sz[1] = rows;
    if (d == 3) { sz[2] = chns; }

    ppv_array_t *A = ppv_array_new(d, sz, maxsmp);
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
                    uint16_t smp = row[x*chns + c];
                    demand(smp <= maxsmp, "invalid sample in image");
                    ppv_set_sample(A, ix, (ppv_sample_t)smp);
                  }
              }
          }
      }
    return A;
  }
