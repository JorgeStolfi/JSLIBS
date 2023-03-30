/* See uint16_image_check_dither.h */
/* Last edited on 2023-03-17 21:14:51 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <jspnm.h>

#include <bool.h>
#include <affirm.h>
#include <uint16_image.h>
#include <uint16_image_check_dither.h>

bool_t uint16_image_check_dither(uint16_image_t *img, bool_t die)
  { 
    int32_t nc = img->chns;
    int32_t nx = img->cols;
    int32_t ny = img->rows;
    assert(nx > 0);
    assert(ny > 0);
    if (nx > (uint16_image_MAX_SAMPLE + 1)/ny)
      { fail_test(die,"dither matrix too big"); }
    uint32_t maxval = img->maxval;
    if (maxval != nx*ny-1)
      { fail_test(die, "bad maxval in dither matrix"); }
    bool_t seen[maxval+1]; /* {seen[s]} is true if sample value {s} occurs. */
    for (int32_t c = 0; c < nc; c++)
      { for (int32_t s = 0; s <= maxval; s++) { seen[s] = FALSE; }
        for (int32_t iy = 0; iy < ny; iy++)
          { uint16_t *smpy = img->smp[iy];
            for (int32_t ix = 0; ix < nx; ix++)
              { uint16_t s = smpy[ix*nc + c];
                assert(s <= maxval);
                if (seen[s])  { fail_test(die,"repeated sample value in dither matrix"); }
                seen[s] = TRUE;
              }
          }
      }
    return TRUE;
  }
