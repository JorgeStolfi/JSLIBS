/* See {frgb_Y_slice.h}. */
/* Last edited on 2023-03-07 01:35:28 by stolfi */ 

#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <frgb.h>
#include <jsmath.h>
#include <frgb_ops.h>
#include <frgb_Y_slice.h>

#define AUTHOR \
  "  Created mar/2023 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP)."

#define YR frgb_YR
#define YG frgb_YG
#define YB frgb_YB
  /* Luminance weights according to European TV standard. */

void frgb_Y_slice_corners(double z, int32_t *nf_P, frgb_t f[])
  { 
    /* Check {z} range: */
    demand((z >= 0) & (z <= 1.0), "invalid brightness value");
    
    /* Intersect the plane {Y = z} with the edges of the RGB unit cube: */
    int32_t nf = 0;
    if (z < 0) 
      { /* No intersection. */ }
    else if (z == 0)
      { /* Just the black (K) corner: */
        f[nf] = (frgb_t){{ 0, 0, 0 }}; nf++;
      }
    else if (z <= YB) /* 0.0~0.1 */
      { /* Triangle; edges are K-R,K-G,K-B: */
        f[nf] = (frgb_t){{ (float)(z/YR), 0, 0 }}; nf++;
        f[nf] = (frgb_t){{ 0, (float)(z/YG), 0 }}; nf++;
        f[nf] = (frgb_t){{ 0, 0, (float)(z/YB) }}; nf++;
      }
    else if (z <= YR) /* 0.1~0.3 */
      { /* Quadrilateral; edges are K-R,K-G,B-C,B-M: */
        f[nf] = (frgb_t){{ (float)(z/YR), 0, 0 }};      nf++;
        f[nf] = (frgb_t){{ 0, (float)(z/YG), 0 }};      nf++;
        f[nf] = (frgb_t){{ 0, (float)((z-YB)/YG), 1 }}; nf++;
        f[nf] = (frgb_t){{ (float)((z-YB)/YR), 0, 1 }}; nf++;
      }
    else if (z < YR+YB) /* 0.3~0.4 */
      { /* Pentagon; edges are R-Y,K-G,B-C,B-M,R-M: */
        f[nf] = (frgb_t){{ 1, (float)((z-YR)/YG), 0 }}; nf++;
        f[nf] = (frgb_t){{ 0, (float)(z/YG), 0 }};      nf++;
        f[nf] = (frgb_t){{ 0, (float)((z-YB)/YG), 1 }}; nf++;
        f[nf] = (frgb_t){{ (float)((z-YB)/YR), 0, 1 }}; nf++;
        f[nf] = (frgb_t){{ 1, 0, (float)((z-YR)/YB) }}; nf++;
      }
    else if (z <= YG) /* 0.4~0.6 */
      { /* Quadrilateral; edges are R-Y,K-G,B-C,M-W: */
        f[nf] = (frgb_t){{ 1, (float)((z-YR)/YG), 0 }};    nf++;
        f[nf] = (frgb_t){{ 0, (float)(z/YG), 0 }};         nf++;
        f[nf] = (frgb_t){{ 0, (float)((z-YB)/YG), 1 }};    nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR-YB)/YG), 1 }}; nf++;
      }
    else if (z < YG+YB) /* 0.6~0.7 */
      { /* Pentago; edges are R-Y,G-Y,G-C,B-C,M-W: */
        f[nf] = (frgb_t){{ 1, (float)((z-YR)/YG), 0 }};    nf++;
        f[nf] = (frgb_t){{ (float)((z-YG)/YR), 1, 0 }};    nf++;
        f[nf] = (frgb_t){{ 0, 1, (float)((z-YG)/YB) }};    nf++;
        f[nf] = (frgb_t){{ 0, (float)((z-YB)/YG), 1 }};    nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR-YB)/YG), 1 }}; nf++;
      }
    else if (z < YR+YG) /* 0.7~0.9 */
      { /* Quadrilateral; edges are R-Y,G-Y,C-W,M-W: */
        f[nf] = (frgb_t){{ (float)((z-YG)/YR), 1, 0 }};    nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR)/YG), 0 }};    nf++;
        f[nf] = (frgb_t){{ (float)((z-YG-YB)/YR), 1, 1 }}; nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR-YB)/YG), 1 }}; nf++;
      }
    else if (z < 1.0) /* 0.9~1.0 */
      { /* Triangle; edges are Y-W,C-W,M-W: */
        f[nf] = (frgb_t){{ (float)((z-YG-YB)/YR), 1, 1 }}; nf++;
        f[nf] = (frgb_t){{ 1, (float)((z-YR-YB)/YG), 1 }}; nf++;
        f[nf] = (frgb_t){{ 1, 1, (float)((z-YG-YR)/YB) }}; nf++;
      }
    else if (z == 1.0)
      { /* Single point, the white (W) corner: */
        f[nf] = (frgb_t){{ 1, 1, 1 }}; nf++;
      }
    else
      { /* No intersection. */ }
      
    if (nf >= 2)
      { /* Eliminate consecutive corners that are repeated by roundoff: */
        int32_t nf_new = 0; /* Number of distinct corners. */
        f[nf_new] = f[0]; nf_new++;
        for (uint32_t k = 1;  k < nf; k++)
          { if (! frgb_eq(&(f[k]), &(f[nf_new-1])))
              { f[nf_new] = f[k]; nf_new++; }
          }
        if (nf_new == 2)
          { /* Must be distinct only by roundoff: */
            nf_new = 1;
          }
        nf = nf_new;
      }
        
    (*nf_P) = nf;
  }
