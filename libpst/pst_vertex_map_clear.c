/* See pst_vertex_map_clear.h */
/* Last edited on 2025-03-06 13:03:15 by stolfi */

#include <math.h>
#include <assert.h>
#include <stdint.h>

#include <bool.h>
#include <i2.h>
#include <vec.h>
#include <affirm.h>
#include <float_image.h>

#include <pst_basic.h>
#include <pst_map_clear.h>

#include <pst_vertex_map_clear.h>

void pst_vertex_map_clear
  ( float_image_t *A,
    int32_t wch,
    i2_vec_t *imin,
    i2_vec_t *imax
  )
  { for (int32_t k = 0; k < imin->ne; k++)
      { i2_t *imink = &(imin->e[k]);
        i2_t *imaxk = &(imax->e[k]);
        pst_map_clear(A, wch, imink->c[0], imaxk->c[0]+1, imink->c[1], imaxk->c[1]+1);
      }
  }
