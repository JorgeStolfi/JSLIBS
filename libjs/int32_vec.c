/* See int32_vec.h. */
/* Last edited on 2025-01-22 19:05:29 by stolfi */

/* Copyright © 2003 Jorge Stolfi, Unicamp. See note at end of file. */
/* Based on Params.m3 by J.Stolfi, DEC-SRC, 1988.  */

#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>

#include <vec.h>
#include <affirm.h>

#include <int32_vec.h>

void int32_vec_regularize(int32_vec_t *v, int32_t NC, int32_t defval)
  { uint32_t KC = v->ne;
    if ((KC != NC) && (KC <= 1))
      { int32_t val = (KC == 1 ? v->e[0] : defval);
        /* Provide default value: */
        int32_vec_expand(v, (vec_index_t)NC);
        for (int32_t c = 0; c < NC; c++) { v->e[c] = val; }
        int32_vec_trim(v, (vec_size_t)NC);
      }
  }
