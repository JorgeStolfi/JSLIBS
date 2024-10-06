/* See {tosl_slice.h} */
/* Last edited on 2024-10-04 10:19:40 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>

#include <haf.h>

#include <tosl.h>

#include <tosl_slice.h>

tosl_slice_t *tosl_slice_new(int32_t NV_max, tosl_coord_t Zp)
  { tosl_slice_t *S = malloc(sizeof(tosl_slice_t));
    S->iarc = malloc(NV_max*sizeof(tosl_arc_id_t));
    S->NV_max = NV_max;
    S->NV = 0;
    S->Z = Zp;
    return S;
  }

void tosl_slice_free(tosl_slice_t *S)
  { free(S->iarc);
    free(S);
  }
