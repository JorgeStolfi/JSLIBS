/* See {hr2_pmap_encode.h}. */
/* Last edited on 2024-09-17 15:00:44 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
 
#include <bool.h>
#include <sign.h>
#include <hr2.h>
#include <r3x3.h>
#include <affirm.h>

#include <hr2_pmap_translation_encode.h>
#include <hr2_pmap_congruence_encode.h>
#include <hr2_pmap_similarity_encode.h>
#include <hr2_pmap_affine_encode.h>
#include <hr2_pmap_generic_encode.h>

#include <hr2_pmap_encode.h>

int32_t hr2_pmap_encode_num_parameters(hr2_pmap_type_t type)
  {
    switch(type)
      { 
        case hr2_pmap_type_IDENTITY:
          return 0;
        case hr2_pmap_type_TRANSLATION:
          return 2;
        case hr2_pmap_type_CONGRUENCE:
          return 3;
        case hr2_pmap_type_SIMILARITY:
          return 4;
        case hr2_pmap_type_AFFINE:
          return 6;
        case hr2_pmap_type_GENERIC:
          return 9;
        case hr2_pmap_type_NONE:
          demand(FALSE, "invalid projective map type");
        default:
          demand(FALSE, "unimplemented map type");
      }
  }

void hr2_pmap_encode(hr2_pmap_t *M, hr2_pmap_type_t type, int32_t ny, double y[])
  {
    /* Assumes that {M} is a positive matrix with {M[0,0] = 1}, near the identity. */
    switch(type)
      { 
        case hr2_pmap_type_IDENTITY:
          { assert(ny == 0);
            /* Nothing to do: */
            break;
          }
        case hr2_pmap_type_TRANSLATION:
          { assert(ny == 2);
            hr2_pmap_translation_encode(M, y);
            break;
          }
        case hr2_pmap_type_CONGRUENCE:
          { assert(ny == 3);
            hr2_pmap_congruence_encode(M, y);
            break;
          }
        case hr2_pmap_type_SIMILARITY:
          { assert(ny == 4);
            hr2_pmap_similarity_encode(M, y);
            break;
          }
        case hr2_pmap_type_AFFINE:
          { assert(ny == 6);
            hr2_pmap_affine_encode(M, y);
            break;
          }
        case hr2_pmap_type_GENERIC:
          { assert(ny == 9);
            hr2_pmap_generic_encode(M, y);
            break;
          }
        case hr2_pmap_type_NONE:
          demand(FALSE, "invalid projective map type");
        default:
          demand(FALSE, "unimplemented map type");
      }
  }

void hr2_pmap_decode(int32_t ny, double y[], hr2_pmap_type_t type, sign_t sgn, hr2_pmap_t *M)
  { 
    /* Assumes that the non-variable elements of {M} are those of the identity matrix. */
    demand((sgn == +1) || (sgn == -1), "invalid {sgn}");
    switch(type)
      { 
        case hr2_pmap_type_IDENTITY:
          { assert(ny == 0);
            (*M) = hr2_pmap_identity();
            break;
          }
        case hr2_pmap_type_TRANSLATION:
          { assert(ny == 2);
            hr2_pmap_translation_decode(y, M);
            break;
          }
        case hr2_pmap_type_CONGRUENCE:
          { assert(ny == 3);
            hr2_pmap_congruence_decode(y, M);
            break;
          }
        case hr2_pmap_type_SIMILARITY:
          { assert(ny == 4);
            hr2_pmap_similarity_decode(y, M);
            break;
          }
        case hr2_pmap_type_AFFINE:
          { assert(ny == 6);
            hr2_pmap_affine_decode(y, M);
            break;
          }
        case hr2_pmap_type_GENERIC:
          { assert(ny == 9);
            hr2_pmap_generic_decode(y, M);
            break;
          }
        case hr2_pmap_type_NONE:
          demand(FALSE, "invalid projective map type");
        default:
          demand(FALSE, "unimplemented map type");
      }
    if (sgn == -1)
      { /* Combine with a {XY} swap: */
        hr2_pmap_t S = hr2_pmap_xy_swap();
        (*M) = hr2_pmap_compose(&S, M);
      }
  }
