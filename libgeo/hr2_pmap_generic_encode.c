/* See {hr2_pmap_generic_encode.h}. */
/* Last edited on 2024-09-17 04:57:32 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <assert.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <affirm.h>

#include <hr2_pmap_generic_encode.h>

void hr2_pmap_generic_encode(hr2_pmap_t *M, double y[])
  { int32_t k = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { y[k] = M->dir.c[i][j]; k++; }
      }
    assert(k == 9);
    double det = r3x3_det(&(M->dir));
    if (det < 0) 
      { for (int32_t j = 0; j < 3; j++) 
          { double tmp = y[3+j]; y[3+j] = y[6+j]; y[6 + j] = tmp; }
      }
  }

void hr2_pmap_generic_decode(double y[], hr2_pmap_t *M)
  { r3x3_t *A = &(M->dir);
    int32_t k = 0;
    for (int32_t i = 0; i < 3; i++)
      { for (int32_t j = 0; j < 3; j++)
          { A->c[i][j] = y[k]; k++; }
      }
    assert(k == 9);
    double det = r3x3_det(&(M->dir));
    if (det < 0)
      { for (int32_t j = 0; j < 3; j++)
          { double tmp = A->c[1][j]; A->c[1][j] = A->c[2][j];  A->c[2][j] = tmp; }
      }

    r3x3_inv(&(M->dir), &(M->inv));
  }
