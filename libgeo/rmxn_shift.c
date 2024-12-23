/* See rmxn_shift.h. */
/* Last edited on 2024-11-23 18:54:49 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <cmp.h>

#include <rmxn_shift.h>

void rmxn_shift_rows(uint32_t m, uint32_t n, double A[], double v[], double M[])
  { for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++) 
          { uint32_t ij = i*n + j; M[ij] = A[ij] + v[j]; }
      }
  }
  
void rmxn_shift_cols(uint32_t m, uint32_t n, double v[], double A[], double M[])
  { for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++) 
          { uint32_t ij = i*n + j; M[ij] = A[ij] + v[i]; }
      }
  }

