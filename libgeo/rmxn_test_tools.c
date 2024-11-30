/* See rmxn_test_tools.h. */
/* Last edited on 2024-11-23 18:55:40 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <rn.h>
#include <rmxn.h>
#include <jsrandom.h>

#include <rmxn_test_tools.h>

#define NO NULL

void rmxn_test_tools_check_all_different(uint32_t m, uint32_t n, double *A, char *msg)
  { 
    for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { double Aij = A[i*n+j];
            /* Check that {A[i,j]} is different from all previous elements: */
            for (uint32_t i1 = 0;  i1 <= i; i1++)
              { for (uint32_t j1 = 0;  j1 < (i == i1 ? j : n); j1++)
                 { demand(A[i1*n+j1] != Aij, msg); }
              }
          }
      }
  }
