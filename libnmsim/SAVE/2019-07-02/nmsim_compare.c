/* See {nmsim_basic.h} */
/* Last edited on 2019-05-28 14:27:36 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
 
#include <nmsim_compare.h>

void nmsim_compare_int64_param(char *name, int64_t vr, int64_t ve)
  { if (vr != ve)
      { fprintf(stderr, "** parameter %s mismatch", name);
        fprintf(stderr, " read = %ld  expected = %ld", vr, ve);
        demand(FALSE, "aborted");
     }
  }

void nmsim_compare_double_param
  ( char *name, 
    double vr, 
    double ve,
    double prec,
    bool_t special_0, 
    bool_t special_1
  )
  { bool_t ok;
    if ((isnan(vr)) || (isnan(ve)))
      { /* Both must be {NAN}: */
        ok = (isnan(vr) & isnan(ve));
      }
    else if ((fabs(vr) == INF) || (fabs(ve) == INF))
      { /* Require identity: */
        ok = (vr == ve);
      }
    else if (special_0 && ((vr == 0.0) || (ve == 0.0)))
      { /* Require identity: */
        ok = (vr == ve);
      }
    else if (special_1 && ((fabs(vr) == 1.0) || (fabs(ve) == 1.0)))
      { /* Require identity: */
        ok = (vr == ve);
      }
    else
      { /* Require approximate equality: */
        ok = (fabs(vr - ve) <= prec);
      }
    if (! ok)
      { fprintf(stderr, "** parameter %s mismatch", name);
        fprintf(stderr, " read = %23.15e  expected = %23.15e", vr, ve);
        fprintf(stderr, " diff = %.15f\n", vr - ve);
        demand(FALSE, "aborted");
     }
  }
