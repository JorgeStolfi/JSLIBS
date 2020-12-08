/* See {nmsim_check.h} */
/* Last edited on 2019-03-21 12:46:02 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>

#include <nmsim_basic.h>

#include <nmsim_check.h>
 
void nmsim_check_int64_value(char *name, int64_t v, int64_t vmin, int64_t vmax)
  { if ((v < vmin) || (v > vmax))
      { fprintf(stderr, "** value of %s = %ld\n", name, v);
        fprintf(stderr, " is out of range {%ld .. %ld}\n", vmin, vmax);
        demand(FALSE, "aborted");
      }
  }

void nmsim_check_double_value(char *name, double v, double vmin, double vmax)
  { if ((isnan(v)) || (v < vmin) || (v > vmax))
      { fprintf(stderr, "** value of %s = %24.16e", name, v);
        fprintf(stderr, " is out of range [%24.16e .. %24.16e]\n", vmin, vmax);
        demand(FALSE, "aborted");
      }
  }
