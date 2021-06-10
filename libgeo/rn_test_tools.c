/* See rn_test_tools.h. */
/* Last edited on 2021-06-09 20:34:54 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <affirm.h>

#include <rn_test_tools.h>

void rn_do_check_eq(double x, double y, int32_t *i, int32_t *j, char *msg, rn_LOCPARMS)
  { if (x != y)
      { fprintf(stderr, " **"); 
        if (i != NULL) { fprintf(stderr, " [%d]", *i); }
        if (j != NULL) { fprintf(stderr, " [%d]", *j); }
        fprintf(stderr, " %20.16e != %20.16e\n", x, y);
        programerror(msg, file, line, func);
      }
  }

void rn_do_check_eps(double x, double y, double eps, int32_t *i, int32_t *j, char *msg, rn_LOCPARMS)
  { if (fabs(x - y) > eps)
      { fprintf(stderr, " **"); 
        if (i != NULL) { fprintf(stderr, " [%d]", *i); }
        if (j != NULL) { fprintf(stderr, " [%d]", *j); }
        fprintf(stderr, " | %20.16e - %20.16e | > %20.16e\n", x, y, eps);
        programerror(msg, file, line, func);
      }
  }

