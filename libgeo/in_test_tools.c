/* See in_test_tools.h. */
/* Last edited on 2024-11-20 12:36:52 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <jswsize.h>
#include <affirm.h>

#include <in_test_tools.h>

void in_do_check_eq(int64_t x, int64_t y, uint32_t *i, uint32_t *j, char *msg, in_LOCPARMS)
  { if (x != y)
      { fprintf(stderr, " **"); 
        if (i != NULL) { fprintf(stderr, " [%d]", *i); }
        if (j != NULL) { fprintf(stderr, " [%d]", *j); }
        fprintf(stderr, (" %" int64_d_fmt " != %" int64_d_fmt "\n"), x, y);
        programerror(msg, file, line, func);
      }
  }
