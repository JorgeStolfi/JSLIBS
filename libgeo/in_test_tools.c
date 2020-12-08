/* See in_test_tools.h. */
/* Last edited on 2011-12-23 22:07:45 by stolfilocal */

#include <in_test_tools.h>

#include <affirm.h>
#include <stdio.h>
#include <stdint.h>
#include <stdint.h>
#include <jswsize.h>

void in_do_check_eq(int64_t x, int64_t y, int *i, int *j, char *msg, in_LOCPARMS)
  { if (x != y)
      { fprintf(stderr, " **"); 
        if (i != NULL) { fprintf(stderr, " [%d]", *i); }
        if (j != NULL) { fprintf(stderr, " [%d]", *j); }
        fprintf(stderr, (" %" int64_d_fmt " != %" int64_d_fmt "\n"), x, y);
        programerror(msg, file, line, func);
      }
  }
