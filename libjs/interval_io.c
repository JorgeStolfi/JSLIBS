/* See {interval_io.h} */
/* Last edited on 2024-11-15 19:13:13 by stolfi */

#include <stdio.h>
#include <math.h>

#include <interval.h>

#include <interval_io.h>

void interval_print(FILE *wr, interval_t *X)
  { interval_gen_print(wr, X, NULL, NULL, NULL, NULL); }

void interval_gen_print(FILE *wr, interval_t *X, char *fmt, char *lp, char *sep, char *rp)
  { if (fmt == NULL) { fmt = "%16.8e"; }
    if (lp == NULL) { lp = "["; }
    if (sep == NULL) { sep = " "; }
    if (rp == NULL) { rp = "]"; }
    fputs(lp, wr);
    fprintf(wr, fmt, X->end[0]);
    fputs(sep, wr);
    fprintf(wr, fmt, X->end[1]);
    fputs(rp, wr);
  }
