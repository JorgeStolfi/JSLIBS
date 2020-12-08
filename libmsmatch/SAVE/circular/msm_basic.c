/* See {msm_basic.h} */
/* Last edited on 2008-05-25 03:25:53 by stolfi */

#define msm_basic_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>

#include <bool.h>
#include <sign.h>
#include <jsrandom.h>
#include <jsfile.h>

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>

void msm_check_malloc(char *tag)
  { 
    fprintf(stderr, "[%s]", tag);
    int N = 10; 
    char *p[N];
    int k;
    for (k = 0; k < N; k++) { p[k] = malloc(2000*(k+1)); }
    for (k = 0; k < N; k++) { free(p[k]); }
    for (k = 0; k < N; k++) { p[k] = malloc(2000*(k+1)); }
    for (k = 0; k < N; k++) { free(p[k]); }
  }

int msm_choose(int lo, int hi, double dProb)
  { int n = hi - lo;
    /* Slow and stupid solution: we simulate a random walk on the
      discrete grid {{lo..hi}×{lo..hi}}, from the lower corner to a
      diagonal goal line {(hi-k,lo+k)}. The Y position is the
      answer. */
    int x = 0, y = 0; /* Displacements from the lower corner: */
    while (x + y < n)
      { double toss = drandom();
        if ((toss -= dProb/2) < 0)
          { x++; }
        else if ((toss -= dProb/2) < 0)
          { y++; }
        else
          { x++; y++; }
      }
    /* If we went through the diagonal goal line, round Y randomly: */
    if (x + y > n)
      { assert(y > 0);
        if (drandom() < 0.5) { y--; }
      }
    return lo + y;
  }

FILE *msm_open_read(char *name, char *tag, char *ext, bool_t verbose)
  { char *fileName = NULL;
    asprintf(&fileName, "%s%s%s", name, tag, ext);
    FILE *rd = open_read(fileName, verbose);
    free(fileName);
    return rd;
  }

FILE *msm_open_write(char *name, char *tag, char *ext, bool_t verbose)
  { char *fileName = NULL;
    asprintf(&fileName, "%s%s%s", name, tag, ext);
    FILE *wr = open_write(fileName, verbose);
    free(fileName);
    return wr;
  }
