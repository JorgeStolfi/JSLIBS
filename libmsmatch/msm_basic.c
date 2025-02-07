/* See {msm_basic.h} */
/* Last edited on 2022-10-20 06:41:43 by stolfi */

#define msm_basic_C_COPYRIGHT \
  "Copyright � 2005  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <stdint.h>

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
    int32_t N = 10; 
    char *p[N];
    int32_t k;
    for (k = 0; k < N; k++) { p[k] = malloc(2000*(k+1)); }
    for (k = 0; k < N; k++) { free(p[k]); }
    for (k = 0; k < N; k++) { p[k] = malloc(2000*(k+1)); }
    for (k = 0; k < N; k++) { free(p[k]); }
  }

int32_t msm_choose(int32_t lo, int32_t hi, double dProb)
  { int32_t n = hi - lo;
    /* Slow and stupid solution: we simulate a random walk on the
      discrete grid {{lo..hi}�{lo..hi}}, from the {(lo,lo)} corner to a
      diagonal goal line {(hi-k,lo+k)}. The Y position is the
      answer. */
    int32_t x = 0, y = 0; /* Displacements from the lower corner: */
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

int64_t msm_round(double x)
  {
    int64_t k = (int64_t)floor(x);
    double f = x - (double)k;
    assert((f >= 0.0) && (f < 1.0));
    if (f < 0.5)
      { return k; }
    else if (f > 0.5)
      { return k+1; }
    else
      { return k + (k & 1); }
  }

FILE *msm_open_read(char *name, char *tag, char *ext, bool_t verbose)
  { char *fileName = NULL;
    char *fileName = jsprintf("%s%s%s", name, tag, ext);
    FILE *rd = open_read(fileName, verbose);
    free(fileName);
    return rd;
  }

FILE *msm_open_write(char *name, char *tag, char *ext, bool_t verbose)
  { char *fileName = NULL;
    char *fileName = jsprintf("%s%s%s", name, tag, ext);
    FILE *wr = open_write(fileName, verbose);
    free(fileName);
    return wr;
  }
