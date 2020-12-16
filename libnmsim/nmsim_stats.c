/* See {nmsim_stats.h} */
/* Last edited on 2020-12-16 00:17:41 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>

#include <nmsim_stats.h>

void nmsim_stats_print
  ( FILE *wr, 
    char *name, 
    nmsim_stats_t *S, 
    double prec, 
    bool_t sgn, 
    bool_t fudge_0, 
    bool_t fudge_1
  )
  {
    auto void wrd(double x);
      /* Calls {nmsim_write_double_value} with {x} and other parameters. */
    
    fprintf(wr, "%s:", name);
    fprintf(wr, " nvalid = %ld", S->nvs);
    fprintf(wr, " avg = "); wrd(S->avg);
    fprintf(wr, " dev = "); wrd(S->dev);
    fprintf(wr, " range = [ "); wrd(S->min); fputs(" _ ", wr); wrd(S->max); fputs(" ]", wr);
    fputc('\n', wr);
    return;
    
    void wrd(double x)
      { nmsim_write_double_value(wr, x, prec, sgn, fudge_0, fudge_1); }
  }

    
void  nmsim_stats_initialize(nmsim_stats_t *S)
  { S->nvs = 0;
    S->min = +INF; S->max = -INF;
    S->avg = 0.0; S->dev = 0.0;
  }
  
void  nmsim_stats_accumulate(nmsim_stats_t *S, double v)
  { S->nvs++;
    S->min = fmin(S->min, v);
    S->max = fmax(S->max, v);
    S->avg += v;
    S->dev += v*v;
  }
  
void  nmsim_stats_finalize(nmsim_stats_t *S)
  { 
    double dn = ((double)S->nvs);
    S->avg = (S->nvs < 1 ? 0.0 : S->avg/dn);
    S->dev = (S->nvs < 2 ? 0.0 : sqrt(fmax(0.0, S->dev - dn*S->avg*S->avg)/(dn - 1.0)));
  }
