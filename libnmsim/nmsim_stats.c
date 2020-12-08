/* See {nmsim_stats.h} */
/* Last edited on 2020-12-07 16:19:09 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

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
