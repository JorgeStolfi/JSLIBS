#ifndef nmsim_stats_H
#define nmsim_stats_H
 
/* Statistical summaries for neuromat network simulation. */
/* Last edited on 2020-12-15 21:37:19 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>

typedef struct nmsim_stats_t
  { int64_t nvs;       /* Number of valid samples. */
    double min;        /* Min valid sample value. */
    double max;        /* Max valid sample value. */               
    double avg;        /* Average of valid sample values. */       
    double dev;        /* Deviation of valid sample values. */   
  } nmsim_stats_t;
  /* Statistical summary of a set of samples.  Samples
    are valid if they are neither {±INF} nor {NAN}. */
    
void nmsim_stats_print
  ( FILE *wr, 
    char *name, 
    nmsim_stats_t *S, 
    double prec, 
    bool_t sgn, 
    bool_t fudge_0, 
    bool_t fudge_1
  );
  /* Writes {S} to {wr}, in one line, prefixed by {name}, in human-readable format.
    The parameters {prec,sgn,fudge_0,fudge_1} as as in {nmsim_write_double_value}. */

/* INCREMENTAL STATS GATHERING 

  An {nmsim_stats_t} record can be computed incrementally by calling
  {nmsim_stats_initialize} just once, then {nmsim_stats_accumulate} zero or
  more times (once for each sample), then {nmsim_stats_finalize} just
  once. Between initialization and finalization, the fields {.avg} and {.dev}
  are temporarily used to hold the sum of samples and the sum of their
  squares, respetively. */
    
void  nmsim_stats_initialize(nmsim_stats_t *S);
  /* Clears {S} preparing it for {nmsim_stats_accumulate}. Namely, sets
    {S.min} to {+INF}, {S.max} to {-INF}, and all other fields to
    zero. */
  
void  nmsim_stats_accumulate(nmsim_stats_t *S, double v);
  /* Updates {S.min,S.max} to include the sameple {v}.  Increments {nvs},
    adds {v} to {avg}, and its square to {dev}. */
  
void  nmsim_stats_finalize(nmsim_stats_t *S, double v);
  /* Converts {S.avg} from sum of sample values to their average.
    Converts {S.dev} from sum of squared samples to their standard
    deviation.
    
    If there were no samples, sets {S.avg,S.dev} to zero (instead of
    {NAN}) and leaves {S.min = +INF}, {S.max = -INF}. If there was only
    one sample, sets {S.dev} to zero (instead of {NAN}). */

#endif 
