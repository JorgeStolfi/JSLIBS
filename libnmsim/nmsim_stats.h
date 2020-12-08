#ifndef nmsim_stats_H
#define nmsim_stats_H
 
/* Statistical summaries for neuromat network simulation. */
/* Last edited on 2020-12-07 16:08:49 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <nmsim_basic.h>

typedef struct nmsim_stats_t
  { nmsim_time_t nvs;  /* Number of valid samples. */
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

#endif 
