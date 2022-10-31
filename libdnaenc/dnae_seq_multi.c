/* See dnae_seq_multi.h */
/* Last edited on 2022-10-31 14:33:00 by stolfi */

#define dnae_seq_multi_C_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)" \
  " and the Federal Fluminense University (UFF)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include <fftw3.h>

#include <vec.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <wt_table.h>
#include <jsfile.h>
#include <hermite3.h>
#include <conv_filter.h>

#include <msm_seq_desc.h>

#include <dnae_nucleic.h>
#include <dnae_sample.h>
#include <dnae_datum.h>

#include <dnae_seq.h>
#include <dnae_seq_multi.h>

/* IMPLEMENTATIONS */

void dnae_seq_multi_free_datums(dnae_seq_t seq[], int32_t maxLevel)
  { int32_t level;
    for (level = 0; level <= maxLevel; level++)
      { dnae_seq_free_datums(&(seq[level])); }
  }

void dnae_seq_multi_free(dnae_seq_t *seq[], int32_t maxLevel)
  { int32_t level;
    for (level = 0; level <= maxLevel; level++)
      { dnae_seq_free(seq[level]); }
  }
   
void dnae_seq_multi_filter
  ( dnae_seq_t *s, 
    int32_t maxLevel, 
    double_vec_t *wt0,
    char *wname0,
    double_vec_t *wt1,
    char *wname1,
    int8_t ek0,
    dnae_seq_t sr[]
  )
  { int32_t k; /* Level of hierarchy. */
    for (k = 0; k <= maxLevel; k++)
      { if (k == 0)
          { sr[k] = dnae_seq_copy(s); }
        else
          { double_vec_t *wt = (k == 1 ? wt0 : wt1);
            char *wname = (k == 1 ? wname0 : wname1);
            int8_t ekk = (int8_t)(k == 1 ? ek0 : 1);    /* Exponent of step for upsampling/downsampling. */
            int8_t ekf = (int8_t)(ekk >= 0 ? ekk : 0);  /* Exponent of step for downsampling. */
            sr[k] = dnae_seq_filter(&(sr[k-1]), wt, ekf, wname);
            if (ekk < 0) 
              { dnae_seq_t tmp = sr[k];
                sr[k] = dnae_seq_interpolate(&tmp, ekk);
                dnae_seq_free_datums(&tmp);
              }
          }
      }
  }
 
void dnae_seq_multi_get_2017_paper_weights
  ( double_vec_t *wt0_P, 
    char **wname0_P, 
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  )
  { 
    int32_t r0 = 6, n0 = 2*r0+1;
    double_vec_t wt0 = double_vec_new(n0);
    wt0.e[r0+0] = 9992;
    wt0.e[r0+1] = 7786;
    wt0.e[r0+2] = 3680;
    wt0.e[r0+3] = 1055;
    wt0.e[r0+4] =  183;
    wt0.e[r0+5] =   19;
    wt0.e[r0+6] =    1;
    for (int32_t i = 1; i <= r0; i++) { wt0.e[r0-i] = wt0.e[r0+i]; }
    char *wname0 = wt_table_make_descr(n0, wt0.e, "%d");
    wt_table_normalize_sum(n0, wt0.e); /* Must be after {wt_table_make_descr}. */
    
    int32_t r1 = 10, n1 = 2*r1+1;
    double_vec_t wt1 = double_vec_new(n1);
    wt1.e[r1+ 0] = 9992;
    wt1.e[r1+ 1] = 9193;
    wt1.e[r1+ 2] = 7161;
    wt1.e[r1+ 3] = 4722;
    wt1.e[r1+ 4] = 2636;
    wt1.e[r1+ 5] = 1245;
    wt1.e[r1+ 6] =  498;
    wt1.e[r1+ 7] =  169;
    wt1.e[r1+ 8] =   48;
    wt1.e[r1+ 9] =   12;
    wt1.e[r1+10] =    2;
    for (int32_t i = 1; i <= r1; i++) { wt1.e[r1-i] = wt1.e[r1+i]; }
    char *wname1 = wt_table_make_descr(n1, wt1.e, "%d");
    wt_table_normalize_sum(n1, wt1.e); /* Must be after {wt_table_make_descr}. */
    
    (*wt0_P) = wt0; (*wname0_P) = wname0;
    (*wt1_P) = wt1; (*wname1_P) = wname1;
  }    
