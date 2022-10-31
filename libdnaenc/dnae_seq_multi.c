/* See dnae_seq_multi.h */
/* Last edited on 2022-10-31 12:12:35 by stolfi */

#define dnae_seq_multi_C_COPYRIGHT \
  "Copyright � 2005  by the State University of Campinas (UNICAMP)" \
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
    double_vec_t *wtb0,
    char *wname0,
    double_vec_t *wtb1,
    char *wname1,
    int8_t ek0,
    dnae_seq_t sr[]
  )
  { int32_t k; /* Level of hierarchy. */
    for (k = 0; k <= maxLevel; k++)
      { if (k == 0)
          { sr[k] = dnae_seq_copy(s); }
        else
          { double_vec_t *wtb = (k == 1 ? wtb0 : wtb1);
            char *wname = (k == 1 ? wname0 : wname1);
            int8_t ekk = (int8_t)(k == 1 ? ek0 : 1);    /* Exponent of step for upsampling/downsampling. */
            int8_t ekf = (int8_t)(ekk >= 0 ? ekk : 0);  /* Exponent of step for downsampling. */
            sr[k] = dnae_seq_filter(&(sr[k-1]), wtb, ekf, wname);
            if (ekk < 0) 
              { dnae_seq_t tmp = sr[k];
                sr[k] = dnae_seq_interpolate(&tmp, ekk);
                dnae_seq_free_datums(&tmp);
              }
          }
      }
  }
 
