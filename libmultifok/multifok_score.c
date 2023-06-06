/* See {multifok_score.h}. */
/* Last edited on 2023-04-18 16:54:48 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <wt_table.h>
#include <affirm.h>
#include <fget.h>
#include <bool.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>

#include <multifok_score.h>

double multifok_score_from_terms(int32_t NT, double wt[], double term[])
  { double score = 0.0;
    for (int32_t kt = 0; kt < NT; kt++)
      { /* Compute and combine term values: */
        score += wt[kt]*term[kt];
      }
    return score;
  }

void multifok_score_read_term_weights_names_get_indices
  ( FILE *rd, 
    int32_t NB, 
    char *belName[],
    int32_t *NT_P, 
    double **wt_P, 
    char ***termName_P,
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P,
    bool_t verbose
  )
  {
    multifok_term_read_weights_and_names(rd, NT_P, wt_P, termName_P, verbose);
    multifok_term_indices_from_names(NB, belName, *NT_P, *termName_P, NP_P, prix_P, verbose);
  }

#define multifok_score_C_COPYRIGHT \
    "© 2022 by the State University of Campinas (UNICAMP)"

