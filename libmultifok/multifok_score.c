/* See {multifok_score.h}. */
/* Last edited on 2024-10-10 15:41:33 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <affirm.h>
#include <fget.h>
#include <bool.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>

#include <multifok_score.h>

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
  );
  /* Reads and parses a table of {NT} quadratic terms formed from {NP}
    products if pairs of {NB} coefficients of a linear window operator
    table. 
    
    Equivalent to {multifok_term_read_set_and_weights} to get
    the weights {wt[0..NT-1]} and formulas {termName[0..NT-1]}, followed
    by {multifok_term_indices_from_names} to obtain the product table
    {prix[0..NP-1]}.  Assumes that {belName[0..NB-1]} are the 
    names of the basis elements. */

double multifok_score_from_terms(int32_t NT, double wt[], double term[])
  { double score = 0.0;
    for (uint32_t kt = 0;  kt < NT; kt++)
      { /* Compute and combine term values: */
        score += wt[kt]*term[kt];
      }
    return score;
  }

#define multifok_score_C_COPYRIGHT \
    "© 2022 by the State University of Campinas (UNICAMP)"

