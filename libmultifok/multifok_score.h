/* Focus detector for multi-focus stereo. */
/* Last edited on 2024-08-02 15:57:28 by stolfi */

#ifndef multifok_score_H
#define multifok_score_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>

/* 
  FOCUS DETECTOR OPERATORS
  
  A /sharpness score/ is a local image operator that estimates the
  sharpness of the image at the applied pixel.  Ideally it is close to 0
  for a very unfocused (blurred) image, and close to 1 for a well-focused (sharp) one.
  
  In practice, it is sufficient that, when applied at the same point of 
  the same image blurred by various amounts, it has the maximum value for 
  the version with least blurring.
  
  The focus detector operators in this module work on a window of {NW}
  by {NW} samples, stored in an array {x[0..NS-1]} of {float}s, where {NS
  = NW*NW}, linearized by rows.  See {multifok_window.h} for 
  mreo details. */ 

typedef double multifok_score_op_t(int32_t NC, double wt[]);
  /* Given an array of window samples ?? */ 

double multifok_score_from_terms(int32_t NT, double wt[], double term[]);
  /* Computes the sharpness score as the weighted sum 
    {SUM{wt[k]*term[k] : k \in 0..NT-1 }}. */ 

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
    
    Equivalent to {multifok_term_read_weights_and_names} to get
    the weights {wt[0..NT-1]} and formulas {termName[0..NT-1]}, followed
    by {multifok_term_indices_from_names} to obtain the product table
    {prix[0..NP-1]}.  Assumes that {belName[0..NB-1]} are the 
    names of the basis elements. */

#endif
