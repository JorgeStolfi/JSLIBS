/* Focus detector for multi-focus stereo. */
/* Last edited on 2023-01-29 14:00:53 by stolfi */

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

double multifok_score_from_terms(int32_t NT, double wt[], double term[]);
  /* Computes the sharpness score as the weighted sum 
    {SUM{wt[k]*term[k] : k \in 0..NT-1 }}. */ 

void multifok_score_read_term_names_and_weights
  ( FILE *rd, 
    int32_t NB, 
    char *belName[],
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P,
    int32_t *NT_P, 
    double **wt_P, 
    char ***termName_P,
    bool_t verbose
  );
  /* Reads from {rd} a a certain number {NT} of terms --
    quadratic local operators -- and the weights of each term in a sharpness score.
    
    Each line must have fields "{kt} {wtk} {namek}" where {kt} is the term index
    in {0..NT-1}, {wtk} is the weight of that term in the score formula, and
    {namek} is the formula of the term, as a sum of products of pairs of 
    basis element names, such as "DX*DX+DY*DY". See {multifok_term_indices_from_names} 
    for the valid format of {namek} and semantics.
    
    The procedure saves {wtk} and {namek} in arrays {wt[0..NT-1]} and
    {termName[0..NT-1]}. It also analyzes the {namek}s, determining the number {NP} of coeff pair 
    products that appear in them and building a table {prix[0..NP-1} as required by
    {multifok_term_values_from_basis_coeffs}.
    
    The values of {NP,prix,NT,wt,termName} are returned in
    {*NP_P,*prix_P,*NT_P,*termName_P}.
    
    Comments starting with "#" and blank lines are ignored. If {verbose}
    is true, also prints the data to stderr. */

void multifok_score_write_term_names_and_weights(FILE *wr, int32_t NT, char *termName[], double wt[]);
  /* Writes to {wr} the weights {wt[0..NT-1]} and names {termName[0..NT-1]} of the terms of a quadratic score
    formula, one per line, in the format described under {multifok_score_read_term_names_and_weights}. */

#endif
