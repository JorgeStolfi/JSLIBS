/* Quadratic operator terms for multi-focus stereo. */
/* Last edited on 2023-04-18 22:36:28 by stolfi */

#ifndef multifok_term_H
#define multifok_term_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <multifok_window.h>
#include <multifok_basis.h>

/* QUADRATIC TERMS FOR SHARPNESS ESTIMATION

  The {multifok} library uses local sharpness estimators
  that are quadratic functions of brightness- and contrast-normalized window
  samples.  
  
  INDEPENDENT QUADRATIC TERMS
  
  If the samples {s[0..NS-1]} are viewed as a row vector, the score formula
  is 
  
    {score(s) = wt0 + \SUM{ t \in 0..NT-1 : wt[t]*s*B*Q[t]*B'*s' }}
    
  where {wt0} and {wt[0..NT-1]} are scalar coefficients, {B} is an {NS×NB} basis
  projection matrix, and each {Q[t]} is a lower triangular {NB×NB} matrix.  Each 
  expression {s*B*Q[t]*B'*s'} is a /term/ of the score formula.  See 
  "papers/graphics/multifocus/2023-sharpness-score/" for the 
  rationale for this decomposition. 
  
  Currently the elements of each term matrix {Q[t]} are supposed to be
  either 1 or 0. Therefore each term is a sum of products of pairs of
  basis coefficients {b[jb1]*b[jb2]}, where {b[0..NB-1]} is the coeff
  vector {b = s*B}, and {0 <= jb1 <= jb2 < NB}.
  
  Moreover we assume that no two matrices {Q[t']}, {Q[t'']} have "1" in
  the same position. Which means that each product of the above form
  appears in at most one term. Which means that the total number of
  products in all terms is at most {NB*(NB+1)/2}.
  
  !!! Is it worth lifting the above restrictions? !!!
  
  */

typedef struct multifok_term_prod_t
  { int32_t jb1; /* Index of first basis elem coeff. */
    int32_t jb2; /* Index of first basis elem coeff. */
    int32_t kt;  /* Index of term which this product belongs to. */
    char *pname; /* Name (formula) of product, e.g. "DX*DX" of "Smm*Spp". */
  } multifok_term_prod_t;
  /* Record describing a product of two coeffs {coeff[jb1]*coeff[jb2]}
    in some local linear basis, that contributes to one indpendent term
    of a quadratic local sharpness estimator.
    
    The product must be in canonical form, meaning {jb1<=jb2}. The name
    shoud be "{el1}*{el2}" where {el1} and {el2} are the names of the
    correspoing basis elements.
    
    As a special case, {jb1} and {jb2} may be {-1}, representing the
    empty product whose value is constant 1. The {pname} then should be
    "1". */

void multifok_term_values_from_basis_coeffs
  ( int32_t NB, 
    double coeff[],
    int32_t NP,
    multifok_term_prod_t prix[],
    int32_t NT,
    double term[]
  );
  /* Assumes that {coeff[0..NB-1]} are the coefficients of the window
    samples in some local operator basis. Computes {NT} quadratic
    operators that are sums of products of pairs of coefficients, as
    specified by the table {prix[0..NT-1]}.
    
    Namely, for each entry {prix[k] = (jb1,jb2,kt,pname)}, computes the
    product {coeff[jb1]*coeff[jb2]} and adds it to {term[kt]}. 
    As a sepcial case, if {jb1=jb2=-1}, adds 1.0 to {term[kt]}. */ 

void multifok_term_write_names
  ( FILE *wr, 
    int32_t NT,
    char *termName[]
  );
  /* Writes to {wr} the the 
    term names (formulas) {termName[0..NT-1]}, one per line. */

void multifok_term_read_weights_and_names
  ( FILE *rd, 
    int32_t *NT_P,
    double **wt_P, 
    char ***termName_P,
    bool_t verbose
  );
  /* Reads from {rd} a table {termName[0..NT-1]} describing a certain
    number {NT} of terms of degree 2 consisting ofproducts of pairs of
    basis coefficients.
  
    The file has one line for each quadratic term. Each line has the
    format "{kt} {wt[kt]} {termName[kt]}" where {kt} is the term index in
    {0..NT-1}, {wt[kt]} is a non-negative fractional weight,
    and {termName[kt]} is the term's name (formula). The
    latter is a string, either just the string "1" or a sum of products
    of names of basis elements like "Smm*Spp+Spm+Spp". These names must
    be valid for {multifok_term_indices_from_names}.
    
    The number of terms {NT} is inferred from the number of lines in the
    file.
    
    The value of {NT} and the arrays {wt} and {termName} are returned in
    {*NT_P,*wt_P,*termName_P}. However, if {wt_P} is {NULL}, the weights are 
    discarded.
    
    Comments starting with "#" and blank lines are ignored. If {verbose}
    is true, also prints the data to stderr. */

void multifok_term_write_weights_and_names
  ( FILE *wr, 
    int32_t NT,
    double wt[],
    char *termName[]
  );
  /* Writes to {wr} the weights {wt[0..NT-1]} and the 
    term names (formulas) {termName[0..NT-1]}, one per line,
    in the format expected by {multifok_term_read_weights_and_names}. 
    
    If {wt} is {NULL}, assumes it is a vector of ones. */

void multifok_term_indices_from_names
  ( int32_t NB, 
    char *belName[],
    int32_t NT, 
    char *termName[], 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P, 
    bool_t verbose
  );
  /* The parameter {termName[0..NT-1]} must be a list of formulas of
    quadratic operator terms. Each formula must be a sum of products of
    pairs of coefficients, like "DX*DY+DXDX*DYDY". Each product in that
    sum must be "{el1}*{el2}" where {el1=belName[jb1]} and
    {el2=belName[jb2]} for some jb1,jb2} with {0 <= jb1,jb2 < NB}. As
    a special case, a product may be just "1".
  
    The procedure determines the total number {NP} of products appearing
    in the {termName[0..NT-1]}, such as "DX*DY" or "Soo*Sop", and
    creates a table {prix[0..NP-1]} that describes the products and the
    terms they belong to. See {multifok_term_values_from_basis_coeffs}
    for the meaning of {prix}.
    
    In each product, the two factors are internally swapped if needed so
    that {jb1 <= jb2}. After this adjustment, each product should appear
    only once. Thus {NP} is at most {NB*(NP+1)/2}. Also each term must
    have at least one product, so {NT} is at most equal to {NP}.
    
    The product count {NP} and the table {prix} are returned in {*NP_P}
    and {*prix_P}. The {pname} fields in the {prix} table will be newly
    allocated strings. */

/* PARSED TERM TABLE I/O */

void multifok_term_read_index_table
  ( FILE *rd, 
    int32_t NB,
    char *belName[],
    int32_t *NP_P,
    multifok_term_prod_t **prix_P,
    int32_t *NT_P,
    char ***termName_P,
    bool_t verbose
  );
  /* Reads from {rd} a table describing how a certain number {NT} of
    terms of degree 2 are to be composed from {NP} products of pairs of
    {NB} basis coefficients.
  
    The file has one line for each coeff product to be evaluated. Each
    line has the format "{ip} {jb1} {jb2} {kt} {pname}" where {ip} is the product
    index in {0..NP-1}, {jb1} and {jb2} are the indices of the two basis
    coeffs, with {0 <= jb1 <= jb2 < NB}; {kt} is the index of the term to which
    this product is part; and {pname} is the textual product
    "{el1}*{el2}", e.g. "Xom*Xpp"; where {el1} must be {belName[jb1]} and 
    {el2} must be {belName[jb2]}.
    
    The pairs {jb1,jb2} must be all distinct  so the
    maximum value of {NP} is {NB*(NB+1)/2}. The quadruples {(jb1, jb2,
    kt, pname)} are saved in a table {prix[0..NP-1]}.
    
    The number of terms {NT}, which is between 1 and {NP}, is inferred
    from the indices {kt}. The procedure also assembles a formula
    {termName[kt]} for each term, denoting the sum of all pairs of
    products assigned to that term, e.g. "Xom*Xop+Xmo*Xpo". It assumes
    that {belName[jb]} is the name of basis element {jb}.
    
    The values of {NP,prix,NT,termName} are returned in
    {*NP_P,*prix_P,*NT_P,*termName_P}.
    
    Comments starting with "#" and blank lines are ignored. If {verbose}
    is true, also prints the data to stderr. */

void multifok_term_write_index_table
  ( FILE *wr, 
    int32_t NP,
    multifok_term_prod_t *prix,
    bool_t verbose
  );
  /* Writes to {wr} the table {prix[0..NP-1]}, in the format expected 
    by {multifok_term_read_index_table}. */

#endif
