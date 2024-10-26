/* Quadratic operator terms for multi-focus stereo. */
/* Last edited on 2024-10-15 17:03:56 by stolfi */

#ifndef multifok_term_H
#define multifok_term_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <multifok_window.h>
#include <multifok_basis.h>

/* QUADRATIC TERMS FOR SHARPNESS ESTIMATION

  The {multifok} library uses local focus sharpness estimators
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
    char *pname; /* Name (formula) of product, e.g. "FX*FX" of "Smm*Spp". */
  } multifok_term_prod_t;
  /* Record describing a product of two coeffs {coeff[jb1]*coeff[jb2]}
    in some local linear basis, that contributes to one indpendent term
    of a quadratic local focus sharpness estimator.
    
    The product must be in canonical form, meaning {jb1<=jb2}. The name
    shoud be "{el1}*{el2}" where {el1} and {el2} are the names of the
    correspoing basis elements.
    
    As a special case, {jb1} and {jb2} may be {-1}, representing the
    empty product whose value is constant 1. The {pname} then should be
    "1". */
typedef struct multifok_term_set_t 
  { /* Basis of linear window operators: */
    multifok_basis_t *basis;     /* Basis of linear window ops. */
    int32_t NT;                  /* Number of quadratic terms to analyze. */
    char **termName;             /* Term names. */
    int32_t NP;                  /* Number of basis coeff products in the quadratic terms. */
    multifok_term_prod_t *prix;  /* Mapping of basis element index pairs to terms. */
  } multifok_term_set_t;
  /* Tables describing a set of quadratic window operator terms that may be
    combined into a focus estimator.
    
    The {NT} quadratic terms are described by the table {prix[0..NP-1]}.
    Specifically, {prix[it]} specifies the indices {jb1,jb2} of two 
    basis elements, and the index {kt} in {0..NT-1} of the quadratic term 
    to which that product is to be added. */

void multifok_term_values_from_basis_coeffs
  ( double coeff[],
    multifok_term_set_t *tset,
    double term[]
  );
  /* Assumes that {coeff[0..NB-1]} are the coefficients of the window
    samples in the linear operator basis {tset.basis}, where {NB =
    tset.basis.NB}. Computes the quadratic operator terms {term[0..NT-1]}
    specified by {tset}, where {NT=tset.NT}.
    
    Namely, initializes {term[0..NT-1]} with zeros; then, for each entry
    {tset.prix[kp] = (jb1,jb2,kt,pname)} with {kp} in {0..tset.NP-1},
    computes the product {coeff[jb1]*coeff[jb2]} and adds it to
    {term[kt]}. As a sepcial case, if {jb1=jb2=-1}, adds 1.0 to
    {term[kt]}. */ 

void multifok_term_set_names_write(FILE *wr, int32_t NT, char *termName[]);
  /* Writes to {wr} the the term names (formulas) {termName[0..NT-1]}, one per line. */
 
void multifok_term_set_names_write_named
  ( char *outPrefix, 
    int32_t NT, 
    char *termName[]
  ); 
  /* Same as {multifok_term_set_write_names}, but 
    to a file file "{outPrefix}-tnames.txt" instead of a {FILE} descriptor. */

multifok_term_set_t *multifok_term_set_read
  ( FILE *rd,
    bool_t weights,
    multifok_basis_t *basis,
    double *wt_P[],
    bool_t verbose
  );
  /* Reads a set of terms {tset} from {rd}.
  
    The number of terms {NT=tset.nt} and the number of products 
    {NP=tset.NP} are inferred from the file's contents.
    
    The file {rd} must contain {NT} data lines, one per term.  Each line
    must have the term index {kt} in {0..NT-1}, followed optionally by a numeric
    weight {wt[kt]} (if {weights} is true) and the term's name {tset.termName[kt]}.
    If {weights} is false, {wt[kt]=1.0} is assumed. 
    
    The term name must be a formula consisting either of the string "1"
    or a sum of products of basis element names, like "Smm*Spp+Spm+Spp".
    These strings are parsed to yield the produc index table
    {tset.prix[0..NP-1]}.
    
    Returns the terms as a newly allcoated {multifok_term_set_t}
    record, and the weight vector {wt} in {*wt_P}. However, if {wt_P} is
    {NULL}, the weights are discarded.
    
    Comments starting with "#" and blank lines are ignored. If {verbose}
    is true, also prints the data to stderr.
    
    */
  
multifok_term_set_t *multifok_term_set_read_named
  ( char *fname,
    bool_t weights,
    multifok_basis_t *basis,
    double **wt_P,
    bool_t verbose
  );
  /* Same as {multifok_term_set_read} but reads
    from file {fname} instead of a {FILE descriptor. */

void multifok_term_set_write(FILE *wr, multifok_term_set_t *tset, bool_t weights, double wt[]);
  /* Writes to {wr} the weights {wt[0..NT-1]} and the 
    term names (formulas) {tset.termName[0..NT-1]}, where {NT=tset.NT}, one per line,
    in the format expected by {multifok_term_set_read(rd,weights,...)}. 
    
    If {weights} is true but {wt} is {NULL}, assumes it is a vector of ones. */

void multifok_term_set_write_named(char *outPrefix, multifok_term_set_t *tset, bool_t weights, double wt[]); 
  /* Same as {multifok_term_set_write}, but 
    to a file file "{outPrefix}-twts.txt" instead of a {FILE} descriptor. */

void multifok_term_set_product_table_write
  ( FILE *wr, 
    int32_t NP,
    multifok_term_prod_t *prix,
    bool_t verbose
  );
  /* Writes to {wr} the table {prix[0..NP-1]}, in the format expected 
    by {multifok_term_read_index_table}. */

void multifok_term_set_product_table_write_named
  ( char *outPrefix, 
    int32_t NP, 
    multifok_term_prod_t prix[],
    bool_t verbose
  );
  /* Same as {multifok_term_write_index_table}, but writes to a file "{outPrefix}
    "{outPrefix}-prix.txt" instead of a {FILE} descriptor. */

#endif
