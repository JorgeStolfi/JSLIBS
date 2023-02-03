/* Focus detector for multi-focus stereo. */
/* Last edited on 2023-01-29 10:02:58 by stolfi */

#ifndef multifok_basis_H
#define multifok_basis_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>

#include <multifok_window.h>

/* LOCAL LINEAR OPERATORS AND BASES

  For the {multifok} modules, /local linear image operator/ (or /linop/
  for short) is a linear fnction of the {NS} samples in a monochrome
  window surrounding a pixel.
  
  The operator is defined by a list of {operator weights} {op[0..NS-1]}.
  Its value at for a window with samples {s[0..NS-1]} is
  {multifok_window_prod(NW,s,op)}; that is, {SUM{ i \in {0..NS-1] :
  s[i]*op[i]}}.
  
  A /local linear image operator basis/ (or /basis/ for short) is a list
  {bas[0..NB-1]} of linops, the /basis elements/.  A basis is represented as an
  array of pointers to the arrays of {NS} doubles, the sample weights of its linops. */
  
typedef enum {
    multifok_basis_type_CANC, /* The canonical basis, 1 at each sample and 0 at others. */
    multifok_basis_type_LAPL, /* The order 2 differential operators: "DX", "DY", "DXX", "DXY", "DYY". */
    multifok_basis_type_CUBE, /* Same as {LAPL} plus the 3rd order saddles "C3", "S3". */
    multifok_basis_type_DIFF, /* The {CUBE} operators plus the checker "Q". */
    multifok_basis_type_HART  /* The Hartley basis waves. */
  } multifok_basis_type_t;
  /* Basis types. */
  
#define multifok_basis_type_FIRST multifok_basis_type_CANC
#define multifok_basis_type_LAST multifok_basis_type_HART 

void multifok_basis_make
  ( multifok_basis_type_t bType,
    int32_t NW, 
    double ws[],
    bool_t ortho,
    int32_t *NB_P, 
    double ***bas_P, 
    char ***belName_P
  );
  /* Returns in {*bas_P} a newly allocated array {bas[0..NB-1][0..NS-1]} with a
    basis to be used to analyze the window samples. Namely {bas[i][j]}
    is the weight of basis element {i} for sample {j}. The number {NB} is chosen by
    the procedure and returned in {*NB_P}.
    
    If {bType} is {multifok_basis_type_CANC}, basis element
    {kb} will be 0 except at sample {kb}, where it will be 1.
    
    If {bType} is {multifok_basis_type_LAPL}, then {NW} must be
    3, and the basis will consist of the five differential operators "DX", "DY",
    "DXX", "DXY", "DYY".
    
    If {bType} is {multifok_basis_type_DIFF}, then {NW} must be
    3, and the basis will consist of the same operators as {LAPL},
    as well as the local image value "F", the two 3-saddles, and the 3x3 checker.
    
    If {bType} is {multifok_basis_type_HART} is false, the
    basis will consists of the Hartley waves whose frequency vectors are
    the coordinates of the samples in the window, relative to the center
    pixel.
    
    if {ws} is not {NULL}, the sample weights in each basis element {bas[i][j]}
    will be scaled by {ws[j]}.  This is useful to make the operators
    of the basis more local and more smootly varying with window position.
    
    If {ortho} is true, the basis will be modified to be orthonormal
    with respect to the the dot product {multifok_basis_prod}; that is, 
    the dot product of
    {bas[ib][0..NS-1]} and {bas[jb][0..NS-1]} is 1 if {ib=jb} and 0
    otherwise.  During this procedure, some elements may be discarded
    for being too close to linear combinations of the previous ones.
    
    If {ortho} is false, each basis element will be scaled so that 
    the dot product with any window whose samples are in {[0_1]} will
    be in the range {[-1 _ +1]} and hit at least one of the two 
    bounds. 
    
    The procedure also returns in {*belName_P} a vector of {NB} newly alocated strings 
    such that {belName[k]} is a scrutable name for basis element {bas[k]}. */
 
void multifok_basis_compute_coeffs(int32_t NW, double s[], int32_t NB, double *bas[], double coeff[]);
  /* Applies a transformation of the window samples {s[0..NS-1]} to coefficients {coeff[0..NB-1]}
    using the basis {bas[0..NB-1][0..NS-1]} and inner product weights {ws[0..NS-1]}.
    
    The basis is assumed to be orthonormal with respect to the dot product with
    sample weights {ws}, so that {coeff[kb]} is the weighted inner product of {bas[kb]}. */

char *multifok_basis_type_to_text(multifok_basis_type_t bType);
  /* Converts the basis type {bType} to text ("LAPL", "DIFF", etc.). 
    Clients should NOT call {free} on the returned strings. */
    
multifok_basis_type_t multifok_basis_type_from_text(char *name, bool_t fail);
  /* Converts the {name} ("LAPL", "DIFF", etc.) to a basis type.  
    If the {name} invalid, fails if {fail} is true, and returns {-1} if false. */

void multifok_basis_free(int32_t NB, double **bas, char **belName);
  /* Frees the storage used by {bas[0..NB-1][0..NS-1]} and {belName[0..NB-1]}. */

/* BASIS AND TERM DATA I/O */

void multifok_basis_write_elem_names(FILE *wr, int32_t NB, char *belName[]);
  /* Writes to {wr} the basis element names {belName[0..NB-1]}, one per line. */

void multifok_basis_print
  ( FILE *wr,
    int32_t NW, 
    int32_t NB, 
    double *bas[],
    char *belName[]
  );
  /* Prints to {wr} the basis {bas[0..NB-1][0..NS-1]}. */

void multifok_basis_ortho_check(FILE *wr, int32_t NW, int32_t NB, double *bas[]);
  /* Checks whether basis {bas[0..NB-1][0..NS-1]} is orthonormal. */

void multifok_basis_check(int32_t NW, multifok_basis_type_t bType, bool_t ortho);
  /* Runs various diagnostics on this module. */

#endif
