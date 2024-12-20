/* Focus detector for multi-focus stereo. */
/* Last edited on 2024-12-05 14:14:08 by stolfi */

#ifndef multifok_basis_H
#define multifok_basis_H

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
    multifok_basis_type_LAPL, /* The differential operators up tp order 2: "F", "FX", "FY", "FXX", "FXY", "FYY". */
    multifok_basis_type_CUBE, /* Same as {LAPL} plus the 3rd order saddles "C3", "S3". */
    multifok_basis_type_DIFF, /* The {CUBE} operators plus the checker "Q". */
    multifok_basis_type_HART  /* The Hartley basis waves. */
  } multifok_basis_type_t;
  /* Basis types. */
  
#define multifok_basis_type_FIRST multifok_basis_type_CANC
#define multifok_basis_type_LAST multifok_basis_type_HART 

typedef struct multifok_basis_t 
  { /* Basis of linear window operators: */
    multifok_basis_type_t type;  /* Basis type. */
    uint32_t NW;                 /* Width and height of window. Must be odd. */
    double *ws;                  /* The window sample weights are {ws[0..NS-1]} where {NS=NW*NW}. */
    bool_t ortho;                /* True iff basis has been orthonormalized. */
    uint32_t NB;                 /* Number of basis elements. */
    double **bas;                /* The pixel values of the basis elements are {bas[0..NB-1][0..NS-1]}. */
    char **belName;              /* The names of basis elements are {belName[0..NB-1]}. */
  } multifok_basis_t ;
  /* Tables describing a basis of linear window operators. The basis
    elements are {bas[0..NB-1][0..NS-1]}, where {NS=NW*NW}, and their
    names are {belName[0..NB-1]}. Specifically, {bas[ip][ks]} is the
    value of window sample {ks} in element {ib}.
    
    If {ortho} is true, the basis elements specified by the {type} have
    been made orthogonal under the weighted dot product with window
    sample weights {ws[0..NS-1]}. */

multifok_basis_t *multifok_basis_make
  ( multifok_basis_type_t type,
    uint32_t NW, 
    double ws[],
    bool_t ortho,
    bool_t verbose
  );
  /* Returns a {basis} record that describes a basis of linear window 
    oerators for a square window of size {NW} (which must be odd)
    with window mask weights {ws[0..NS-1]} wgere {NS=NW*NW}. 
    
    The record {basis} and the internal tables are all newly allocated
    by the procedure.
    
    The {type} parameter specifies the size and nature of the initial basis:
    
      * {multifok_basis_type_CANC}: {NW} must be 3. Each basis element will be
      0 except at one sample position {ks}, where it will be 1. 
      The elements will be sorted by increasing distance 
      of the sample from the window center.    
      The element names  will be "S{X}{Y}" where {X} and {Y} are the 
      indices of each sample relative to the center sample, encoded as in
      {multifok_window_sample_names}.

      * {multifok_basis_type_LAPL}: {NW} must be 3. The basis will
      consist of the six differential operators up to order 2: "F",
      "FX", "FY", "FXX", "FXY", "FYY".

      * {multifok_basis_type_DIFF}: {NW} must be 3. The basis will
      consist of the same operators as {LAPL}, as well as the two 
      3-saddles, and the 3x3 checker.

      * {multifok_basis_type_HART}: The basis will consists of the {NS}
      Hartley waves. Element {ks} will be the wave with frequency vector
      {(fx,fy)},] where {fx} and {fy} are the indices of sample {ks} of
      the window relative to the center sample, in {-HW..+HW}. The
      element name will be "H{FX}{FY}" where {FX} and {FY} are the
      frequencies {fx} and {fy}, encoded as in
      {multifok_window_sample_names}.

    If {ortho} is true, the basis sepcified by {type}, as above, will be
    modified to be orthonormal with respect to the the dot product
    {multifok_basis_prod} with weights {ws[0..NS-1]}; that is, the dot
    product of {bas[ib][0..NS-1]} and {bas[jb][0..NS-1]} is 1 if {ib=jb}
    and 0 otherwise. During this procedure, some elements may be
    discarded for being too close to linear combinations of the previous
    ones.
    
    If {ortho} is false, each basis element will be scaled so that 
    the dot product with any window whose samples are in {[0_1]} will
    be in the range {[-1 _ +1]} and hit at least one of the two 
    bounds. 
    
    The procedure also returns in {*belName_P} a vector of {NB} newly alocated strings 
    such that {belName[k]} is a scrutable name for basis element {bas[k]}. */
 
void multifok_basis_compute_coeffs(double s[], multifok_basis_t *basis, double coeff[]);
  /* Applies a transformation of the window samples {s[0..NS-1]}, where
    {NS=NW*NW} and {NW=basis.NW}, to coefficients {coeff[0..NB-1]} using
    the basis {basis->bas[0..NB-1][0..NS-1]} and inner product weights
    {basis->ws[0..NS-1]}.
    
    The basis must be orthonormal with respect to the dot product with
    those sample weighst, so that {coeff[kb]} is just the weighted inner
    product of {basis->bas[kb]} and {s}. */

char *multifok_basis_type_to_text(multifok_basis_type_t type);
  /* Converts the basis type {type} to text ("LAPL", "DIFF", etc.). 
    Clients should NOT call {free} on the returned strings. */
    
multifok_basis_type_t multifok_basis_type_from_text(char *name, bool_t fail);
  /* Converts the {name} ("LAPL", "DIFF", etc.) to a basis type.  
    If the {name} invalid, fails if {fail} is true, and returns {-1} if false. */

void multifok_basis_free(multifok_basis_t *basis);
  /* Frees the storage used by {basis}, 
    including all internal tables and strings (except the weight table {basis->ws),
    and the record {basis} itself. */

/* BASIS DATA I/O */

void multifok_basis_elem_names_write(FILE *wr, uint32_t NB, char *belName[]);
  /* Writes to {wr} the basis element names {belName[0..NB-1]}, one per line. */

void multifok_basis_elem_names_write_named(char *outPrefix, multifok_basis_t *basis);
  /* Writes the basis element names {basis.belName[0..basis.NB-1]} to a file file
    "{outPrefix}-bnames.txt" one per line. */

void multifok_basis_print(FILE *wr, multifok_basis_t *basis);
  /* Prints to {wr} the given {basis}. */

void multifok_basis_ortho_check(FILE *wr, multifok_basis_t *basis);
  /* Checks whether the given {basis} is orthonormal.  Prints 
    to {wr} the moment matrix (pairwise element dot products). */

void multifok_basis_module_check(uint32_t NW, multifok_basis_type_t type, bool_t ortho, bool_t verbose);
  /* Runs various diagnostics on this module. */

#endif
