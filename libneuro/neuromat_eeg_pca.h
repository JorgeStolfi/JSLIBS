#ifndef neuromat_eeg_pca_H
#define neuromat_eeg_pca_H

/* Principal componnet analysis toos for NeuroMat. */
/* Last edited on 2021-08-29 11:44:46 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

/* COMPONENT EXTRACTION AND COMBINATION

  In the following procedures, {P} is a list of {np} potential patterns
  for {ne} electrodes.
  
  More precisely, {P} is a linearized array with {np} rows and {ne}
  columns, such that {P[ip*ne + ie]} is the value of electrode {ie} in
  pattern {ip}.
  
  The array {vin} must contain potentials for {nt} frames and {ne}
  electrodes. More precisely, the procedure assumes that {vin[it][ie]},
  is the potential of electrode {ie} in frame {it}, for {it} in
  {0..nt-1} and {ie} in {0..ne-1}.
  
  The array {Q}, when given, must be  a linearized matrix {Q} with {np} rows
  and {np} columns, which is the inverse of {P*P'} where {P'} is the
  transpose of {P}. */

int32_t neuromat_eeg_pca_eigen_decomp(int32_t ne, double *A, double minMag, double *Ev, double emag[]);
  /* Computes the eigendecompostion of the matrix {A}, assumed
    to be square of size {ne × ne}, symmetric, and stored by rows.

    Returns the number {nv} of eigenvectors actually found, discarding
    those whose eigenvalues are less than {minMag^2}. This number
    will be between 0 and {ne}, inclusive.
    
    For {iv} in {0...nv-1}, the eigenvector number {iv} is stored in
    {Ev[iv*ne + je]} for {je} in {0..ne-1}. The square root of the
    corresponding eigenvalue is stored in {emag[iv]}. */

void neuromat_eeg_pca_compute_fitting_matrix(int32_t np, int32_t ne, double *P, double *Q);
  /* Stores into {Q} the inverse of {P*P'}.  Also checks whether the inverse is 
    valid. */

void neuromat_eeg_pca_fit_patterns
  ( int32_t nt, 
    int32_t ne, 
    double **val, 
    int32_t np, 
    double *P, 
    double *Q, 
    double **coeff,
    double **vpara,
    double **vperp
  );
  /* Computes {np} coefficients of {np} given electrode patterns, for 
    each frame in {val}.  
    
    For each input frame {val[it]}, the program considers its first {ne}
    samples to be a row vector {vi[0..ne-1]}, and decomposes {vi} into a
    linear combination {wi} of rows of {P}, and a residue {ui} that is
    perpendicular to all those rows. 
    
    If {coeff} is not {NULL}, stores the {np}
    coefficients of that linear combination into {ci[0..np-1] =
    coeff[it][0..np-1]}.
    
    If {vpara} is not {NULL}, stores the pattern combination {wi} itself into
    {vpara[it][0..ne-1]}. 

    If {vperp} is not {NULL}, stores the residual {ui} orthogonal to all
    patterns in {vperp[it][0..ne-1]}.
    
    Either {vpara} or {vperp} may be the same as {val}.
    
    The matrix {Q} should be an inverse of {P*P'} as computed by
    {neuromat_eeg_pca_compute_fitting_matrix}. The coefficients
    {ci[0..np-1]} are computed by the product {Q*P*vi'}. */

void neuromat_eeg_pca_combine_patterns
  ( int nt, 
    int ne, 
    int np, 
    double *P, 
    double **coef, 
    int mp, 
    int ip[], 
    double **vout
  );
  /* Creates an EEG datset from a selected subset of patterns from {P}
    and their coefficients from {coef}.

    Assumes that {ip[0..mp-1]} are integers in the range {0..np-1}.
    
    The procedure stores into channels {0..ne-1} of each frame
    {vout[it]} the EEG electrode potentials corresponding to the
    linear combination of rows {ip[0..mp-1]} of {P}, with their
    respective coefficients from {coef[it]}.  That is
    
      {vout[it][ie] = SUM{ coef[it][ip[k]]*P[ip[k]*ne + ie] : k \in 0..mp-1}}
      
    for {it} in {0..nt-1} and {ie} in {0..ne-1}.  Any excess channels in {vout}
    are left unchanged. */
          
#endif
