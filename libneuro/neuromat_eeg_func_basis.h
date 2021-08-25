#ifndef neuromat_eeg_func_basis_H
#define neuromat_eeg_func_basis_H

/* NeuroMat tools to interpolate values over flat scalp image. */
/* Last edited on 2021-08-22 09:03:28 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <r2.h>
#include <r3.h>
#include <float_image.h>

#include <neuromat_eeg.h>
    
typedef void neuromat_eeg_func_basis_eval_t (int32_t ne, double bval[], r3_t *p); 
  /* Type of a funtion that evaluates a function basis with {ne} elements at the point {p},
    ans returns the result in {bval[0..ne-1]}. */

void neuromat_eeg_func_basis_eval
  ( int32_t ne, 
    double bval[],
    r3_t *p3D,
    neuromat_eeg_func_basis_eval_t mother,
    double L[],
    bool_t norm
  );
  /* Computes an approximation/interpolation basis {bas[0..ne-1]} for {ne}
    electrodes, evaluates them at the point {p3D}, and stores the
    results in {bval[0..ne-1]}.
    
    See {neuromat_eeg_geom.h} for the definition of the /schematic/ and
    /idealized/ coordinate systems. The idealized 3D scalp is assumed to
    have radius {srad.c[k]} along axis {k}.
    
    Each basis function {bas[k]} is a real function of idealized 3D scalp coordinates.
    and it is intended to give the contribution of the potential measured
    at electrode {k} on the interpolated value at eavery other point of the scalp.
    ]at each point of the  
    
    The basis functions are derived from the /mother functions/ {mof[0..ne-1]}.
    The function {mother(ne,bval,p)} should set {bval[k]} to 
    where {mof[k](p)} for all {k}.
    
    Then, if {L} is not {NULL}, {bval} (viewed as a column vector) is replaced by {L * bval}.
    The matrix {L} can be used, for example, to create an interpolating function basis 
    out of an arbitrary one.
    
    Then, if {norm} is true, the values {bval[0..ne-1]} are scaled so that their sum is 1.
    This will ensure that the the interpolation basis will reproduce constant function: that is,
    if all electrode potentials are the same value {V}, the interpolated potentials will be
    {v} everywhere. */
    
/* UTILITIES */

double *neuromat_eeg_func_basis_nearest_dists(int32_t ne, r3_t pos3D[]);
  /* Given {ne} idealized 3D positions {pos3D[0..ne-1]}, returns a vector
    {erad[0..ne-1]} such that {erad[ie]} is the Euclidean distance from
    electrode {ie} to its nearest neighbor. If there is only one
    electrode, sets {erad[0]} to {+INF}. */
    
double *neuromat_eeg_func_basis_lagrangian_matrix(int32_t ne, neuromat_eeg_func_basis_eval_t mother, r3_t pos3D[]);
  /* Returns an {ne} by {ne} matrix {L} that converts an arbitrary function basis  {mof[0..ne-1]} into 
    a function basis {bas[0..ne-1]} that is interpolating at specified electrode position.
    
    The electrode positions on the idealized 3D scalp are assumed to be
    {pos3D[0..ne-1]}.
    
    The input basis {mof} is defined by the {mother} procedure, namely
    {mof[k](p) = mother(k,p)}. These functions must be linearly
    independent when restricted to the points {pos3D[0..ne-1]}; that is,
    the array {M} such that {M[i,j] = mof[i](pos3D[j])} should be
    non-singular.
    
    The new basis {bas} is defined as {bas[i](p) = \sum{j=0}^{ne-1}
    L[i,j]*mof[j](p)}. It will be such that {bas[k](pos3D[i]])} will be
    1 if {k=i} and 0 if {k!=i}.  The matrix {L} will be linearized
    by rows; that is, {L[i,j]} will be actually {L[i*ne + j]}.. */
    
/* SOME MOTHER FUNCTIONS 
    
    The folowing functions may be useful when defining mother functions.
    
    Some mother functions below require computing distances between points in
    the scalp. Since the actual positions of the electrodes on the
    scalp are unknown, the distances are computed from their
    idealized 3D positions. */

double neuromat_eeg_func_basis_shepard_weight(r3_t *p, r3_t *ctr, double rho, double sigma, double ord);
  /* A modified Shepard radial element, basically {f(p) = 1/z^ord} where {z = dist(p,ctr)/rho}. 
    The support of the function is trimmed by multiplying it by a Gaussian bell with
    radius {sigma} and by fudging it so that it is still finite (but very large) when 
    {p} is {ctr} or very close to it.  Thus after normalization, the basis will be interpolating.
    The elements are positive everywhere. */
    
double neuromat_eeg_func_basis_gauss_bell(r3_t *p, r3_t *ctr, double sigma);
  /* A Gaussian bell {f(p) = exp(-z^2/2)} where {z = dist(p,ctr)/sigma}. */

double neuromat_eeg_func_basis_mexican_hat(r3_t *p, r3_t *ctr, double sigma, double tau);
  /* A Gaussian-modified radial {sinc} function, namely {f(p) = exp(-z^2/2)*sin(pi*t/2)}
    where {z = dist(p,ctr)/tau} and {t = dist(p,ctr)/sigma}. */

double neuromat_eeg_func_basis_voronoi_ind(r3_t *p, int32_t ie, int32_t ne, r3_t pos3D[]);
  /* The Voronoi region indicator, namely 1 if the nearest electrode to {p} is {pos3D[ie]},
    0 otherwise. */

#endif
