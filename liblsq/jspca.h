#ifndef jspca_H
#define jspca_H

/* JS tools for principal component analysis (PCA). */
/* Last edited on 2024-12-05 11:55:22 by stolfi */

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

/* COMPONENT EXTRACTION AND COMBINATION

  In the following procedures, all arrays that are conceptually bidimensional
  are stored as unidimensional arrays, linearized by rows.  
  
  That is, an array {A} conceptually with {m} rows and {n} columns is
  stored as an array of {m*n} {double}s such that the conceptual element
  {A[i,j]} is actually stored in {A[i*n+j]} for {i} in {0..m} and {j} in
  {0..n}.*/

uint32_t jspca_compute_components
  ( uint32_t nd, 
    uint32_t nv, 
    double D[], 
    double w[], 
    double d[], 
    double E[], 
    double e[],
    bool_t verbose
  );
  /* Computes the principal components of a set of {nd} /data points/ each
    consisting of {nv} /variables/.
    
    The input data points should be stored as the {nd} rows of a the array
    {D}, linearized by rows. That is, {D[id*nv + jv]} is conceptually
    {D[id,jv]}, the value of variable {jv} in data point {id}.
    
    The input parameter {w}, if not {NULL}, must be a vector of {nd} non-negative
    weights.  Basically, the data point {id} is implicitly replicated {w[id]}
    times in the PCA analysis.  Only the ratios of the weights matter.  If {w}
    is {NULL}, all weights are assumed to be 1.

    The output parameter {d} should be allocated as a vector of {nv}
    elements. Element {d[jv]} will be set to the average value of
    {D[id,jv]} over all rows {id} in {0..nd-1}, each weighted by
    {w[id]}.

    The output parameter {E} should be allocated as an array of {nv}
    rows and {nv} columns, but only the first {ne} rows will be used,
    for some number {ne} to be determined. 
    
    The procedure will store in these rows the {ne} principal components
    of the input data. These will be {ne} vector of {nv} components
    which are pairwise orthogonal with length 1.
    
    MATHEMATICS 
    
    The principal components, the rows of {E}, will be {ne} eigenvectors of the
    weighted covariance matrix {(D-u*d)'*W*(D-u*d)}, where {d} is
    interpreted as a 1 by {nv} row vector, {u} is the {nd} by 1 column
    vector whose elements are all 1, where {W} is the {nd} by {nd}
    diagonal matrix with the weights {w[0..nd-1]} in the diagonal. 
    
    The output parameter {e} should be allocated with {nv} elements, but
    only the first {ne} will be used. The procedure will set{e[ie]} to
    the square root of the eigenvalue associated with eigenvector {ie}.
    
    INTUITION
    
    If the data set {D} is viewed as a cloud of {nd} point masses in
    {\R^nv}, then the returned vector {d} will be the barycenter of the
    cloud. The {ne} rows of {E} will be the {ne} orthogonal directions
    in {\R^nd} along which the mass will be spread out to the largest
    extent, relative to the barycenter {d}.   Then {e[ie]} will be the
    weighted standard deviation of the points, minus {d}, projected
    along direction {ie}. */

void jspca_decompose_data
  ( uint32_t nd,  /* Number of data points. */
    uint32_t nv,  /* Number of variables per data point. */
    double D[],   /* The data points. */
    double d[],   /* Barycenter of points. */
    uint32_t ne,  /* Number of principal components. */
    double E[],   /* Principal component vectors. */
    double C[],   /* (OUT) Eigenvector coeff matrix. */
    double P[],   /* (OUT) Projections of data points in row space of {E}. */
    double R[],   /* (OUT) Residuals of data points, orthogonal to {E}. */
    bool_t verbose
  );
  /* Computes the coefficients of a set of principal components
    in the  original data points.
  
    The array {D}, conceptually with {nd} rows by {nv} columns, should contain a set
    of {nd} data points with {nv} variables each.
    
    The input parameter {d} should be an array of {nv} elements.
    
    The input parameter {E} must be a matrix of {ne} rows and {nv}
    columns, whose rows are pairwise orthogonal with length 1.
    
    The output parameter {C}, if not {NULL}, must be an array of {nd}
    rows and {ne} columns. The parameters {U} and {R}, if not {NULL},
    must each be an array with {nd} rows and {nv} columns.
    
    For each data point index {id} in {0..nd-1}, let {Di[0..nv-1]},
    {Ci[0..ne-1]}, {Pi[0..nv-1]} and {Ri[0..nv-1]} be row {id} of
    {D,C,P,R}, respectively. That is, {Ci[ie] = C[id*ne + ie]}, {Di[jv]
    = D[id*nv + jv]}, and the same for {P} and {R}.
    
    The procedure decomposes each relative data vector {Di - d} into a linear
    combination {Pi[0..nv-1]} of the first {ne} rows of {E} with
    coefficients {Ci[0..ne-1}, and a residue {Ri[0..nv-1]} that is
    perpendicular to all those rows.  
    
    That is, on return we will have {C = (D - u*d)*E'}, {P = C*E}, {D =
    P + R + u*d}, {R*E' = z}; where {z} is the {nd} by {ne} zero matrix.
    
    If either of {C}, {P} or {R} is {NULL}, the corresponding array is
    simply not returned. Either {P} or {R} may be the same as {D}.
    
    The arrays {D}, {d} and {E} need not be the result of
    {jspca_compute_components}, although they usually are. */

/* DEBUGGING */

void jspca_prv(char *name, uint32_t n, double v[], char *fmt);
  /* Prints to [stderr} the row vector {v[0..n-1]}. */
  
void jspca_prm(char *name, uint32_t m, uint32_t n, double A[], char *fmt);
  /* Prints to [stderr} the array {A}, assumed {m} by {n}.  Each element is printed with the format {fmt}. */

void jspca_prm2(char *names, uint32_t m, uint32_t n1, double A1[], uint32_t n2, double A2[], char *fmt);
  /* Prints to [stderr} the arrays {A1}, assumed {m} by {n1}, and {A2},
    assumed {m} by {n2}, side by side. Each element is printed with the
    format {fmt}. */

void jspca_prm3(char *names, uint32_t m, uint32_t n1, double A1[], uint32_t n2, double A2[], uint32_t n3, double A3[], char *fmt);
  /* Prints to [stderr} the arrays {A1}, assumed {m} by {n1}, {A2},
    assumed {m} by {n2}, and {A3}, assumed {m} by {n3}, side by side.
    Each element is printed with the format {fmt}. */
          
#endif
