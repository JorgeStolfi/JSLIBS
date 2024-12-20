/* See {jspca.h}.  */
/* Last edited on 2024-12-05 11:56:31 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
/* #include <string.h> */
#include <assert.h>
#include <math.h>

#include <vec.h>
#include <rmxn.h>
#include <rn.h>
#include <affirm.h>
#include <sym_eigen.h>

#include <jspca.h>

void jspca_compute_barycenter(uint32_t nd, uint32_t nv, double D[], double w[], double d[], bool_t verbose);
  /* For each {jv} in {0..nv-1}, stores into {d[jv]} the weighted average of {D[id,jv]} 
    over all {id} in {0..nd-1} */

double *jspca_compute_covariance_matrix(uint32_t nd, uint32_t nv, double D[], double w[], double d[], bool_t verbose);
  /* Returns the covariance matrix {(D - u*d)'*W*(D - u*d)/det(W)} */

uint32_t jspca_eigen_decomp(uint32_t nv, double A[], double E[], double e[], bool_t verbose);
  /* Computes the eigendecompostion of the matrix {A}, assumed
    to be square of size {nv × nv}, symmetric, non-negative definite,
    and stored by rows.
    
    Returns the number {ne} of eigenvectors actually found. This
    number will be between 0 and {nv}, inclusive.
    
    The eigenvectors are returned as the first {ne} rows of the array {E}, which
    should be allocated with {nv} rows of {nv} columns and is linearized
    by rows. Namely, for each {ie} in {0...nv-1} and each {jv} in {0..nv-1},
    component {jv} of eigenvector number {ie} is stored in {E[ie*nv + jv]}.
    
    The output vector {e} should be allocated with {nv} elements. The procedure
    stores in {e[0..ne-1]} is the square roots of the eigenvalue
    corresponding to eigenvector {ie}.  The eigenvalues are expected 
    to be all non-negative. */

void jspca_check_ortho(uint32_t ne, uint32_t nv, double E[]);
  /* Sanity check that the rows of the {ne} by {nv} matrix {E} are orthonormal. */ 

uint32_t jspca_compute_components
  ( uint32_t nd, 
    uint32_t nv, 
    double D[], 
    double w[], 
    double d[], 
    double E[], 
    double e[],
    bool_t verbose
  )
  { 
    jspca_compute_barycenter(nd, nv, D, w, d, verbose);
    double *A = jspca_compute_covariance_matrix(nd, nv, D, w, d, verbose);
    uint32_t ne = jspca_eigen_decomp(nv, A, E, e, verbose);
    free(A);
    return ne;
  }
    
void jspca_compute_barycenter(uint32_t nd, uint32_t nv, double D[], double w[], double d[], bool_t verbose)
  { if (verbose) { fprintf(stderr, "computing the data vector barycenter {b}...\n"); }
    for (int32_t jv = 0;  jv < nv; jv++) { d[jv] = 0; }
    double sumW = 0;
    for (int32_t id = 0;  id < nd; id++)
      { double *Di = &(D[id*(int32_t)nv]);
        double wi = w[id];
        demand(wi >= 0, "invalid weight {w[id]}");
        for (int32_t jv = 0;  jv < nv; jv++) 
          { d[jv] += wi*Di[jv]; }
        sumW += wi;
      }
    demand(sumW > 0, "total weight is zero");
    for (int32_t jv = 0;  jv < nv; jv++) { d[jv] /= sumW; }
    if (verbose) { rn_gen_print(stderr, nv, d, "%+14.8f", "  d = [ ", " ", " ]\n"); }
  }
    
double *jspca_compute_covariance_matrix(uint32_t nd, uint32_t nv, double D[], double w[], double d[], bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "computing the covariance matrix {A}...\n"); }
    double *A = rmxn_alloc(nv,nv); /* The array {(D-u*d)'*W*(D-u*d)}. */
    for (int32_t jv1 = 0;  jv1 < nv; jv1++)
      { for (int32_t jv2 = 0;  jv2 <= jv1; jv2++)
          { A[jv1*(int32_t)nv + jv2] = 0; }
      }
    double sumW = 0;
    for (int32_t id = 0;  id < nd; id++)
      { double *Di = &(D[id*(int32_t)nv]);
        double vi[nv];
        rn_sub(nv, Di, d, vi);
        double wi = w[id];
        assert(wi >= 0);
        for (int32_t jv1 = 0;  jv1 < nv; jv1++) 
          { for (int32_t jv2 = 0;  jv2 <= jv1; jv2++)
              { A[jv1*(int32_t)nv + jv2] += vi[jv1]*wi*vi[jv2]; }
          }
        sumW += wi;
      }
    /* Scale the array {A} by {1/sumW} to get the covariances, and fill teh upper half: */
    assert(sumW > 0);
    for (int32_t jv1 = 0;  jv1 < nv; jv1++)
      { for (int32_t jv2 = 0;  jv2 <= jv1; jv2++)
          { double Ajj = A[jv1*(int32_t)nv + jv2] / sumW;
            A[jv1*(int32_t)nv + jv2] = Ajj;
            if (jv2 < jv1) { A[jv2*(int32_t)nv + jv1] = Ajj; }
          }
      }
    if (verbose) { rmxn_gen_print(stderr, nv, nv, A, "%+14.8f", "  A = [\n", "", "  ]\n", "    [ ", " ", " ]\n"); }
    return A;
  }

uint32_t jspca_eigen_decomp(uint32_t nv, double A[], double E[], double e[], bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "computing the eigenvectors of {A}...\n"); }
    double *ET = rmxn_alloc(nv,nv); /* Unsorted eigenvectors. */
    double *eT = rn_alloc(nv); /* Unsorted eigenvalues. */
    uint32_t ne; /* Number of eugenvalues actually computed: */
    sym_eigen(nv, A, eT, ET, &ne);
    if ((ne < nv) || verbose) { fprintf(stderr, "  eigendecomposition found %d eigenvectors out of %d\n", ne, nv); }
    if (verbose) { rmxn_gen_print2(stderr, ne, nv, ET, 1, eT,  "%+14.8f", "  E,e2 = [\n", "\n", "  ]\n", "    [ ", " ", " ]", "  "); }

    /* Copy those {ne} eigenpairs in *decreasing* magnitude order discarding small ones: */
    uint32_t ie = 0;
    for (uint32_t ke = 0; ke < ne; ke++)
      { double e2k = eT[ke];         /* Eigenvalue of eigenvector {ke}. */
        assert(e2k > -1.0e-300); /* Damn minus zero! */
        double ek = sqrt(fmax(e2k, 0.0));
        if (verbose) { fprintf(stderr, "eigenvector %3d eigenvalue = %18.10f magnitude = %18.10f", ke, e2k, ek); }
        /* copy eigenpair to row {ie} of {e,E}: */
        if (verbose) { fprintf(stderr, " renumbered %3d\n", ie); }
        e[ie] = ek;
        double *Rk = &(ET[ke*nv]); /* Row {ke} of {ET}. */
        double *Ei = &(E[ie*nv]); /* Row {ie} of {E}. */
        for (int32_t jv = 0;  jv < nv; jv++) { Ei[jv] = Rk[jv]; }
        ie++;
        /* A covariance matrix should not have negative eigenvalues: */
        if (e2k < 0) { fprintf(stderr, " ** negative eigenvalue - should not happen"); }
      }
    if (verbose) { rmxn_gen_print2(stderr, ne, nv, E, 1, e,  "%+14.8f", "  E,e = [\n", "\n", "  ]\n", "    [ ", " ", " ]", "  "); }
    
    /* !!! Move to test programs !!! */
    jspca_check_ortho(ne, nv, E);

    /* Free the temporary storage: */
    free(eT);
    free(ET);
    return ne;
  }
  
void jspca_check_ortho(uint32_t ne, uint32_t nv, double E[])
  {
    for (uint32_t ie1 = 0;  ie1 < ne; ie1++)
      { double *e1 = &(E[ie1*nv]);
        double dot11 = rn_dot(nv, e1, e1);
        assert(fabs(dot11 - 1) < 1.0e-8);
        for (uint32_t ie2 = 0;  ie2 < ie1; ie2++)
          { double *e2 = &(E[ie2*nv]);
            double dot12 = rn_dot(nv, e1, e2);
            assert(fabs(dot12) < 1.0e-8); 
          }
      }
  }
    
void jspca_decompose_data
  ( uint32_t nd,  /* Number of data points. */
    uint32_t nv,  /* Number of variables per data point. */
    double D[],  /* The data points. */
    double d[],  /* Barycenter of points. */
    uint32_t ne,  /* Number of principal components. */
    double E[],  /* Principal component vectors. */
    double C[],  /* (OUT) Eigenvector coeff matrix. */
    double P[],  /* (OUT) Projections of data points in row space of {E}. */
    double R[],  /* (OUT) Residuals of data points, orthogonal to {E}. */
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "computing the decompositon of {D} onto {E}...\n"); }
    if (verbose) { fprintf(stderr, "nd = %d  nv = %d  ne = %d\n", nd, nv, ne); }
    
    if ((C == NULL) && (P == NULL) && (R == NULL)) { /* Silly, nothing to do: */ return; }

    demand(nd >= 0, "invalid {nd}");
    demand(nv >= 0, "invalid {nv}");
    demand(ne >= 0, "invalid {ne}");
    demand(ne <= nv, "more components than variables");
    
    /* Work vectors for eacg row of {D}: */
    double vi[nv]; /* Vector from barycenter to point. */
    double ci[ne]; /* Coefficients of data vector projected on basis {E}. */
    double pi[nv]; /* Linear combination of data vectors. */
    for (uint32_t id = 0;  id < nd; id++)
      { double *Di = &(D[id*nv]); /* Data point {id}. */
        /* Subtract barycenter: */
        rn_sub(nv, Di, d, vi);
        /* Compute fitting coefficients {ci[0..nd-1]}: */
        rmxn_map_col(ne, nv, E, vi, ci);
        if (C != NULL) { double *Ci = &(C[id*ne]); rn_copy(ne, ci, Ci); }
        if ((P != NULL) || (R != NULL))
          { /* Compute projection of {vi} on row subspace of {E}: */
            rmxn_map_row(ne, nv, ci, E, pi);
            if (P != NULL) { double *Pi = &(P[id*nv]); rn_copy(nv, pi, Pi); }
            if (R != NULL) { double *Ri = &(R[id*nv]); rn_sub(nv, vi, pi, Ri); }
          }
      }
  }
   
void jspca_prv(char *name, uint32_t n, double v[], char *fmt)
  { 
    fprintf(stderr, "  %s (%d) = ", name, n);
    rn_gen_print(stderr, n, v, fmt, "[ ", " ", " ]\n");
  }
   
void jspca_prm(char *name, uint32_t m, uint32_t n, double A[], char *fmt)
  { 
    fprintf(stderr, "  %s (%dx%d) =\n", name, m, n);
    rmxn_gen_print
      ( stderr, m, n, A,
        fmt, 
        "  [\n    ", "\n    ", "\n  ]\n", 
        "[ ", " ", " ]"
      );
  }
 
void jspca_prm2(char *names, uint32_t m, uint32_t n1, double A1[], uint32_t n2, double A2[], char *fmt)
  { 
    fprintf(stderr, "  %s (%dx%d, %dx%d) =\n", names, m, n1, m, n2);
    rmxn_gen_print2
      ( stderr, m, n1, A1, n2, A2,  
        fmt, 
        "  [\n    ", "\n    ", "\n  ]\n", 
        "[ ", " ", " ]", 
        "  "
      );
  }
 
void jspca_prm3(char *names, uint32_t m, uint32_t n1, double A1[], uint32_t n2, double A2[], uint32_t n3, double A3[], char *fmt)
  { 
    fprintf(stderr, "  %s (%dx%d, %dx%d, %dx%d) =\n", names, m, n1, m, n2, m, n3);
    rmxn_gen_print3
      ( stderr, m, n1, A1, n2, A2, n3, A3,  
        fmt, 
        "  [\n    ", "\n    ", "\n  ]\n", 
        "[ ", " ", " ]", 
        "  "
      );
  }
