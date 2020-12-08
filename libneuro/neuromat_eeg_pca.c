/* See {neuromat_eeg_pca.h}.  */
/* Last edited on 2017-10-18 06:13:13 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
/* #include <string.h> */
#include <assert.h>
#include <math.h>

#include <vec.h>
#include <rmxn.h>
#include <rn.h>
#include <affirm.h>
#include <sym_eigen.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_pca.h>

int neuromat_eeg_pca_eigen_decomp(int ne, double *A, double minMag, double *Ev, double emag[])
  {
    demand(minMag >= 0, "bad {minMag}");
    
    /* Convert {A} to symmetric tridiagnonal form {T}: */
    double *dT = rn_alloc(ne); /* Diagonal of {T}, then eigenvalues. */
    double *eT = rn_alloc(ne); /* Subdiagonal of {T}. */
    double *R = rmxn_alloc(ne,ne); /* Rotation matrix, then eigenvectors. */
    syei_tridiagonalize(ne, A, dT, eT, R);

    /* Compute the eigensystem for {T}: */
    int mv; /* Number of eugenvalues actually computed: */
    int absrt = 0; /* Sort eigenvalues by signed value. */
    syei_trid_eigen(ne, dT, eT, R, &mv, absrt);
    if (mv < ne) 
      { fprintf(stderr, "eigendecomposition found only %d eigenvectors of %d\n", mv, ne); }

    /* Copy those {mv} eigenpairs in *decreasing* magnitude order discarding small ones: */
    int nv = 0;
    int kv;
    for (kv = mv-1; kv >= 0; kv--)
      { double evk = dT[kv];
        double emagk = sqrt(fabs(evk)); 
        fprintf(stderr, "eigenvector %3d eigenvalue = %18.10f magnitude = %18.10f", kv, evk, emagk);
        /* Note that sign of {evk} is important here: */
        if (evk >= minMag*minMag)
          { /* Eigenvalue is non-negative and large enough: copy eigenpair to {emag,Ev}: */
            int iv = nv; /* Next free position. */
            fprintf(stderr, " renumbered %3d", iv);
            emag[iv] = emagk;
            double *Rk = &(R[kv*ne]);   /* Row {kv} of {R}. */
            double *Evi = &(Ev[iv*ne]); /* Row {iv} of {Ev}. */
            int je;
            for (je = 0; je < ne; je++) { Evi[je] = Rk[je]; }
            nv++;
          }
        else
          { fprintf(stderr, " discarded"); }
        if (evk < 0) { fprintf(stderr, " ** negative eigenvalue"); }
        fprintf(stderr, "\n");
      }

    /* Free the temporary storage: */
    free(eT);
    free(dT);
    free(R);
    return nv;
  }

void neuromat_eeg_pca_compute_fitting_matrix(int np, int ne, double *P, double *Q)
  {
    double *R = rmxn_alloc(np,np); /* The array {P*P'}. */
    rmxn_mul_tr(np, np, ne, P, P, R);
    rmxn_inv(np,R,Q);
    /* Check the inverse: */
    double *S = rmxn_alloc(np,np); /* Should be the identity. */
    rmxn_mul(np,np,np,R,Q,S);
    int i, j;
    for (i = 0; i < np; i++)
      { for (j = 0; j < np; j++)
          { double Sij = S[i*np + j];
            double Iij = (i == j ? 1.0 : 0.0);
            if (fabs(Sij - Iij) > 1.0e-6)
              { fprintf(stderr, "  S[%d,%d] = %24.16e\n", i, j, Sij);
                demand(FALSE, "inversion failed, aborted");
              }
          }
      }
    free(R);
    free(S);
  }
    
void neuromat_eeg_pca_fit_patterns
  ( int nt, 
    int ne, 
    double **val, 
    int np, 
    double *P, 
    double *Q, 
    double **coeff,
    double **vpara,
    double **vperp
  )
  {
    if ((coeff == NULL) && (vpara == NULL) && (vperp == NULL))
      { /* Silly, nothing to do: */
        return;
      }

    demand(ne >= 0, "invalid {ne}");
    demand(np >= 0, "invalid {np}");
    demand(np <= ne, "more components than electrodes");
    
    /* Per-frame work vectors: */
    double bi[np]; /* Independent vector of linear system. */
    double ci[np]; /* Coefficients of patterns. */
    double wi[ne]; /* Linear combination of patterns. */
    for (int it = 0; it < nt; it++)
      { double *vi = val[it]; /* Original electrode values {vi[0..ne-1]}. */
        /* Compute fitting coefficients {ci[0..np-1]}: */
        rmxn_map_col(np, ne, P, vi, bi);
        rmxn_map_col(np, np, Q, bi, ci);
        if (coeff != NULL) { rn_copy(np, ci, coeff[it]); }
        if ((vpara != NULL) || (vperp != NULL))
          { /* Compute linear combination {wi]0..ne-1]} of patterns: */
            rmxn_tr_mul(np, ne, 1, P, ci, wi);
            if (vpara != NULL) { rn_copy(ne, wi, vpara[it]); }
            if (vperp != NULL) { rn_sub(ne, vi, wi, vperp[it]); }
          }
      }
  }

