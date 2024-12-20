/* See {lsq_robust_mf.h} */
/* Last edited on 2024-12-05 10:33:16 by stolfi */

#define lsq_robust_mf_C_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <rmxn.h>
#include <jsmath.h>
#include <gausol_solve.h>
#include <rn.h>

#include <lsq.h>
#include <lsq_array.h>
#include <lsq_robust_mf.h>

void lsq_robust_mf_fit
  ( int32_t nx,       /* Number of independent variables. */
    int32_t nf,       /* Number of dependent variables (functions to fit). */
    int32_t nt,       /* Number of cases to generate. */
    double X[],   /* Sample values of the independent variables ({nt} by {nx}). */
    double F[],   /* Corresponding measured values of the dependent variables ({nt} by {nf}. */
    double W[],   /* Corresponding weights ({nt} elements). */
    int32_t maxiter,  /* Max iteration count. */
    double U[],   /* (OUT) Fitted linear transformation matrix. */
    double K[],   /* (OUT) Assumed principal axes of normal errors, or NULL. */
    double L[],   /* (OUT) Assumed principal variances of normal errors, or NULL. */
    double P[],   /* (OUT) Outlier probabilities, or NULL. */
    lsq_robust_mf_report_t *report,
    bool_t verbose
  )
  {
    /* bool_t debug = (report != NULL); */
    bool_t debug = TRUE;
    bool_t verbacc = FALSE;
    
    auto void gen_case(int32_t it, int32_t nx, double xi[], int32_t nf, double fi[], double *wiP);
      /* Returns row {it} of {X} in {xi}, row {it} of {F} in {fi}, and {W[it]} in {*wiP}. */
      
    double *A = rmxn_alloc(nx,nx); /* Moment matrix. */
    double *B = rmxn_alloc(nx,nf); /* Right-hand side. */

    /* Compute the system matrices {A,B}: */
    lsq_compute_matrix_and_rhs(nx, nf, nt, gen_case, A, B, verbacc);
    
    /* The matrix does not change: */
    rmxn_inv(nx, A, A);

    /* Solve first without correction: */
    rmxn_mul(nx, nx, nf, A, B, U);
    if (report != NULL) { report(0, U, F, NULL); }

    if (maxiter > 0)
      { 
        /* Apply the outlier exclusion method: */
        double *Fa = rmxn_alloc(nt,nf); /* Current approximation. */
        double *Fc = rmxn_alloc(nt,nf); /* Deviations first, fake values later. */

        /* Estimated statistics: */
        double *Pc = (P != NULL ? P : rn_alloc(nt)); /* Probability of each datum to be inlier. */
        double P_gud, P_bad;         /* A priori probability of being inlier or outlier. */
        double *E_gud = rn_alloc(nf);  /* Mean error of inliers. */
        double *E_bad = rn_alloc(nf);  /* Mean error of outliers. */
        double *K_gud = (V != NULL ? V : rmxn_alloc(nf,nf));  /* Covariance matrix of inliers. */
        double *K_bad = rmxn_alloc(nf,nf);  /* Covariance matrix of outliers. */
        
        int32_t iter;
        rn_all(nt, 0.5, Pc); /* A priori, inliers and outliers are equally likely: */
        for (iter = 1; iter <= maxiter; iter++)
          { 
            if (debug) { fprintf(stderr, "    iteration %d\n", iter); }

            /* Set {Fa[0..n-1]} to current approximation: */
            rmxn_mul(nt, nx, nf, X, U, Fa);
            
            /* Set {Fc} to the deviations from approimation: */
            rmxn_mix(nt, nx, 1.0, F, -1.0, Fa, Fc);
            
            /* Recompute parameters of inlier deviations: */
            lsq_robust_mf_compute_stats(nt, nf, Fc, W, Pc, TRUE,  NULL,  E_gud, K_gud, L_gud, &P_gud, debug);
            
            /* Recompute parameters of outlier values: */
            lsq_robust_mf_compute_stats(nt, nf, F,  W, Pc, FALSE, L_gud, E_bad, K_bad, L_bad, &P_bad, debug);

            assert(fabs(P_bad + P_gud - 1.0) < 0.0001);
           
            /* Recompute data point inlier/outlier probabilities and adjusted data {Fc[0..n-1]}: */
            int32_t k;
            double E_gud_k[nf]; /* Assumed average of inlier error distribution for each dependent variable. */
            for (k = 0; k < nt; k++) 
              { /* Decide the probability {Pc[k]} of each data record {X[k,*},F[k,*]} being an inlier: */
                /* Grab the sample value vector: */
                double *Fk = &(F[k*nf]);
                /* Compute its assumed average if inlier: */
                double *Fak = &(Fa[k*nf]);
                rn_sub(nf, Fak, E_gud, E_gud_k);
                /* Compute the probability of {Fk} being inlier: */
                double pk = lsq_robust_mf_bayes(nf, Fk, P_gud, E_gud_k, K_gud, L_gud, E_bad, K_bad, L_bad, verbose);
                Pc[k] = pk;
                /* Set {Fc[k]} to the expected value of the inlier part of {F[k]}: */
                double *Fck = &(Fc[k*nf]);
                rn_mix(nf, pk, Fk, 1-pk, Fak, Fck);
              }
              
            /* Recompute coeffs using {Fc} instead of {F}: */
            lsq_array_compute_rhs(nx, nf, nt, X, Fc, W, B);
            lsq_solve_system(nx, nf, A, B, 0,NULL,NULL, U,NULL, verbose); 
            rmxn_mul(nx, nx, nf, A, U, Fa);

            if (report != NULL) { report(iter, U, Fa, Pc); }
          }
        free(Fa);
        free(Fc);
        if (Pc != P) { free(Pc); }
        free(E_gud);
        free(E_bad);
        if (V_gud != V) { free(V_gud); }
        free(V_bad);
      }
    
    void gen_case(int32_t it, int32_t nx, double xi[], int32_t nf, double fi[], double *wiP)
      { double *Xi = &(X[it*nx]);
        double *Fi = &(F[it*nf]);
        int32_t j;
        for (j = 0; j < nx; j++) { xi[j] = Xi[j]; }
        for (j = 0; j < nf; j++) { fi[j] = Fi[j]; }
        (*wiP) = W[it];
      }
  }

void lsq_robust_mf_compute_stats
  ( int32_t nt,
    int32_t nf,
    double Y[], 
    double W[], 
    double P[], 
    double Q[],
    bool_t good, 
    double E[], 
    double V[],
    double K[],
    double L[],
    double *priP,
    bool_t verbose
  )
  {
    /* Compute the average and overall probability: */
    int32_t i, j;
    for (i = 0; i < nf; i++) { E[i] = 0; }
    double sum_wp = 0;
    double sum_w = 0;
    int32_t k;
    for (k = 0; k < nt; k++) 
      { double* Yk = &(Y[k*nf]);
        double wk = (W == NULL ? 1.0 : W[k]);
        double pk = (P == NULL ? 0.5 : (good ? P[k] : 1.0 - P[k]));
        for (i = 0; i < nf; i++) 
          { double yki = Yk[i];
            E[i] += wk*pk*yki;
          }
        sum_wp += wk*pk;
        sum_w += wk;
      }
    sum_wp += 1.0e-300; /* To avoid division by zero. */
    sum_w += 1.0e-300; /* To avoid division by zero. */
    for (i = 0; i < nf; i++) { E[i] /= sum_wp; assert(! isnan(E[i])); }
    double pri = sum_wp/sum_w;      
    assert(! isnan(pri));
    
    /* Clear the {V} matrix (upper half): */
    for (i = 0; i < nf; i++) 
      { double* Vi = &(V[i*nf]);
        for (j = i; j < nf; j++) { Vi[j] = 0; }
      }
    
    /* Accumulate the covariance sums in the {V} matrix (upper half): */
    for (k = 0; k < nt; k++) 
      { double* Yk = &(Y[k*nf]);
        double wk = (W == NULL ? 1.0 : W[k]);
        double pk = (P == NULL ? 0.5 : (good ? P[k] : 1.0 - P[k]));
        for (i = 0; i < nf; i++) 
          { double dki = Yk[i] - E[i];
            double* Vi = &(V[i*nf]);
            for (j = i; j < nf; j++)
              { double dkj = Yk[j] - E[j];
                Vi[j] += wk*pk*dki*dkj;
              }
          }
      }
    
    /* Compute the covariances in the upper half of {V} and copy to the lower half: */
    for (i = 0; i < nf; i++) 
      { double* Vi = &(V[i*nf]);
        for (j = i; j < nf; j++)
          { Vi[j] /= sum_wp; 
            assert(! isnan(Vi[j]));
            if (j > i) { V[j*nf + i] = Vi[j]; }
          }
      }
    
    if (verbose) 
      { char* tag = ((char*[2]){ "bad", "gud" })[good];
        lsq_robust_mf_debug_distr(stderr, tag, nf, E, V, pri);
      }

    if (good)
      { /* Scale up the covariance matrix to compensate the shrink effect: */
        double alpha = 1.1; /* !!! Figure out the right factor !!! */
        rmxn_scale(nf, nf, alpha, V, V);
      }
    
    /* Ensure invertibility of {V} and outlier broader than inlier: */
    double beta = 2.0; /* Difference factor. !!! Figure out !!! */
    demand ((Q != NULL) == (! good), "the {Q} matrix is used only for the bad distribution");
    for (i = 0; i < nf; i++) 
      { double *Vii = &(V[i*nf + i]);
        (*Vii) += 1.0e-200; /* Make sure it is nonzero. */
        if (! good)
          { double Vmin = beta * Q[i*nf + i];
            if ((*Vii) < Vmin) { (*Vii) = Vmin; }
          }
      }
      
    /* !!! Make {V} internal. !!! */

    if (nf > 1)
      { /* Find the Eigenvalue decomposition of {V_gud,V_bad}: */
        double *Vc = rmxn_alloc(nf, nf);
        rmxn_copy(nf, nf, V, Vc);
        double et[nf];
        syei_tridiagonalize(nf, Vc, L, et, K);
        int32_t rank;
        syei_trid_eigen(nf, L, et, K, &rank, 0);
        assert(rank == nf);
        free(Vc);
      }
    else
      { K[0] = 1.0;
        L[0] = V[0];
      }
    
    



    /* Return: */
    if (priP != NULL) { (*priP) = pri; }
  }

void lsq_robust_mf_fudge_covariance_matrix(int32_t nf, double alpha, double V[], double beta, double Q[])
  {
    /* !!! Rethink !!! */
    int32_t i;
  }

double lsq_robust_mf_bayes
  ( int32_t nf, 
    double F[], 
    double P_gud,
    double E_gud[], 
    double K_gud[],
    double L_gud[], 
    double E_bad[], 
    double K_bad[],
    double L_bad[],
    bool_t verbose
  )
  {
    /* Compute the squared discrepancy relative to the good distr: */
    double d_gud[nf]; /* Displacement of {F} from inlier average. */
    rn_sub(nf, F, E_gud, d_gud);
    rmxn_mul_col(nf, nf, K_gud, d_gud, d_gud);
    rn_unweigh(nf, L_gud, d_gud);
    double d2_gud = rn_norm_sqr(nf, d_gud);
    if (d2_gud > 36.0) { return 0.0; }
    double m_gud = rn_prod(nf, L_gud); /* Prob. scaling factor. */
    
    /* Compute the squared discrepancy relative to the bad distr: */
    double d_bad[nf]; /* Displacement of {F} from inlier average. */
    rn_sub(nf, F, E_bad, d_bad);
    rmxn_mul_col(nf, nf, K_bad, d_bad, d_bad);
    rn_unweigh(nf, L_bad, d_bad);
    double d2_bad = rn_norm_sqr(nf, d_bad);
    if (d2_bad > 36.0) { return 0.0; }
    double m_bad = rn_prod(nf, L_bad); /* Prob. scaling factor. */

    /* Apply Bayes's formula: */
    double P_bad = 1 - P_gud;
    double Pc_gud = exp(-d2_gud/2)/m_gud; /* Prob of func values times {sqrt(2*PI)}, assuming good. */
    double Pc_bad = exp(-d2_bad/2)/m_bad; /* Prob of func values times {sqrt(2*PI)}, assuming bad. */ 
    double PP_gud = P_gud*Pc_gud;
    double PP_bad = P_bad*Pc_bad;
    return PP_gud/(PP_gud + PP_bad);
  }

void lsq_robust_mf_debug_distr(FILE *wr, char *tag, int32_t nf, double E[], double V[], double pri)
  {
    fprintf(wr, "--- distribution of independent variables for \"%s\" records ---\n", tag);
    fprintf(wr, "prior prob = %9.4f\n", pri);
    int32_t i, j;
    for (i = 0; i < nf; i++)
      { fprintf(wr, "%4d %+16.8f | ", E[i]);
        for (j = 0; j < nf; j++) { fprintf(wr, " %+20.16f", V[i*nf + j]); }
        fprintf(wr, "\n");
      }
  }     
