/* See {float_image_hdyn.h}. */
/* Last edited on 2013-10-20 23:06:11 by stolfilocal */

#define _GNU_SOURCE
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
 
#include <bool.h>
#include <jsmath.h>
#include <gauss_elim.h>
#include <affirm.h>
#include <sample_conv.h>
#include <sample_conv_hdyn.h>
#include <float_image.h>
#include <float_image_hdyn.h>

/* INTERNAL PROTOTYPES */

int float_image_hdyn_sample_class(double V, double Vmin, double Vmax);
  /* The classification of a sample with value {V}. Returns {-1} when
    {V <= Vmin} (underexposed pixel), {+1} when {V >= Vmax}
    (overexposed pixel), an 0 otherwise (valid pixel). */

bool_t float_image_hdyn_check_convergence
  ( int NI,
    double new[],
    double old[],
    double tol,
    bool_t log_scale, 
    bool_t verbose,
    char *name
  );
  /* Checks convergence of a parameter vector. Namely, compares
    {new[0..NI-1]} with {old[0..NI-1]}; returns TRUE if the changes
    are all less than {tol} in absolute value. If {logscale} is true,
    compares their logarithms instead. If {verbose} is TRUE, prints
    out the {name} and the largest change seen. */

void float_image_hdyn_compute_avg_corr_and_weight
  ( int c,               /* Channel. */
    float_image_t *imgA, /* First image. */
    double vminA,        /* Underexposed value of {imgA} */
    double vmaxA,        /* Overexposed value of {imgA} */
    double sigmaA,       /* Assumed deviation of noise in {imgA}. */
    float_image_t *imgB, /* Second image. */
    double vminB,        /* Underexposed value of {imgB} */
    double vmaxB,        /* Overexposed value of {imgB} */
    double sigmaB,       /* Assumed deviation of noise in {imgB}. */
    bool_t verbose,      /* TRUE to print out detailed diagnostics. */
    double *MAP,         /* (OUT) Mean value of {imgA} image on common pixels. */
    double *MBP,         /* (OUT) Mean value of {imgB} image on common pixels. */
    double *CP,          /* (OUT) Regression coef of {imgA} rel to {imgB} on common pixels. */
    double *WP           /* (OUT) Fuzzy count of common pixels. */
  );
  /* This procedure assigns to each pixel {x,y} a weight that is 1 iff
    the values of {imgA,imgB} are well away from 0 and 1, assuming
    noise with deviations {sigmaA,sigmaB}, respectively.
    
    Computes the weighted averages {*MAP} and {*MBP} of the values
    {imgA} and {imgB}, respectively. It also returns in {*WP} the sum
    of those weights, and in {*CP} the weighted linear regression
    coefficient {C} for the model {(imgA-MA) = C*(imgB-MB)}, which 
    may be 0, {INF}, or {NAN}. */

void float_image_hdyn_print_system(FILE *wr, int m, int n, double A[], int p, double B[]);
  /* Prints to {wr} the matrix {A} with {m} rows and {n} columns, side 
    by side with the matrix {B} with {m} rows and {p} columns, both linarized by rows. */ 
  
typedef double float_image_hdyn_diff_func_t(int i, int j);
/* Type of function that defines the difference between two variables {X[i]} and {X[j]}.
   See {float_image_hdyn_solve_diff_eqs}. */
    
void float_image_hdyn_solve_diff_eqs
  ( int N, 
    double W[], 
    float_image_hdyn_diff_func_t *diff, 
    bool_t verbose,
    double X[]
  );
  /* Solves a system of {N*(N-1)/2} linear equations on {N} unknowns
    {X[0..N-1]}. The equations are {X[i] - X[j] = diff(i,j)} for all
    {i} in {0..N-1} and all {j} in {0..N-1} with {i != j}, where
    {diff} is an arbitrary function.
    
    The system is underscontrained, since adding a constant to every
    {X[i]} preserves every difference {X[i] - X[j]}. To remove the
    ambiguity, the procedure adds the requirement that the sum of all
    {X[i]} is zero.
    
    The system is also overconstrained since there are {N*(N-1)}
    equations and only {N-1} degrees of freedom. Therefore it is
    solved by weighted least squares. The weight of each equation is
    {W[i,j]} where {W} is an {N} by {N} matrix of finite non-negative
    numbers, stored by rows in a vector with {NI^2} elements. In
    particular, if {W[i,j]} is zero, the equation {i,j} is ignored.
    The equation is ignored also if {diff(i,j)} returns {NAN}. If too
    many equations are ignored, the procedure may fail. */

/* IMPLEMENTATIONS */

void float_image_hdyn_estimate_gains_offsets_sigmas_values
  ( int c,                /* Channel. */
    int NI,               /* Number of input images. */
    float_image_t *img[], /* Input images. */
    double vmin[],        /* Value of underexposed pixels. */
    double vmax[],        /* Value of overexposed pixels. */
    bool_t verbose,       /* TRUE for diagnostics. */
    double gain[],        /* (OUT) Estimated gain of each image. */
    double offset[],      /* (OUT) Estimated value offset of each image. */
    double sigma[],       /* (OUT) Estimated noise deviation. */
    float_image_t *omg    /* (OUT) Estimated true values. */
  )
  {
    /* Initial guess for {sigmas}: */
    int i;
    for (i = 0; i < NI; i++) { sigma[i] = 0.003; /* About right for 8-bit quantization */ }
    int max_iterations = 1;
    double tol_gain = 0.0001;
    double tol_offset = 0.0001;
    double tol_sigma = 0.0001;
    bool_t converged = FALSE;
    double gain_old[NI];
    double offset_old[NI];
    double sigma_old[NI];
    int iter = 0;
    while (TRUE)
      { if (verbose)
          { /* Print current state: */
            fprintf(stderr, "iteration %d\n", iter);
            if (iter > 0)
              { fprintf(stderr, "  gain =   ");
                for (i = 0; i < NI; i++) { fprintf(stderr, " %10.6f", gain[i]); }
                fprintf(stderr, "\n");
                fprintf(stderr, "  offset = ");
                for (i = 0; i < NI; i++) { fprintf(stderr, " %10.6f", offset[i]); }
                fprintf(stderr, "\n");
              }
            fprintf(stderr, "  sigma =  ");
            for (i = 0; i < NI; i++) { fprintf(stderr, " %10.6f", sigma[i]); }
            fprintf(stderr, "\n");
          }
        /* Check for convergence: */
        if ((iter >= max_iterations) || converged) { break; }
        /* Save current {gain,sigma}: */
        for (i = 0; i < NI; i++) { gain_old[i] = gain[i]; offset_old[i] = offset[i]; sigma_old[i] = sigma[i]; }
        /* Re-estimate {gain} assuming current {sigma}: */
        float_image_hdyn_estimate_gains_offsets(c, NI, img, vmin, vmax, sigma, verbose, /*OUT*/ gain, offset);
        /* Estimate {Y} values assuming {gain,offset,sigma}: */
        float_image_hdyn_estimate_values(c, NI, img, vmin, vmax, gain, offset, sigma, verbose, /*OUT*/ omg);
        /* Re-estimate {sigma} assuming current {gain,offset} and values: */
        float_image_hdyn_estimate_sigmas(c, NI, img, vmin, vmax, gain, offset, omg, verbose, /*OUT*/ sigma);
        /* Make sure that {Vmax - Vmin >= 2*sigma}: */
        for (i = 0; i < NI; i++)
          { double sigma_max = (vmax[i] - vmin[i])/2.0001;  /* Max valid sigma. */
            if (sigma[i] > sigma_max) { sigma[i] = sigma_max; }
          }
        /* Did one more iteration: */
        iter++;
        /* Find max change in {gain,offset,sigma}: */
        if (verbose){ fprintf(stderr, "  max changes:"); }
        converged = TRUE;
        converged &= float_image_hdyn_check_convergence(NI, gain,   gain_old,   tol_gain,   TRUE,  verbose, "gain");
        converged &= float_image_hdyn_check_convergence(NI, offset, offset_old, tol_offset, FALSE, verbose, "offset");
        converged &= float_image_hdyn_check_convergence(NI, sigma,  sigma_old,  tol_sigma,  FALSE, verbose, "sigma");
      }
    if (! converged) { fprintf(stderr, "  failed to converge in %d iterations\n", iter); }
  }

bool_t float_image_hdyn_check_convergence
  ( int NI,
    double new[],
    double old[],
    double tol,
    bool_t log_scale, 
    bool_t verbose,
    char *name
  )
  {
    double max_d = 0.0;
    int i;
    for (i = 0; i < NI; i++)
      { double di;
        if (log_scale)
          { di = log(new[i]) - log(old[i]); }
        else
          { di = new[i] - old[i]; }
        if (fabs(di) > fabs(max_d)) { max_d = di; }
      }
    if (verbose)
      { if (log_scale)
          { fprintf(stderr, "  in log(%s) = ", name); }
        else
          { fprintf(stderr, "  in %s = ", name); }
        fprintf(stderr, "%+12.8f\n", max_d);
      }
    return fabs(max_d) <= tol;
  }

int float_image_hdyn_sample_class(double V, double Vmin, double Vmax)
  { 
    if (V <= Vmin)
      { return -1; }
    else if (V >= Vmax)
      { return +1; }
    else
      { return 0; }
  }

void float_image_hdyn_estimate_gains_offsets
  ( int c,                /* Channel. */
    int NI,               /* Number of input images. */
    float_image_t *img[], /* Input images. */
    double vmin[],        /* Value of underexposed pixels. */
    double vmax[],        /* Value of overexposed pixels. */
    double sigma[],       /* Assumed noise deviation of each image. */
    bool_t verbose,       /* TRUE for diagnostics. */
    double gain[],        /* (OUT) Estimated gain of each image. */
    double offset[]       /* (OUT) Estimated value offset of each image. */
  )
  {
    bool_t debug = FALSE;
    
    /* For each pair of distinct images
      {i,j}, compute the weighted average {M[i,j]} of the values of image {i},
      the linear regression coefficient {C[i,j]} of image {i} versus image {j},
      and a fuzzy pixel count {W[i,j]}, over all pixels that are valid in 
      both images.  The coefficient {C[i,j]} assumes the equation
      {(V[i] - M[i,j]) = C[i,j]*(V[j] - M[j,i]).
    */
    
    /* Matrices for image pair statistics: */
    if (verbose) { fprintf(stderr, "collecting image statistics...\n"); }
    int NI2 = NI*NI;
    double *M = notnull(malloc(NI2*sizeof(double)), "no mem"); /* Average values. */
    double *C = notnull(malloc(NI2*sizeof(double)), "no mem"); /* Correlation coefficients. */
    double *W = notnull(malloc(NI2*sizeof(double)), "no mem"); /* Weights. */
    int i, j;
    for (i = 0; i < NI; i++) 
      { for (j = 0; j <= i; j++) 
          { double Wij, Cij, Cji, Mij, Mji;
            float_image_hdyn_compute_avg_corr_and_weight
              ( c, 
                img[i], vmin[i], vmax[i], sigma[i], 
                img[j], vmin[j], vmax[j], sigma[j], 
                debug & (i == 2) && (j == 0),
                &Mij, &Mji, &Cij, &Wij
              );
            if (! isnan(Cij))
              { /* Ensure {C[i,j]} is non-negative: */
                if (Cij < -0.00001)
                  { fprintf(stderr, "correlation = %24.15e for images %d %d\n", Cij, i, j);
                    /* assert(FALSE); */
                  }
                if (Cij < 0) { Cij = 0; }
                Cji = (Wij == 0 ? 0 : 1/Cij);
              }
            else
              { Cji = NAN; }
            if (j != i)
              { int kij = i*NI + j;
                C[kij] = Cij;
                M[kij] = Mij;
                W[kij] = Wij;
                /* Fill the upper half of the matrices {C,W}: */
                int kji = j*NI + i;
                C[kji] = Cji;
                M[kji] = Mji; 
                W[kji] = W[kij];
              }
            else
              { int kii = i*NI + i;
                /* Consistency checks: */
                if (fabs(log(Cij)) > log(1.00001))
                  { fprintf(stderr, "self-correlation = %24.15e for image %d\n", Cij, i); }
                if (fabs(Mij - Mji) > 1.0e-10*(fabs(Mij)+fabs(Mji)))
                  { fprintf(stderr, "discrepant means = %24.15e %24.15e for image %d\n", Mij, Mji, i); }
                C[kii] = 1.0;
                M[kii] = Mij;
                W[kii] = Wij;
              }
          }
      }
    if (verbose)
      { fprintf(stderr, "Correlations, means, and weights for image pairs:\n");
        
        fprintf(stderr, "%3s", "");
        for (j = 0; j < NI; j++) { fprintf(stderr, " %12d", j); }
        fprintf(stderr, "\n");
        
        for (i = 0; i < NI; i++) 
          { fprintf(stderr, "%3d", i);
            for (j = 0; j < NI; j++) { fprintf(stderr, " %+12.6f", C[i*NI+j]); }
            fprintf(stderr, "\n");
            
            fprintf(stderr, "%3s", "");
            for (j = 0; j < NI; j++) { fprintf(stderr, " %12.6f", M[i*NI+j]); }
            fprintf(stderr, "\n");
            
            fprintf(stderr, "%3s", "");
            for (j = 0; j < NI; j++) { fprintf(stderr, " %12.1f", W[i*NI+j]); }
            fprintf(stderr, "\n");
            
            fprintf(stderr, "\n");
          }
          
        fprintf(stderr, "\n");
      }
            
    /* Let {g[i]} be the gain {dV/dY} of each image {i}. For each pair
      of distinct images {i,j}, we assume a basic equation {g[i]/g[j] =
      C[i,j]}, which in log scale become {X[i] - X[j] = log(C[i,j])}
      where {X[i] = log(g[i]). We solve these equations by weighted
      least squares, where the weight of each equation is 
      proportional to {W[i,j]}.
      
      The basic equations are indeterminate by a common additive term.
      We add the linear constraint {sum{X[i]} = 0} to remove that
      indeterminacy. */
      
    auto double diff_gains(int i, int j);
    
    double diff_gains(int i, int j)
      { assert(i != j);
        double Cij = C[i*NI + j];
        if (isfinite(Cij) && (Cij != 0))
          { return log(Cij); }
        else
          { return NAN; }
      }
    
    double *X = notnull(malloc(NI*sizeof(double)), "no mem");
    if (verbose) { fprintf(stderr, "solving the gain equations...\n"); }
    float_image_hdyn_solve_diff_eqs(NI, W, &diff_gains, verbose, X);

    /* Recover the raw gains {dV/dY}: */
    if (verbose) { fprintf(stderr, "extracting the raw gains from the solution...\n"); }
    double min_valid_gain = 1.0e-4; /* Min allowed {gain[i]}. */
    double max_valid_gain = 1.0e+4; /* Max allowed {gain[i]}. */
    double log_min_valid_gain = log(min_valid_gain);
    double log_max_valid_gain = log(max_valid_gain);
    double max_gain = 0.0; /* Max gain among all images: */
    for (i = 0; i < NI; i++)
      { assert(! isnan(X[i]));
        if (X[i] < log_min_valid_gain)
          { fprintf(stderr, "gain too low = exp(%24.15e) for image %d\n", X[i], i);
            gain[i] = 0.0;
          }
        else if (X[i] > log_max_valid_gain)
          { fprintf(stderr, "gain too high = exp(%24.15e) for image %d\n", X[i], i);
            gain[i] = +INF;
          }
        else
          { gain[i] = exp(X[i]);
            if (gain[i] > max_gain) { max_gain = gain[i]; }
          }
      }
    if (verbose) 
      { fprintf(stderr, "raw gains     = ");
        for (i = 0; i < NI; i++) { fprintf(stderr, " %12.6f", gain[i]); }
        fprintf(stderr, "\n");
      }
    
    /* Now we estimate the {offset}s.  For each pair {i,j}
      we have an equation {K[i] - K[j] = M[i,j]/g[i] - M[j,i]/g[j]}
      where {K[i] = offset[i]/gain[i]}. We solve these equations too 
      by weighted least squares, where the weight of each equation is 
      proportional to {W[i,j]}. 
      
      Again, the basic equations are indeterminate by a common additive term.
      We add the linear constraint {sum{K[i]} = 0} to remove that
      indeterminacy. */
    
    auto double diff_offsets(int i, int j);
    
    double diff_offsets(int i, int j)
      { assert(i != j);
        double gi = gain[i];
        double gj = gain[j];
        if (isfinite(gi) && (gi != 0) && isfinite(gj) && (gj != 0))
          { double Mij = M[i*NI + j];
            double Mji = M[j*NI + i];
            if (isfinite(Mij) && isfinite(Mji))
              { return Mij/gi - Mji/gj; }
            else
              { return NAN; }
          }
        else
          { return NAN; }
      }
    
    if (verbose) { fprintf(stderr, "solving the offset equations...\n"); }
    float_image_hdyn_solve_diff_eqs(NI, W, &diff_offsets, verbose, X);
    if (verbose) 
      { fprintf(stderr, "raw Y offsets = ");
        for (i = 0; i < NI; i++) { fprintf(stderr, " %+12.6f", X[i]); }
        fprintf(stderr, "\n");
      }
    
    /* Recover the raw value offsets: */
    for (i = 0; i < NI; i++)
      { assert(! isnan(X[i]));
        if (isfinite(gain[i]) && (gain[i] != 0))
          { offset[i] = X[i]*gain[i]; }
        else
          { offset[i] = 0.0; }
      }
    if (verbose) 
      { fprintf(stderr, "raw V offsets = ");
        for (i = 0; i < NI; i++) { fprintf(stderr, " %+12.6f", offset[i]); }
        fprintf(stderr, "\n");
      }

    /* Compute the range of luminosities implied by offsets and gains: */
    double Ymin_raw = +INF; /* Min raw luminosity among all images: */
    double Ymax_raw = -INF; /* Max raw luminosity among all images: */
    if (verbose) { fprintf(stderr, "raw Y ranges\n"); }
    for (i = 0; i < NI; i++)
      { assert(! isnan(X[i]));
        if (isfinite(gain[i]) && (gain[i] != 0))
          { /* Compute raw luminosity range {[ylo_yhi]} implied by valid values in image {i}: */
            double Ylo = (vmin[i] - offset[i])/gain[i];
            double Yhi = (vmax[i] - offset[i])/gain[i];
            if (verbose) { fprintf(stderr, "  %2d [ %+12.6f _ %+12.6f ]\n", i, Ylo, Yhi); }
            if (Ylo < Ymin_raw) { Ymin_raw = Ylo; }
            if (Yhi > Ymax_raw) { Ymax_raw = Yhi; }
          }
      }
    if (verbose) { fprintf(stderr, "raw Y range =   [ %+12.6f _ %+12.6f ]\n", Ymin_raw, Ymax_raw); }
      
    /* Adjust gains and offsets so that the valid lum range is about {[0_1]}: */
    double delta = Ymax_raw - Ymin_raw; /* Gain adjustment factor. */
    double theta = Ymin_raw;   /* Offset adjustment shift.*/
    for (i = 0; i < NI; i++) 
      { if (isfinite(gain[i]) && (gain[i] != 0))
          { offset[i] = offset[i] + gain[i]*theta;
            gain[i] = gain[i]*delta;
          }
      }

    if (verbose) 
      { fprintf(stderr, "fixed gains   = ");
        for (i = 0; i < NI; i++) { fprintf(stderr, " %12.6f", gain[i]); }
        fprintf(stderr, "\n");
        
        fprintf(stderr, "fixed offsets = ");
        for (i = 0; i < NI; i++) { fprintf(stderr, " %+12.6f", offset[i]); }
        fprintf(stderr, "\n");
      }
    
    free(X);
    free(W);
    free(C);
    free(M);
  }
    
void float_image_hdyn_solve_diff_eqs
  ( int N, 
    double W[], 
    float_image_hdyn_diff_func_t *diff, 
    bool_t verbose,
    double X[]
  )
  {  
    /* Assemble the least squares system: */
    if (verbose) { fprintf(stderr, "assembling the least squares system...\n"); }
    double *A = notnull(malloc(N*N*sizeof(double)), "no mem");
    double *B = notnull(malloc(N*sizeof(double)), "no mem");
    /* Clear the arrays: */
    int i, j;
    for (i = 0; i < N; i++)
      { B[i] = 0;
        for (j = 0; j < N; j++) { A[i*N + j] = 0; }
      }
    /* Accumulate the spring force terms: */
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { if (i != j) 
            { int kii = i*N + i;
              int kjj = j*N + j;
              int kij = i*N + j;
              int kji = j*N + i;
              double Wij = W[kij];
              demand(isfinite(Wij), "invalid weight");
              demand(Wij >= 0, "negative weight");
              if (Wij != 0)
                { double Dij = diff(i, j);
                  if (! isnan(Dij))
                    { demand(isfinite(Dij), "invalid difference");
                      A[kii] += +Wij;
                      A[kjj] += +Wij;
                      A[kij] += -Wij;
                      A[kji] += -Wij;
                      B[i] += +Wij*Dij;
                      B[j] += -Wij*Dij;
                    }
                }
            }
          }
      }
    /* Normalize each equation so that the diagonal elements are 1, and add the zero-sum force: */
    for (i = 0; i < N; i++)
      { int kii = i*N + i;
        double S = A[kii];
        /* If {A[i,i]} is zero, the whole row must be zero. */
        if (S == 0) { S = 1.0; }
        for (j = 0; j < N; j++) 
          { int kij = i*N + j;
            A[kij] = (i == j ? 1.0 : A[kij]/S) + 1.0;
          }
        B[i] = B[i] / S;
      }
    
    /* Solve the system: */
    if (verbose) { fprintf(stderr, "solving the linear system...\n"); }
    gsel_solve(N, N, A, 1, B, X, 0.0);
    
    free(B);
    free(A);
  }

void float_image_hdyn_print_system(FILE *wr, int m, int n, double A[], int p, double B[])
  {
    int i, j, k;
    fprintf(wr, "\n");

    fprintf(wr, "%3s", "");
    for (j = 0; j < n; j++) { fprintf(wr, " %16d", j); }
    fprintf(wr, "%3s", "");
    for (k = 0; k < p; k++) { fprintf(wr, " %16d", k); }
    fprintf(wr, "\n");

    for (i = 0; i < m; i++) 
      { fprintf(wr, "%3d", i);
        for (j = 0; j < n; j++) { fprintf(wr, " %+16.6f", A[i*n+j]); }
        fprintf(wr, "%3s", "");
        for (k = 0; k < p; k++) { fprintf(wr, " %+16.6f", B[i*p+k]); }
        fprintf(wr, "\n");
      }

    fprintf(wr, "\n");
  }

void float_image_hdyn_estimate_values
  ( int c,                /* Channel. */
    int NI,               /* Number of input images. */
    float_image_t *img[], /* Input images. */
    double vmin[],        /* Value of underexposed pixels. */
    double vmax[],        /* Value of overexposed pixels. */
    double gain[],        /* Assumed gain of each image. */
    double offset[],      /* Assumed value offset of each image. */
    double sigma[],       /* Assumed noise deviation of each image. */
    bool_t verbose,       /* TRUE for diagnostics. */
    float_image_t *omg    /* (OUT) Estimated luminances. */
  )
  {
    if (verbose) { fprintf(stderr, "estimating the true image values...\n"); }
    /* Get image dimensions: */
    int NX = (int)img[0]->sz[1];
    int NY = (int)img[0]->sz[2];
    
    /* Scan pixels: */
    int x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { int nlo = 0; /* Number of images where this pixel is underexposed. */
            int nhi = 0; /* Number of images where this pixel is overexposed. */
            int nin = 0; /* Number of images where this pixel is well-exposed. */
            int nun = 0; /* Number of bad images. */
            double sum_W = 0.0;
            double sum_WY = 0.0;
            int i;
            for (i = 0; i < NI; i++) 
              { if (isfinite(gain[i]) && (gain[i] != 0))
                  { double Vi = float_image_get_sample(img[i], c, x, y);
                    double Yi = (Vi - offset[i])/gain[i];
                    double sigYi = sigma[i]/gain[i];
                    int class = float_image_hdyn_sample_class(Vi, vmin[i], vmax[i]);
                    if (class < 0) 
                      { nlo++; }
                    else if (class > 0)
                      { nhi++; }
                    else
                      { nin++;
                        double Vrel = (Vi - vmin[i])/(vmax[i] - vmin[i]);
                        double W_val = 0.5*(1 - cos(2*M_PI*Vrel));
                        double W_sig = 1/(sigYi*sigYi);
                        double Wi = W_val*W_sig;
                        sum_WY += Wi*Yi;
                        sum_W += Wi;
                      }
                  }
                else
                  { nun++; }
              }
            double Y;
            if (nhi >= NI)
              { /* Seems overexposed in all good images: */
                Y = +INF;
              }
            else if (nlo >= NI)
              { /* Seems underexposed in all good images: */
                Y = -INF; 
              }
            else if (nin == 0)
              { /* Seems either under- or over-exposed in all good images: */
                Y = NAN; 
              }
            else if ((nun == NI) || (sum_W == 0))
              { /* There is no good image: */
                Y = NAN;
              }
            else
              { Y = sum_WY/sum_W; }
            float_image_set_sample(omg, c, x, y, (float)Y);
          }
      }
  }
  
void float_image_hdyn_estimate_sigmas
  ( int c,                /* Channel. */
    int NI,               /* Number of input images. */
    float_image_t *img[], /* Input images. */
    double vmin[],        /* Value of underexposed pixels. */
    double vmax[],        /* Value of overexposed pixels. */
    double gain[],        /* Assumed gain of each image. */
    double offset[],      /* Assumed value offset of each image. */
    float_image_t *omg,   /* Assumed luminances. */
    bool_t verbose,       /* TRUE for diagnostics. */
    double sigma[]        /* (OUT) Estimated noise deviation. */
  )
  {
    /* Get image dimensions: */
    int NX = (int)img[0]->sz[1];
    int NY = (int)img[0]->sz[2];
    int i;
    if (verbose) { fprintf(stderr, "estimating the noise deviations...\n"); }
    for (i = 0; i < NI; i++) 
      { if (verbose) { fprintf(stderr, "  sigma[%2d] =", i); }
        double Nvalid = NAN;
        double sum_W = 1.0e-200;
        double sum_WD2 = 0;
        int x, y;
        for (y = 0; y < NY; y++)
          { for (x = 0; x < NX; x++)
              { double V_obs = float_image_get_sample(img[i], c, x, y);
                int class = float_image_hdyn_sample_class(V_obs, vmin[i], vmax[i]);
                if (class == 0)
                  { double W = 1.0;
                    double Y = float_image_get_sample(omg, c, x, y);
                    double V_exp = gain[i]*Y + offset[i];
                    double D = V_obs - V_exp;
                    sum_WD2 += W*D*D;
                    sum_W += W;
                  }
              }
          }
        sigma[i] = sqrt(sum_WD2/sum_W);
        fprintf(stderr, " %10.6f", sigma[i]);
        Nvalid = sum_W;
        if (verbose) { fprintf(stderr, " ~%.2f%% valid pixels", 100*Nvalid/(NX*NY)); }
        fprintf(stderr, "\n");
      }
  }

void float_image_hdyn_compute_avg_corr_and_weight
  ( int c,               /* Channel. */
    float_image_t *imgA, /* First image. */
    double vminA,        /* Underexposed value of {imgA} */
    double vmaxA,        /* Overexposed value of {imgA} */
    double sigmaA,       /* Assumed deviation of noise in {imgA}. */
    float_image_t *imgB, /* Second image. */
    double vminB,        /* Underexposed value of {imgB} */
    double vmaxB,        /* Overexposed value of {imgB} */
    double sigmaB,       /* Assumed deviation of noise in {imgB}. */
    bool_t verbose,      /* TRUE to print out detailed diagnostics. */
    double *MAP,         /* (OUT) Mean value of {imgA} image on common pixels. */
    double *MBP,         /* (OUT) Mean value of {imgB} image on common pixels. */
    double *CP,          /* (OUT) Regression coef of {imgA} rel to {imgB} on common pixels. */
    double *WP           /* (OUT) Fuzzy count of common pixels. */
  )
  {
    /* Get image dimensions: */
    int NX = (int)imgA->sz[1]; demand(((int)imgB->sz[1]) == NX, "incompat cols");
    int NY = (int)imgA->sz[2]; demand(((int)imgB->sz[2]) == NY, "incompat rows");

    /* Compute weight and means: */
    double sum_W = 1.0e-200;
    double sum_WVA = 0;
    double sum_WVB = 0;
    int x, y;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { double VA = float_image_get_sample(imgA, c, x, y);
            assert(VA >= 0);
            int classA = float_image_hdyn_sample_class(VA, vminA, vmaxA);

            double VB = float_image_get_sample(imgB, c, x, y);
            assert(VB >= 0);
            int classB = float_image_hdyn_sample_class(VB, vminB, vmaxB);
            
            if ((classA == 0) && (classB == 0))
              { double WAB = 1.0;
                sum_WVA += WAB*VA;
                sum_WVB += WAB*VB;
                sum_W += WAB;
              }
          }
      }
    double MA = sum_WVA/sum_W;
    double MB = sum_WVB/sum_W;
    double W = sum_W;

    /* Compute correlation: */
    double sum_WdVAdVB = 0;
    double sum_WdVB2 = 0;
    for (y = 0; y < NY; y++)
      { for (x = 0; x < NX; x++)
          { double VA = float_image_get_sample(imgA, c, x, y);
            int classA = float_image_hdyn_sample_class(VA, vminA, vmaxA);

            double VB = float_image_get_sample(imgB, c, x, y);
            int classB = float_image_hdyn_sample_class(VB, vminB, vmaxB);
            
            if ((classA == 0) && (classB == 0))
              { double WAB = 1.0;
                double dVA = VA - MA;
                double dVB = VB - MB;
                if (verbose) { fprintf(stderr, "  %5d %5d %+12.8f %+12.8f\n", x, y, dVA, dVB); }
                sum_WdVAdVB += WAB*dVA*dVB;
                sum_WdVB2 += WAB*dVB*dVB;
              }
          }
      }
    if (verbose) { fprintf(stderr, "  sums AB = %+12.8f BB = %+12.8f\n", sum_WdVAdVB, sum_WdVB2); }
    double C = sum_WdVAdVB/sum_WdVB2; /* May be 0, {INF} or {NAN}. */
    
    (*MAP) = MA;
    (*MBP) = MB;
    (*CP) = C;
    (*WP) = W;
  }
