/* See rmxn_throw.h. */
/* Last edited on 2024-11-23 18:56:12 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <cmp.h>

#include <rmxn_throw.h>

#define Pr fprintf
#define Er stderr

void rmxn_throw_near_singular_pair
  ( uint32_t n,
    double *A,
    double *B,
    double detRef,
    bool_t isMax
  );
  /* Does {rmxn_throw_almost_singular_pair} with {detMax=detRef} 
    if {isMax} is true, and {rmxn_throw_non_singular_pair} 
    with {detMIn=detRef} if {isMax} is false. */

void rmxn_throw_matrix(uint32_t m, uint32_t n, double *A)
  { for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < n; j++) 
          { A[i*n + j] += dabrandom(-1.0, +1.0); }
      }
  }

void rmxn_rot2(double c, double s, double *x, double *y);
  /* Rotates the vector {(x,y)} by an angle {theta}, given
     {x}, {y}, {c = cos(theta)}, and {s = sin(theta)}. */

void rmxn_throw_ortho(uint32_t n, double *A)
  { /* Start with the identity matrix: */
    rmxn_ident(n, n, A);
    /* Now apply mirroring ops to it: */
    double s[n];
    uint32_t flip = 0;
    for (uint32_t k = 0; k < n; k++)
      { /* Now the first {k} rows and cols of {A} are a random orthonormal {k×k} matrix. */
        double *Ak = &(A[n*k]); /* Row {k} of {A} */
        /* Pick a random direction {v} in {R^{k+1}}, store it in row {k} of {A}: */
        rn_throw_dir(k+1, Ak);
        /* We now apply to rows {0..k-1} a reflection that takes {u_k} to {v}. */
        /* Compute the unit vector {s[0..k]} normal to the bisector of {u_k} and {v}: */
        double s2 = 0;
        for (uint32_t j = 0; j <= k; j++)
          { double sj = Ak[j] - (j == k ? 1 : 0); 
            s[j] = sj; 
            s2 += sj*sj;
          }
        if (s2 != 0)
          { /* Reflection is not trivial. */
            double sm = sqrt(s2);
            for (uint32_t j = 0;  j <= k; j++) { s[j] /= sm; }
            /* Apply reflection along {s} to each row: */
            for (uint32_t i = 0;  i < k; i++)
              { double *Ai = &(A[i*n]);
                (void)rn_mirror(k+1, Ai, s, Ai);
                flip = 1 - flip;
              }
          }
      }
    /* Now flip the first row if needed to keep the determinant positive: */
    if (flip != 0) { for (uint32_t j = 0;  j < n; j++) { A[j] = -A[j]; } }
  }

void rmxn_throw_ortho_complement(uint32_t n, uint32_t ma, double *A, uint32_t mb, double *B)
  {
    demand(ma + mb <= n, "invalid {ma+mb}");
    if (mb == 0) { return; }
    for (uint32_t k = 0;  k < mb; k++)
      { double *Bk = &(B[k*n]); /* Row {k} o f{B}. */
        while (TRUE)
          { /* Pick a random unit vector in {\RR^n}: */
            rn_throw_dir (n, Bk);
            /* Project it onto the space orthogonal to {A} and {B}: */
            for (uint32_t r = 0;  r < ma; r++)
              { double *Ar = &(A[r*n]);
                double s = rn_dot(n, Ar, Bk);
                rn_mix_in(n, -s, Ar, Bk);
              }
            for (uint32_t r = 0;  r < k; r++)
              { double *Br = &(B[r*n]);
                double s = rn_dot(n, Br, Bk);
                rn_mix_in(n, -s, Br, Bk);
              }
            /* Normalize and accept if not too small: */
            double len = rn_dir(n, Bk, Bk);
            if (len >= 1.0e-4) { break; }
          }
      }
  }

void rmxn_throw_directions(uint32_t m, uint32_t n, double U[])
  {
    demand((n >= 2) || (m == n), "invalid {m,n} combination");
    for (uint32_t i = 0; i < m; i++)
      { double *Ui = &(U[i*n]); 
        if (i < n)
          { rn_axis(n, i, Ui); }
        else
          { double minMaxCos = +INF;
            double u[n];
            uint32_t maxTry = 10;
            uint32_t nTry = 0;
            while(TRUE)
              { nTry++;
                /* Generate a random direction {u}: */
                rn_throw_dir (n, u);
                /* Check whether {u} is far enough from previous directions: */
                double maxCos = 0.0; /* Max {fabs(cos(u,Uk))} for previous rows {Uk}. */
                for (uint32_t k = 0; k < i; k++)
                  { double *Uk = &(U[k*n]); 
                    double cos = fabs(rn_dot(n, u, Uk));
                    if (cos > maxCos) { maxCos = cos; }
                  }
                if (maxCos < minMaxCos)
                  { rn_copy(n, u, Ui); 
                    minMaxCos = maxCos;
                    if ((nTry >= maxTry) || (minMaxCos < cos(M_PI/6))) { break; }
                  }
              }
          }
      }
  }
  
void rmxn_throw_almost_singular_pair(uint32_t n, double *A, double *B, double detMax)
  {
    if (n <= 1) { demand(detMax >= 1.0, "detMax must be at least 1 if n<=1"); }
    bool_t isMax = TRUE;
    rmxn_throw_near_singular_pair(n, A, B, detMax, isMax);
  }

void rmxn_throw_non_singular_pair(uint32_t n, double *A, double *B, double detMin)
  {
    if (n <= 1) { demand(detMin <= 1.0, "detMin must be at most 1 if n<=1"); }
    bool_t isMax = FALSE;
    rmxn_throw_near_singular_pair(n, A, B, detMin, isMax);
  }

void rmxn_throw_near_singular_pair(uint32_t n, double *A, double *B, double detRef, bool_t isMax)
  {
    bool_t debug = FALSE;
    if (debug) { Pr(Er, "    > --- %s ---\n", __FUNCTION__); }
      
    if (debug) { Pr(Er, "    n = %d  det%s = %24.16e\n", n, (isMax ? "Max" : "Min"), detRef); }
    
    if (n == 0) { return; }
    
    if (n == 1) { A[0] = 1.0; B[0] = 1.0; return; }
    
    demand(n >= 2, "invalid matrix size");
    
    if (drandom() < 0.5)
      { /* Swap pointers {A,B} to make the result statistically symmetrical */
        double *T = A; A = B; B = T;
      }

    /* Say that a matrix is normal if it has unit RMS element, and
      normalization consists in scaling a matrix by a homogeneous
      positive factor so that it is normalized.
      
      Let {A} be an {n×n} matrix obtained from a normal singular matrix
      {S} by adding to each element a random number in the range {[-eps
      _ +eps]}, for small {eps}, and then normalizing the result. Let
      {B} be the inverse of {A}, normalized.
      
      Tests show that the determinant of {A} will have zero mean and
      deviation {~mulA*(n)*eps}, and that of {B} will have deviation
      {~mulB*(n)*eps^{n-1}}.
      
      Therefore, if we use
      {eps=max(epsA,epsB)} where {epsA=detRef/mulA} and
      {epsB=(detRef/mulB)^{1/(n-1)}}, we have a good chance that 
      both determinants are at least {detRef}, 
      and also a good chance that at least one of them is 
      at most {detRf}. 
      
      Even if the desired condition is not satisfied with that 
      {eps}, at each subsequent try we should decrease it if {isMax} is true,
      and increase it if {isMax} is false. */
    
    demand(detRef >= 1.0e-13, "{detRef} is too small");
    bool_t ABok;
    double mulA = 0.5; /* !!! Fix !!! */
    double epsA = detRef/mulA;
    double mulB = 1.0; /* !!! Fix !!! */
    double epsB = pow(detRef/mulB, 1.0/(n-1));
    double eps = fmax(epsA, epsB);
    uint32_t iter = 0;
    do
      { if (debug) { Pr(Er, "      iteration %d\n", iter); }
        
        if (debug) { Pr(Er, "      filling {A} with a singular array ...\n"); }
        rmxn_throw_singular(n, A);
        
        if (debug) { Pr(Er, "      normalizing {A} to unit RMS elem ...\n"); }
        double Senorm = rmxn_norm(n, n, A)/n; /* RMS of {A} elems. */
        if (debug) { Pr(Er, "      norm(A)/n = %24.16e\n", Senorm); }
        assert(isfinite(Senorm) && (Senorm > 1.0e-200)); /* Almost certain. */
        rmxn_scale(n, n, 1.0/Senorm, A, A);
            
        if (debug) { Pr(Er, "      perturbing {A} with eps = %24.16e\n", eps); }
        rmxn_perturb_unif(n, n, eps, 0.0, A);

        ABok = TRUE; /* Turns false when the iteration failed. */
        
        double Aenorm;
        if (ABok) 
          { if (debug) { Pr(Er, "      checking if {A} is far enough from zero ...\n"); }
            Aenorm = rmxn_norm(n, n, A)/n; /* RMS of {A} elems. */
            if (debug) { Pr(Er, "      norm(A)/n = %24.16e\n", Aenorm); }
            if ((! isfinite(Aenorm)) || (fabs(Aenorm) < 1.0e-100)) { ABok = FALSE; }
          }
        double detA;
        if (ABok) 
          { if (debug) { Pr(Er, "      normalizing {A} to unit RMS elem ...\n"); }
            rmxn_scale(n,n, 1/Aenorm, A, A);
            detA = rmxn_det(n, A);
            if (debug) { Pr(Er, "      det(A) = %24.16e\n", detA); }
            if (! isMax)
              { /* We want both dets to be at least {detRef}, so: */
                if (debug) { Pr(Er, "      ensuring {det(A) >= detRef} ...\n"); }
                if (fabs(detA) < detRef) { ABok = FALSE; }
              }
          }
        /* Now get the inverse {B} of {A} and make sure it is OK: */
        double Benorm;
        if (ABok)
          { if (debug) { Pr(Er, "      computing {B} the inverse of {A} ...\n"); }
            rmxn_inv(n, A, B);
            Benorm = rmxn_norm(n, n, B)/n;
            if (debug) { Pr(Er, "      norm(B)/n = %24.16e\n", Benorm); }
            if ((! isfinite(Benorm)) || (fabs(Benorm) < 1.0e-100)) { ABok = FALSE; }
          }
        double detB;
        if (ABok)
          { if (debug) { Pr(Er, "      normalizing {B} to unit RMS elem ...\n"); }
            rmxn_scale(n, n, 1/Benorm, B, B);
            detB = rmxn_det(n, B);
            if (debug) { Pr(Er, "      det(B) = %24.16e\n", detB); }
            if (! isMax)
              { /* We want both dets to be at least {detRef}, so: */
                if (debug) { Pr(Er, "      ensuring {det(B) >= detRef} ...\n"); }
                if (fabs(detB) < detRef) { ABok = FALSE; }
              }
          }
        if (ABok && isMax)
          { /* We want at least one determinant to be leq {detRef}, so: */
            if (debug) { Pr(Er, "      ensuring {det(A) <= detRef} or {det(B) <= detRef} ...\n"); }
            if (fmin(fabs(detA), fabs(detB)) > detRef) { ABok = FALSE; }
          }
        iter++;
        demand (iter < 10, "exceeded max iterations");
        eps = (isMax ? eps/2 : eps*2);
        if (debug) { Pr(Er, "\n"); }
      }
    while (! ABok);
  }

void rmxn_throw_LT_matrix(uint32_t m, double *Lmm)
  { for (uint32_t i = 0;  i < m; i++)
      { for (uint32_t j = 0;  j < m; j++) 
          { Lmm[m*i + j] = (j <= i ? 2.0 * drandom() - 1.0 : 0.0); }
      }
  }
  
void rmxn_throw_singular(uint32_t n, double *A)
  { bool_t debug = FALSE;
    demand(n >= 2, "size must be at least {2×2}");

    /* Clear out row {n-1}: */
    double *Alast = &(A[(n-1)*n]);
    for (uint32_t j = 0;  j < n; j++) { Alast[j] = 0.0; }

    /* Fill rows {0..n-2} with random elems, mix into row {n-1}: */
    for (uint32_t i = 0;  i < n-1; i++) 
      { double ri = 2*drandom() - 1;
        for (uint32_t j = 0;  j < n; j++)
          { A[i*n + j] = (2*drandom() - 1); 
            Alast[j] += ri * A[i*n + j];
          }
      }
     double Aenorm = rmxn_norm(n, n, A)/n; /* RMS of each element. */
     if (debug) { Pr(Er, "      Aenorm = %24.16e\n", Aenorm); }
     rmxn_scale(n, n, 1.0/Aenorm, A, A);
  }

