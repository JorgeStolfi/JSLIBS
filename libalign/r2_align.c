/* See {r2_align.h}. */
/* Last edited on 2023-04-01 04:29:59 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <sym_eigen.h>
#include <affirm.h>

#include <r2_align.h>

/* INTERNAL PROTOTYPES */

bool_t r2_align_coord_is_variable (r2_t arad[], int32_t i, int32_t j);
  /* Returns {TRUE} iff {arad[i].c[j]} is nonzero. */

i2_t r2_align_count_variable_coords (int32_t ni, r2_t arad[]);
  /* Returns a integer pair {nv} such that {nv.c[j]} is the number of
    indices {i} such that {arad[i].c[j]} is not zero. */

void r2_align_plot_mismatch_lines
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    double step,  /* Grid step. */
    r2_t u0[],    /* Direction across lines. */
    int32_t n0,   /* Number of lines on each side of center. */
    r2_t u1[],    /* Direction along lines. */
    int32_t n1,   /* Max number of points along each line on each side of center. */
    r2_align_mismatch_t *F2
  );
  /* Similar to {r2_align_plot_mismatch}, but samples the function only
    for {n0} lines spaced by {step} along direction {u0}, Each line will
    have {n1} points spaced by {step} in direction {u1}. For each
    sampled point {p}, writes a line with two independent coordinates
    {s0,s1} along {u0,u1}, and the value of {F2(ni,p)} . */

/* ---------------------------------------------------------------------- */

void r2_align_throw_ortho_disp_vector(int32_t ni, r2_t arad[], int32_t k, r2_t U[]);
  /* Assumes that the rows of {U} are {k} delta vectors with {ni} points
    each, orthonormal, balanced, and conformal to {arad}, and {U} has space
    for at least one more row. Stores into row {k} of {U} another alignment
    vector that is normalized, balanced, conformal to {arad}, and
    orthogonal to all of them.
    
    Specifically, assumes that the delta vector basis element {u[kb]} with index {kb} is
    points {u[kb][i] = U[kb*ni + i]} for {i} in {0..ni-1}. 
    
    The parameter {ni} must be at least 2, and {k} must be less than {nd}
    where {nv} is the number of non-zero coordinates in {arad}. */


/* IMPLEMENTATIONS */

i2_t r2_align_count_variable_coords(int32_t ni, r2_t arad[])
  { 
    i2_t nv = (i2_t){{ 0, 0 }};
    for (int32_t i = 0; i < ni; i++) 
      { for (int32_t j = 0; j < 2; j++)
          { double arij = arad[i].c[j];
            demand (isfinite(arij) && (arij >= 0), "invalid {arad}");
            if (arij != 0) { nv.c[j]++; }
          }
      }
    return nv;
  }

int32_t r2_align_count_degrees_of_freedom(i2_t nv, bool_t bal)
  { int32_t nd;
    if (bal) 
      { nd = (nv.c[0] >= 1 ? nv.c[0]-1 : 0) + (nv.c[1] >= 1 ? nv.c[1]-1 : 0); }
    else
      { nd = nv.c[0] + nv.c[1]; }
    return nd;
  }

void r2_align_compute_search_ellipsoid
  ( int32_t ni,
    r2_t arad[],
    bool_t bal, 
    int32_t *nd_P,
    r2_t U[], 
    double urad[]
  )
  { 
    bool_t debug = TRUE;
    
    if (debug) { fprintf(stderr, "--> entering %s\n", __FUNCTION__); }

    demand(ni >= 0, "invalid {ni}");
    i2_t nv = r2_align_count_variable_coords(ni, arad); /* Variable coords count per axis of {\R^2}. */
    int32_t nh = nv.c[0] + nv.c[1]; /* Total count of variable coords of {\RC}. */
    int32_t nd == r2_align_count_degrees_of_freedom(nv, bal);
    if (nd == 0) { return; }
   
    /* If there are no balancing constraints, the rows of the search
      basis {U} are simply the {nh} canonical unit vectors of {\RC}
      corresponding to the variable coordinates (those with positive
      {arad}), and {urad[0..nd-1]} are the non-zero elelements of
      {arad[0..ni-1].c[0..1]}. Otherwise we must cut the unconstrained
      ellipsoid defined by {arad} with the balancing equations. */
    
    if (bal && (nd < nh))
      { if (debug) { fprintf(stderr, "  ... lienarizing {arad}\n"); } 
        double hrad[2*ni]; /* The radii {arad[0..ni-1].c[0..1]}, linearized. */
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++)
              { hrad[2*i + j] = arad[i].c[j]; }
          }

        if (debug) { fprintf(stderr, "  ... computing the blancing constraints\n"); } 
        int32_t nc_bal = nh - nd;   /* Number of balancing constraints. */
        double A[nc_bal*2*ni];  /* Constraint equation coefficients on linearized {\RC} coords. */
        int32_t kc = 0;   /* Next free row of {A}. */
        for (int32_t j = 0; j < 2; j++)
          { if (nv.c[j] >= 2)
              { /* Append the the balancing constraint on coordinate {j} of {\RR^2}: */
                assert(kc < nc_bal);
                double *Ak = &(A[kc*2*ni]);
                for (int32_t i = 0; i < ni; i++)
                  { double arij = arad[i].c[j];
                    if (arij != 0) { Ak[2*i + j] = 1.0; }
                  }
                kc++;
              }
          }
        assert(kc == nc_bal);
        
        if (debug) { fprintf(stderr, "  ... orthonormalizing balancing and conform constraints\n"); } 
        int32_t nc; /* Number of independent constraints. */
        double *C = NULL; /* Orthonormalized constraints over {\RR^{2*ni}}. */
        rmxn_ellipsoid_normalize_constraints(2*ni, hrad, nc_bal, A, verbose_norm, &nc, &C);
        int32_t nc_zer = 2*ni - nh; /* Number of conformation constraints. */
        assert(nc == nc_bal + nc_zer); /* Balancing and conform constraints should be indep. */
          
        /* Cut the ellipsoid in the linearized space {\RR^{2*ni}}: */
        double V[nd*2*ni]; /* Linearized search basis {U}. */
        rmxn_ellipsoid_cut(2*ni, hrad, nc, C, nd, V, urad);
        
        /* Repack basis {V} as alignment vectors in {\RC}: */
        for (int32_t kd = 0; kd < nd; kd++)
          { r2_t *Uk = &(U[k*ni]);
            double Vk = &(V[k*2*ni]);
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++)
                  { double Vkij = Vk[2*i + j];
                    /* All {U} vectors must be conformal: */
                    if (arad[i].c[j] == 0) { assert(Vkij == 0); }
                    Uk[i].c[j] = Vkij; 
                  }
              }
          }
      }
    else
      { int32_t kd = 0; /* Count nonzero elements of {arad}. */
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++)
              { double arij = arad[i].c[j];
                if (arij != 0) 
                  { assert(kd < nd);
                    /* Store into row {kd} the {i,j} canonical vector of {\RC}: */
                    r2_t *Uk = &(U[kd*ni]);
                    for (int32_t i1 = 0; i1 < ni; i1++)
                      { for (int32_t j1 = 0; j1 < 2; j1++)
                         { Uk[i1].c[j1] = 0.0; }
                      }
                    Uk[i].c[j] = 1.0;
                    urad[kd] = arij;
                    kd++;
                  }
              }
          }
      }
    assert(kd = nd;
  }
  
        
        
    

    if (debug) { fprintf(stderr, "... finding an orthonormal basis {H} for {\\RU} ...\n"); }

    r2_t H[nd*ni];
    for (int32_t k = 0; k < nd; k++)
      { r2_t *hk = &(H[k*ni]);
        r2_align_throw_ortho_disp_vector(ni, arad, k, H);
        if (debug) { r2_align_print_vector(stderr, ni, "h", k, hk, FALSE); }
      }
      
    if (debug) { fprintf(stderr, "... computing the metric matrix {M} for {\\RF} ...\n"); }
    double M[nd*nd];
    for (int32_t r = 0; r < nd; r++)
      { r2_t *hr = &(H[r*ni]);
        for (int32_t s = 0; s <= r; s++)
          { r2_t *hs = &(H[s*ni]);
            double sum = 0.0;
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++)
                  { double hrij = hr[i].c[j];
                    double hsij = hs[i].c[j];
                    double rij = arad[i].c[j];
                    if (rij == 0)
                      { assert((hrij == 0) && (hsij == 0)); }
                    else
                      { sum += hrij*hsij/(rij*rij); }
                  }
              }
            M[r*nd + s] = sum;
            M[s*nd + r] = sum; /* Diag is assigned twice, but OK. */
          }
      }
    
    if (debug) { fprintf(stderr, "... computing the eigen decomp {S,d} of {M} ...\n"); }
    double Q[nd*nd]; /* Eigenvector matrix. */
    double d[nd]; /* Eigenvalues. */
    { /* Convert {M} to tridiag with diagonal {d[0..nd-1]} and subdiagonal {e[1..nd-1]}: */
      double e[nd]; /* Sub-diagonal elements of temporary tridiagonal matrix. */
      syei_tridiagonalize(nd, M, d, e, Q);
      /* Compute eigenvalues and eigenvectors from {Q} and tridiag matrix: */
      int32_t p; /* Number of eigenvalues computed. */
      int32_t absrt = 0; /* Sort eigenvalues by signed value. */
      syei_trid_eigen(nd, d, e, Q, &p, absrt);
      /* Check that all eigenvalues are positive: */
      demand(p == nd, "failed to determine eigenvalues of {M}");
    }
    /* Compute the radii {urad[0..nd-1]} of {\RF}: */
    for (int32_t k = 0; k < nd; k++)
      { demand(d[k] > 0.0, "non-positive eigenvalue");
        urad[k] = 1.0/sqrt(d[k]);
      }
    
    if (debug) { fprintf(stderr, "... computing the basis {U = Q H} aligned with axes of {\\RF} ...\n"); }
    for (int32_t r = 0; r < nd; r++)
      { double *qr = &(Q[r*nd]);
        r2_t *ur = &(U[r*ni]);
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++)
              { double sum = 0.0;
                for (int32_t s = 0; s < nd; s++)
                  { r2_t *hs = &(H[s*ni]);
                    double qrs = qr[s];
                    double hsij = hs[i].c[j];
                    sum += qrs*hsij;
                  }
                ur[i].c[j] = sum;
              }
          }
        if (debug) { r2_align_print_vector(stderr, ni, "u", r, ur, FALSE); }
        if (debug) { fprintf(stderr, "radius {urad[%d]} = %.8f\n", r, urad[r]); }
      }
  }

void r2_align_throw_ball_vector(int32_t ni, double rmin, double rmax, r2_t p[])
  { demand ((0 <= rmin) && (rmin <= 0.80*rmax), "bad norm range");
    while (TRUE) 
      { /* Throw a random alignment vector in the cube of side {2*rmax} centered at the origin, compute its squared norm: */
        double sum2 = 0.0; /* Sum of squares of all coordinates. */
        for (int32_t i = 0; i < ni; i++) 
          { for (int32_t j = 0; j < 2; j++) 
              { double pij = rmax*(2*drandom() - 1.0);
                p[i].c[j] = pij;
                sum2 += pij*pij;
              }
          }
        /* Reject if the norm is outside the requested range:: */
        if ((sum2 <= rmax*rmax) && (sum2 >= rmin*rmin)) { return; }
      }
  }

double r2_align_norm_sqr(int32_t ni, r2_t p[])
  { 
    double sum2 = 0.0; /* Sum of squares of all coordinates. */
    for (int32_t i = 0; i < ni; i++) 
      { for (int32_t j = 0; j < 2; j++) 
          { double pij = p[i].c[j];
            sum2 += pij*pij;
          }
      }
    return sum2;
  }
            
double r2_align_dot(int32_t ni, r2_t p[], r2_t q[])
  { 
    double sdot = 0.0;
    for (int32_t i = 0; i < ni; i++) 
      { for (int32_t j = 0; j < 2; j++) 
          { sdot += p[i].c[j]*q[i].c[j]; }
      }
    return sdot;
  }

double r2_align_dist_sqr(int32_t ni, r2_t p[], r2_t q[])
  { 
    double sum2 = 0.0; /* Sum of squares of all coordinates. */
    for (int32_t i = 0; i < ni; i++) 
      { for (int32_t j = 0; j < 2; j++) 
          { double dij = p[i].c[j] - q[i].c[j];
            sum2 += dij*dij;
          }
      }
    return sum2;
  }
            
double r2_align_rel_disp_sqr(int32_t ni, r2_t p[], r2_t q[], r2_t arad[])
  { double d2 = 0.0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            double pij = p[i].c[j];
            double qij = (q == NULL ? 0.0 : q[i].c[j]);
            if (rij == 0)
              { if (pij != qij) { return +INF; } }
            else
              { double dij = (pij - qij)/rij;
                d2 += dij*dij;
              }
          }
      }
    return d2;
  }
    
bool_t r2_align_coord_is_variable(r2_t arad[], int32_t i, int32_t j)
  { double rij = arad[i].c[j];
    demand(rij >= 0, "invalid search radius");
    return (rij > 0);
  }

void r2_align_points_to_vars(int32_t ni, r2_t p[], r2_t arad[], r2_t ctr[], int32_t nv, double y[])
  { int32_t k = 0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            if (rij != 0)
              { double c0 = (ctr != NULL ? ctr[i].c[j] : 0.0);
                demand(k < nv, "{nv} too small");
                y[k] = p[i].c[j] - c0;
                k++;
              }
          }
      }
    demand(k == nv, "wrong {nv}");
  }

void r2_align_vars_to_points(int32_t nv, double y[], int32_t ni, r2_t arad[], r2_t ctr[], r2_t p[])
  { int32_t k = 0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            if (rij != 0) 
              { double c0 = (ctr != NULL ? ctr[i].c[j] : 0.0);
                demand(k < nv, "{nv} too small");
                p[i].c[j] = c0 + y[k];
                k++;
              }
          }
      }
    demand(k == nv, "wrong {nv}");
  }

void r2_align_throw_ortho_disp_vector(int32_t ni, r2_t arad[], int32_t k, r2_t H[]) 
  { 
    bool_t debug = TRUE;
    if (debug) { fprintf(stderr, "  computing H[%d]\n", k); }

    r2_t *hk = &(H[ni*k]);
    int32_t nit = 0; /* Counts iterations, for safety. */
    while (TRUE)
      { nit++;
        assert(nit < 1000); /* Should never happen if {k < nd} */
      
        /* Generate a vector {hk} uniformly distributed in the unit ball that is not too short: */
        double rmin = 0.5; /* To reduce the effect of roundoff noise on normalization. */
        r2_align_throw_ball_vector(ni, rmin, 1.0, hk);
    
        /* Project {hk} perpendicular to coordinates where {arad} is zero: */
        for (int32_t i = 0; i < ni; i++) 
          { for (int32_t j = 0; j < 2; j++) 
              { double rij = arad[i].c[j];
                demand (rij >= 0, "invalid {arad}");
                if (rij == 0.0) { hk[i].c[j] = 0.0; }
              }
          }
        
        /* Project {hk} perpendicular to the all-ones vectors, preserving conformity: */
        for (int32_t j = 0; j < 2; j++) 
          { double sum = 0.0;
            int32_t nv = 0;
            for (int32_t i = 0; i < ni; i++) 
              { double rij = arad[i].c[j];
                if (rij != 0.0) { sum += hk[i].c[j]; nv++; }
              }
            if (nv > 0)
              { double avg = sum/nv;
                for (int32_t i = 0; i < ni; i++) 
                  { double rij = arad[i].c[j];
                    if (rij != 0.0) { hk[i].c[j] -= avg; }
                  }
              }
          }
        
        /* Project {hk} perpendicular to the previous delta vectors. */
        /* Since the previous vectors are conformal to {arad} and balanced, 
          the projection preserves these properties: */
        for (int32_t r = 0; r < k; r++)
          { r2_t *hr = &(H[ni*r]);
            double sdot = r2_align_dot(ni, hk, hr);
            for (int32_t i = 0; i < ni; i++) 
              { for (int32_t j = 0; j < 2; j++) 
                  { hk[i].c[j] -= sdot*hr[i].c[j]; }
              }
          }
          
        /* Check if norm is still large enough: */
        double sum2 = r2_align_norm_sqr(ni, hk);
        if (sum2 >= rmin*rmin)
          { /* Normalize and return this vector: */
            double norm = sqrt(sum2);
            for (int32_t i = 0; i < ni; i++) 
              { for (int32_t j = 0; j < 2; j++) 
                  { hk[i].c[j] /= norm; }
              }
            return;
          }
          
        /* Remaining vector {hk} was too short. Try again. */
      }
  }

void r2_align_print_vector(FILE *wr, int32_t ni, char *name, int32_t ix, r2_t p[], bool_t ctvars)
  { 
    int32_t nz = 0;
    for (int32_t i = 0; i < ni; i++)
      { fprintf(wr, "  %s", name);
        if (ix >= 0) { fprintf(wr, "[%d]", ix); }
        fprintf(wr, "[%d] = (", i);
        for (int32_t j = 0; j < 2; j++)
          { double pij = p[i].c[j];
            fprintf(wr, " %12.8f", pij); 
            if (pij == 0) { nz++; }
          }
        fprintf(wr, " )\n");
      }
    if (ctvars)
      { int32_t nc = 2*ni;
        i2_t nv = r2_align_count_variable_coords(ni, p);
        int32_t nd = r2_align_count_degrees_of_freedom (ni, p);
        fprintf(wr, "  nc = %d nz = %d nv = %d+%d = %d nd = %d\n", nc, nz, nv.c[0], nv.c[1], nv.c[0]+nv.c[1], nd);
      }
  }

void r2_align_plot_mismatch
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    double step,
    r2_align_mismatch_t *F2
  )
  {
    bool_t debug = TRUE;

    /* Choose two conformal balanced orthonormal delta vectors: */
    int32_t nd = r2_align_count_degrees_of_freedom(ni, arad);
    demand(nd >= 2, "not enough degrees of freedom to plot");
    r2_t U[nd*ni];
    double urad[nd];
    r2_align_compute_search_ellipsoid (ni, arad, nd, U, urad);
    
    /* Choose the two plot axis vectors: */
    auto void choose_plot_axis (int32_t k, r2_t **ukP, int32_t *nkP);
      /* Chooses the unit displacement {*ukP} and the number of steps {*nkP}
        for plot axis {k} in {0..1}. */
    
    r2_t *u0 = NULL; int32_t n0 = -1; choose_plot_axis (0, &u0, &n0);
    r2_t *u1 = NULL; int32_t n1 = -1; choose_plot_axis (1, &u1, &n1);

    r2_align_plot_mismatch_lines(wr, ni, ctr, arad, step, u0, n0, u1, n1, F2);

    fflush(wr);
    return;
    
    void choose_plot_axis (int32_t k, r2_t **ukP, int32_t *nkP)
      { r2_t *uk = &(U[k*ni]); 
        double uradk = urad[k];
        int32_t nk = (int32_t)floor(uradk/step);
        if (debug)
          { fprintf(stderr, "plot axis %d length = %.8f steps = %d\n", k, uradk, nk);
            r2_align_print_vector(stderr, ni, "u", k, uk, FALSE);
          }
        (*ukP) = uk;
        (*nkP) = nk;
      }
   }
    
void r2_align_plot_mismatch_lines
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    double step,  /* Grid step. */
    r2_t u0[],    /* Direction across lines. */
    int32_t n0,   /* Number of lines on each side of center. */
    r2_t u1[],    /* Direction along lines. */
    int32_t n1,   /* Max number of points along each line on each side of center. */
    r2_align_mismatch_t *F2
  )
  { 
    /* Plot on a grid along {u0,u1}: */
    r2_t psmp[ni]; /* Probe alignment vector. */
    for (int32_t s0 = -n0; s0 <= +n0; s0++)
      { double d0 = s0*step;
        for (int32_t s1 = -n1; s1 <= +n1; s1++)
          { double d1 = s1*step;
            /* Compute the probe vector {psmp[0..ni-1]}. */
            for (int32_t i = 0; i < ni; i++)
              { r2_mix(d0, &(u0[i]), d1, &(u1[i]), &(psmp[i])); 
                r2_add(&(ctr[i]), &(psmp[i]), &(psmp[i])); 
              }
            /* Evaluate the function and plot: */
            double rdist2 = r2_align_rel_disp_sqr(ni, ctr, psmp, arad);
            if (rdist2 <= 1.0) 
              { double F2val = F2(ni, psmp);
                fprintf(wr, "%+9.6f %+9.6f  %12.6f\n", d0, d1, F2val);
              }
          }
        /* Blank line between scanlines, for {gnuplot}: */
        fprintf(wr, "\n"); 
      }
   }

void r2_align_throw_arad(int32_t ni, r2_t zfrac, double rmin, double rmax, r2_t arad[])
  { 
    for (int32_t j = 0; j < 2; j++)
      { /* Decide how many coordinates {arad[0..ni-1].c[j]} will be set to zero. */
        int32_t nzgoal = (int32_t)floor(zfrac.c[j]*ni + 0.5); 
        if (ni >= 2)
          { /* Cook {nzgoal} to be {0} or {ni} iff {zfrac} is 0.0 or 1.0, resp.: */
            if ((nzgoal ==  0) && (zfrac.c[j] > 0.0)) { nzgoal = 1; }
            if ((nzgoal == ni) && (zfrac.c[j] < 1.0)) { nzgoal = ni-1; }
          }
        assert(nzgoal <= ni);
        int32_t nc = ni;      /* Number of coordinates remaining. */
        int32_t nz = nzgoal;  /* Number of coordinates still to be set to zero. */
        for (int32_t i = 0; i < ni; i++)
          { assert(nz <= nc);
            double rij;
            if ((nz == nc) || (nc*drandom() < nz))
              { rij = 0.0; nz--; }
            else
              { rij = rmin + (rmax - rmin)*drandom(); }
            arad[i].c[j] = rij;
            nc--;
          }
      }
  }

