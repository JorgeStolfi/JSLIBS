/* See {r2_align.h}. */
/* Last edited on 2024-11-08 00:03:29 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <r2.h>
#include <i2.h>
#include <rmxn.h>
#include <rmxn_ellipsoid.h>
#include <sym_eigen.h>

#include <r2_align.h>

/* INTERNAL PROTOTYPES */

bool_t r2_align_coord_is_variable (r2_t arad[], int32_t i, int32_t j);
  /* Returns {TRUE} iff {arad[i].c[j]} is nonzero. */

void r2_align_plot_mismatch_lines
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    double step,  /* Grid step. */
    r2_t u0[],    /* Direction across lines. */
    double urad0, /* Max extent in direction {u0}. */
    r2_t u1[],    /* Direction along lines. */
    double urad1, /* Max extent in direction {u1}. */
    r2_align_mismatch_t *F2
  );
  /* Similar to {r2_align_plot_mismatch}, but samples the function only
    for {n0} straight paths spaced by {step} along direction {u0}. Each path will
    have {n1} points spaced by {step} in direction {u1}. For each
    sampled point {p}, writes to the file a line with two independent coordinates
    {s0,s1} along {u0,u1}, and the value of {F2(ni,p)} . */

/* ---------------------------------------------------------------------- */

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
  { int32_t nd = nv.c[0] + nv.c[1];
    if (bal) 
      { if (nv.c[0] >= 2) { nd--; }
        if (nv.c[1] >= 2) { nd--; }
      }
    return nd;
  }

void r2_align_compute_search_ellipsoid
  ( int32_t ni,
    r2_t arad[],
    bool_t bal, 
    int32_t nd,
    r2_t U[], 
    double urad[]
  )
  { 
    bool_t debug = FALSE;
    bool_t verbose_norm = TRUE;
    
    if (debug) { fprintf(stderr, "--> entering %s\n", __FUNCTION__); }

    demand(ni >= 0, "invalid {ni}");
    i2_t nv = r2_align_count_variable_coords(ni, arad); /* Variable coords count per axis of {\R^2}. */
    int32_t nh = nv.c[0] + nv.c[1]; /* Total count of variable coords of {\RC}. */
    int32_t nd2 = r2_align_count_degrees_of_freedom(nv, bal);
    demand(nd2 == nd, "incorrect search space dimension {nd}");

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
        
        if (debug) { fprintf(stderr, "  ... orthonormalizing the balancing and conformation constraints\n"); } 
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
          { r2_t *Uk = &(U[kd*ni]);
            double *Vk = &(V[kd*2*ni]);
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
      { /* Simple optimization, each variable coord is an indep minimization variable: */
        int32_t kd = 0; /* Count nonzero elements of {arad}. */
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
        assert(kd == nd);
      }
  }

bool_t r2_align_coord_is_variable(r2_t arad[], int32_t i, int32_t j)
  { double rij = arad[i].c[j];
    demand(rij >= 0, "invalid search radius");
    return (rij > 0);
  }

void r2_align_print_vector(FILE *wr, int32_t ni, char *name, int32_t ix, r2_t p[])
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
  }

void r2_align_plot_mismatch
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    bool_t bal,
    double step,
    r2_align_mismatch_t *F2
  )
  {
    /* Compute degrees of freedom: */
    i2_t nv = r2_align_count_variable_coords (ni, arad);
    fprintf(stderr, "num of variable coords nv = (%d,%d)\n", nv.c[0], nv.c[1]);
    int32_t nd = r2_align_count_degrees_of_freedom(nv, bal);
    fprintf(stderr, "search dimensions nd = %d\n", nd);
    demand(nd >= 2, "not enough degrees of freedom to plot");
    
    /* Compute the search ellipsoid: */
    r2_t U[nd*ni];
    double urad[nd];
    r2_align_compute_search_ellipsoid(ni, arad, bal, nd, U, urad);
    
    /* Choose the two plot axis vectors: */
    r2_t *u0 = &(U[0*ni]); double urad0 = urad[0];
    r2_t *u1 = &(U[1*ni]); double urad1 = urad[1];

    r2_align_plot_mismatch_lines(wr, ni, ctr, arad, step, u0, urad0, u1, urad1, F2);

    fflush(wr);
    return;
  }
    
void r2_align_plot_mismatch_lines
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    double step,  /* Grid step. */
    r2_t u0[],    /* Direction across lines. */
    double urad0, /* Max extent in direction {u0}. */
    r2_t u1[],    /* Direction along lines. */
    double urad1, /* Max extent in direction {u1}. */
    r2_align_mismatch_t *F2
  )
  { 
    bool_t debug = FALSE;
    auto int32_t num_points(int32_t k, r2_t *uk, double uradk);
      /* Chooses the unit displacement {*ukP} and the number of steps {*nkP}
        for plot axis {k} in {0..1}. */
    
    /* Get the number of points on each axis: */
    int32_t n0 = num_points(0, u0, urad0);
    int32_t n1 = num_points(1, u1, urad1);

    /* Plot on a grid along {u0,u1}: */
    r2_t psmp[ni]; /* Probe alignment vector. */
    for (int32_t s0 = -n0; s0 <= +n0; s0++)
      { double d0 = s0*step;
        for (int32_t s1 = -n1; s1 <= +n1; s1++)
          { double d1 = s1*step;
            /* Check if inside the ellipsoid: */
            double du0 = d0/urad0, du1 = d1/urad1;
            if (du0*du0 + du1*du1 <= 1.0)
              { /* Compute the probe vector {psmp[0..ni-1]}. */
                for (int32_t i = 0; i < ni; i++)
                  { r2_mix(d0, &(u0[i]), d1, &(u1[i]), &(psmp[i])); 
                    r2_add(&(ctr[i]), &(psmp[i]), &(psmp[i])); 
                  }
                /* Evaluate the function and plot: */
                double F2val = F2(ni, psmp);
                fprintf(wr, "%+9.6f %+9.6f  %12.6f\n", d0, d1, F2val);
              }
          }
        /* Blank line between scanlines, for {gnuplot}: */
        fprintf(wr, "\n"); 
      }
    return;
    
    /* INTERNAL FUNC IPLEMENTATIONS */
    
    int32_t num_points(int32_t k, r2_t *uk, double uradk)
      { int32_t nk = (int32_t)floor(uradk/step);
        if (debug)
          { fprintf(stderr, "plot axis %d length = %.8f steps = %d\n", k, uradk, nk);
            r2_align_print_vector(stderr, ni, "u", k, uk);
          }
        return nk;
      }
   }

void r2_align_throw_vector (int32_t ni, double rmin, double rmax, r2_t p[])
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
        /* Return {p} if its norm is within the requested range:: */
        if ((sum2 <= rmax*rmax) && (sum2 >= rmin*rmin)) { return; }
      }
  }

void r2_align_throw_ctr (int32_t ni, double cmax, r2_t ctr[], bool_t verbose)
  { 
    r2_align_throw_vector(ni, 0.0, cmax, ctr);
    if (verbose) 
      { fprintf(stderr, "center alignment {ctr}:\n");
        r2_align_print_vector(stderr, ni, "ctr", -1, ctr);
      }
    return;
  }
  
void r2_align_throw_arad (int32_t ni, double rmax, int32_t nvmin, r2_t arad[], bool_t verbose)
  { 
    double rmin = 0.5*rmax;
    /* Choose the ideal fraction {zfrac.c[j]} of non-zeros in each coordinate {j}: */
    r2_t vfrac = (ni <= 2 ? (r2_t){{ 1.00, 1.00 }} : (r2_t){{ 0.25, 0.75 }} );
    for (int32_t j = 0; j < 2; j++)
      { /* Decide how many coordinates {arad[0..ni-1].c[j]} will be non-zero. */
        int32_t nvgoal = (int32_t)floor(vfrac.c[j]*ni + 0.5); 
        if ((nvgoal == ni) && (vfrac.c[j] < 1.0) && (ni >= 2)) { nvgoal = ni-1; }
        if (nvgoal < nvmin) { nvgoal = nvmin; }
        if (nvgoal > ni) { nvgoal = ni; }
        int32_t nc = ni;      /* Number of coordinates remaining. */
        int32_t nv = nvgoal;  /* Number of coordinates still to be set to non-zero. */
        for (int32_t i = 0; i < ni; i++)
          { assert(nv <= nc);
            double rij;
            if ((nv == nc) || (nc*drandom() < nv))
              { rij = rmin + (rmax - rmin)*drandom(); nv--; }
            else
              { rij = 0.0; }
            arad[i].c[j] = rij;
            nc--;
          }
      }
    if (verbose) 
      { fprintf(stderr, "adjustment search radius {arad}:\n");
        r2_align_print_vector(stderr, ni, "arad", -1, arad);
      }
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
    
