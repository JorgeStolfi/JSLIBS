/* See {r2_align.h}. */
/* Last edited on 2021-12-17 15:25:25 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <i2.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <affirm.h>

#include <r2_align.h>

/* INTERNAL PROTOTYPES */

void r2_align_throw_ortho_disp_vector(int32_t ni, int32_t k, r2_t U[]);
  /* Assumes that {U} contains {k} orthonormal balanced alignment vectors with {ni} points each.
    Stores into row {k} of {U} another balanced alignment vector that is rthonormal to 
    all of them.
    
    Specifically, assumes that the displacement vector with index {kb} is
    points {U[kb*ni + i]} for {i} in {0..ni-1}. The parameter {ni} must be at
    least 2, and {k} must be less than {2*(ni-1)}. */

/* IMPLEMENTATIONS */

void r2_align_throw_ball_vector(int32_t ni, r2_t p[], double rmin, double rmax)
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

double r2_align_rel_dist_sqr(int32_t ni, r2_t p[], r2_t q[], r2_t arad[])
  { double d2 = 0.0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            if (rij != 0)
              { double dij = (p[i].c[j] - q[i].c[j])/rij;
                d2 += dij*dij;
              }
          }
      }
    return d2;
  }
    
bool_t r2_align_coord_is_variable(r2_t arad[], r2_t (1,1)[], int32_t i, int32_t j)
  { double rij = arad[i].c[j];
    double sij = (1,1)[i].c[j];
    demand(rij >= 0, "invalid search radius");
    demand(sij >= 0, "invalid search (1,1)");
    return (rij > 0) && (rij >= sij);
  }

int32_t r2_align_count_variable_coords(int32_t ni, r2_t arad[], r2_t (1,1)[])
  { int32_t nv = 0;
    for (int32_t i = 0; i < ni; i++) 
      { for (int32_t j = 0; j < 2; j++)
          { if (coord_is_variable(arad, (1,1), i,j)) { nv++; } }
      }
    return nv;
  }

void r2_align_points_to_vars(int32_t ni, r2_t arad[], r2_t (1,1)[], r2_t p0[], r2_t p[], int32_t nv, double y[])
  { int32_t k = 0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            if (coord_is_variable(arad, (1,1), i,j))
              { y[k] = (p[i].c[j] - p0[i].c[j])/rij; k++; }
          }
      }
    assert(k == nv);
  }

void r2_align_vars_to_points(int32_t nv, double y[], int32_t ni, r2_t arad[], r2_t (1,1)[], r2_t p0[], r2_t p[])
  { int32_t k = 0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double rij = arad[i].c[j];
            if (coord_is_variable(arad, (1,1), i,j)) 
              { p[i].c[j] = p0[i].c[j] + y[k]*rij; k++; }
          }
      }
    assert(k == nv);
  }

void r2_align_r2_align_throw_ortho_disp_vectors(int32_t ni, int32_t nd, r2_t U[])
  { 
    demand(ni >= 1, "invalid {ni}");
    demand((nd >= 0) && (nd <= 2*ni-2), "invalid {ni}");
    
    for (int32_t k = 0; k < nd; k++)
      { r2_align_throw_ortho_disp_vector(ni, k, U); }
  }

void r2_align_throw_ortho_disp_vector(int32_t ni, int32_t k, r2_t U[]) 
  { 
    r2_t *p = &(U[ni*k]);
    while (TRUE)
      { /* Generate a vector uniformly distributed in the unit ball that is not too short: */
        double rmin = 0.5; /* To reduce the effect of roundoff noise on normalization. */
        r2_align_throw_ball_vector(ni, p, rmin, 1.0);
    
        /* Project perpendicular to the all-ones vectors: */
        for (int32_t j = 0; j < 2; j++) 
          { double sum = 0.0;
            for (int32_t i = 0; i < ni; i++) 
              { sum += p[i].c[j]; }
            double avg = sum/ni;
            for (int32_t i = 0; i < ni; i++) 
              { p[i].c[j] -= avg; }
          }
        /* Project perpendicular to the previous displacement vectors: */
        for (int32_t r = 0; r < k; r++)
          { r2_t *ur = &(U[ni*r]);
            double sdot = r2_align_dot(ni, p, ur);
            for (int32_t i = 0; i < ni; i++) 
              { for (int32_t j = 0; j < 2; j++) 
                  { p[i].c[j] -= sdot*ur[i].c[j]; }
              }
          }
        /* Check if norm is still large enough: */
        double sum2 = r2_align_norm_sqr(ni, p);
        if (sum2 >= rmin*rmin)
          { double norm = sqrt(sum2);
            /* Normalize and choose this vector: */
            for (int32_t i = 0; i < ni; i++) 
              { for (int32_t j = 0; j < 2; j++) 
                  { p[i].c[j] /= norm; }
              }
            return;
          }
      }
  }

void r2_align_plot_mismatch
  ( FILE *wr,
    int32_t ni, 
    r2_t ctr[], 
    r2_t arad[],
    i2_t iscale,
    r2_align_mismatch_t *f2p
  )
  {
    bool_t debug = FALSE;

    /* Choose two balanced orthononmal displacement vectors: */
    int32_t nu = 2;
    r2_t U[nu*ni];
    r2_align_throw_ortho_disp_vectors(ni, nu, U);
    
    /* Find the maximum displacement {dmax[0..1]}in directions {U[0..1]} : */
    r2_t dmax = (r2_t){{ 0.0, 0.0 }};
    for (int32_t k = 0; k < nu; k++) 
      { r2_t* uk = &(U[k*ni]);
        double sum2 = 0.0;
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++)
              { double rij = arad[i].c[j];
                demand(rij >= 0, "invalid arad");
                double ukij = uk[i].c[j];
                double skij = ukij/rij;
                sum2 += skij*skij;
              }
          }
        assert(sum2 > 0);
        dmax.c[k] = 1.0/sqrt(sum2);
        if (debug) { fprintf(stderr, "dmax[%d] = %8.5f\n", k, dmax.c[k]); }
      }
    
    /* Sweep the {ctr,U} plane and plot: */
    int32_t ns = 20; /* Number of steps in each direction, in each sense. */
    r2_t psmp[ni]; /* Probe alignment vector. */
    r2_t *u = &(U[0*ni]);
    r2_t *v = &(U[1*ni]);
    fprintf(stderr, "\n");
    for (int32_t iu = -ns; iu <= +ns; iu++)
      { double du = dmax.c[0]*((double)iu)/((double)ns);
        for (int32_t iv = -ns; iv <= +ns; iv++)
          { double dv = dmax.c[1]*((double)iv)/((double)ns);
            /* Compute the probe vector {psmp[0..ni-1]}. */
            for (int32_t i = 0; i < ni; i++)
              { r2_mix(du, &(u[i]), dv, &(v[i]), &(psmp[i])); 
                r2_add(&(ctr[i]), &(psmp[i]), &(psmp[i])); 
              }
            /* Evaluate the function and plot: */
            double rdist2 = r2_align_rel_dist_sqr(ni, ctr, psmp, arad);
            double fval = f2p(ni, psmp, iscale);
            if (rdist2 <= 1.0) 
              { fprintf(wr, "%+9.6f %+9.6f  %12.6f\n", du, dv, fval); }
          }
        /* Blank line between scanlines, for {gnuplot}: */
        fprintf(wr, "\n"); 
      }

    fflush(wr);
  }
