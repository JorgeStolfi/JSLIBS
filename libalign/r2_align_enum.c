/* See {r2_align_enum.h}. */
/* Last edited on 2023-09-07 18:11:05 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <jsmath.h>
#include <affirm.h>
#include <r2_align.h>

#include <r2_align_enum.h>

/* INTERNAL PROTOTYPES */

void r2_align_enum_grid
  ( int32_t ni,               /* Number of images to align. */
    r2_align_mismatch_t *F2,  /* Function that evaluates the mismatch between the images. */
    r2_t ctr[],               /* (IN/OUT) Corresponding points in each image. */
    r2_t arad[],              /* Max delta vector coordinates for each image. */
    int32_t nd,               /* Degrees of freedom (dimension of {\RF}). */
    r2_t U[],                 /* Adjstment vectors parallet to axes of {\RF}. */
    double urad[],            /* Radii of {\RF} along those directions. */
    double step,              /* Grid step. */
    r2_t p[],                 /* (OUT) Optimum alignment vector. */
    double *F2val_P           /* (OUT) Value of {F2} for the optimum alignment vector. */
  );
  /* Enumerates every point {psmp} that is {ctr} plus a linear combination of the
    lattice vectors {U[k]} for {k} in {0..nv-1} with coefficients that are {step} times
    integers in the range {-tmax..+tmax}; where {U[k][i] = U[k*ni + i]} for 
    {i} in {0..ni-1}.   
    
    The {arad} parameter is used only to skip coordinates where {\RE} has size zero. 
    
    Then evaluates {F2(ni, psmp)} at every alignment vector {psmp} 
    that is inside the ellipsoid with center {ctr} and radius vector {arad}.
    Sets {*F2val_P} to the minimum of those values, and sets {p[0..ni-1]}
    to the alignment vector {psmp} that realizes the minimum. */

/* IMPLEMENTATIONS */

void r2_align_enum
  ( int32_t ni,               /* Number of images to align. */
    r2_align_mismatch_t *F2,  /* Function that evaluates the mismatch between the images. */
    r2_t arad[],              /* Max delta vector coordinates for each image. */
    bool_t bal,              /* True if alignment vector adjustments should be balanced. */
    double tol,               /* Desired precision. */
    r2_t p[],                 /* (IN/OUT) Initial and optimized alignment vector. */
    double *F2val_P           /* (OUT) Mismatch for the computed alignment vector. */
  )
  {
    demand(tol > 0, "invalid {tol}");
    i2_t nv = r2_align_count_variable_coords (ni, arad);
    int32_t nd = r2_align_count_degrees_of_freedom(nv, bal);

    /* Save the initial guess: */
    r2_t ctr[ni]; /* Center of ellipsoid (aved initial guess). */
    for (int32_t i = 0; i < ni; i++) { ctr[i] = p[i]; }

    /* Compute the axes and radii of the search ellipsoid {\RF}: */
    r2_t U[nd*ni]; /* Basis of conformal balanced delta vectors. */
    double urad[nd];
    r2_align_compute_search_ellipsoid (ni, arad, bal, nd, U, urad);

    /* !!! Convert {U} to the packed lattice basis !!! */

    r2_align_enum_grid(ni, F2, ctr, arad, nd, U, urad, tol, p, F2val_P);
  }
  
void r2_align_enum_grid
  ( int32_t ni,               /* Number of images to align. */
    r2_align_mismatch_t *F2,  /* Function that evaluates the mismatch between the images. */
    r2_t ctr[],               /* (IN/OUT) Corresponding points in each image. */
    r2_t arad[],              /* Max delta vector coordinates for each image. */
    int32_t nd,               /* Degrees of freedom (dimension of {\RF}). */
    r2_t U[],                 /* Adjstment vectors parallet to axes of {\RF}. */
    double urad[],            /* Radii of {\RF} along those directions. */
    double step,              /* Grid step. */
    r2_t p[],                 /* (OUT) Optimum alignment vector. */
    double *F2val_P            /* (OUT) Value of {F2} for the optimum alignment vector. */
  )    
  {
    bool_t debug = FALSE;
    
    /* Compute the number of samples along each main axis of {\RF}: */
    int32_t n[nd];     /* Number of sample points on each side of 0 along each axis of {\RF}. */
    for (int32_t k = 0; k < nd; k++) 
      { n[k] = (int32_t)floor(urad[k]/step);
        if (debug)
          { fprintf(stderr, "plot axis %d length = %.8f steps = %d\n", k, urad[k], n[k]);
            r2_t *uk = &(U[k*ni]);
            r2_align_print_vector(stderr, ni, "u", k, uk);
          }
      }
    
    /* Enumerate all integer tuples {t[0..nd-1]} where {t[k]} ranges in {-n[k]..+n[k]}: */
    (*F2val_P) = +INF;  /* Minimum mismatch found so far. */
    int32_t t[nd];     /* Enumeration variables. */
    for (int32_t k = 0; k < nd; k++) { t[k] = -n[k]; }
    r2_t psmp[ni];     /* Sampling point. */
    int32_t knext = -1; /* Next tuple elem to be incremented is {t[knext]}. */
    while (knext < nd)
      { if (debug)
          { fprintf(stderr, "  t = ( ");
            for (int32_t k = 0; k < nd; k++) { fprintf(stderr, " %d", t[k]); }
            fprintf(stderr, " )\n");
          }
        
        /* !!! Should avoid generating tuples outside {\RF} !!! */
        
        /* Compute the {urad}-relative squared norm {sum2} of the displacement: */
        double sum2 = 0; /* Squared norm of {(p-ctr)/arad}. */
        for (int32_t k = 0; k < nd; k++)  { double ek = t[k]*step/urad[k]; sum2 += ek*ek; }
        
        if (sum2 <= 1.0 + 1.0e-8)
          { /* Sample point will be inside {\RF}. */
        
            /* Build the sample point {psmp}: */
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++)
                  { psmp[i].c[j] = ctr[i].c[j];
                    double rij = arad[i].c[j];
                    if (rij != 0)
                      { double dij = 0;
                        for (int32_t k = 0; k < nd; k++) 
                          { r2_t *uk = &(U[k*ni]);
                            dij += t[k]*step*uk[i].c[j];
                          }
                        psmp[i].c[j] += dij;
                      }
                  }
              }

            /* Evaluate the mismatch function at {psmp}: */
            double F2val = F2(ni, psmp);
            
            if (F2val < (*F2val_P))
              { /* Update the current optimum: */
                for (int32_t i = 0; i < ni; i++) { p[i] = psmp[i]; }
                (*F2val_P) = F2val;
              }
          }

        /* Get the next tuple {t[0..nd-1]}: */
        if ((knext < 0) || (t[knext] >= n[knext]))
          { /* Can't increment this, try next one: */
            knext++;
          }
        else
          { /* Increment {t[knext]} and clear all previous {t[k]}: */
            t[knext]++;
            while (knext > 0) { knext--; t[knext] = -n[knext]; }
          }
      }
  }

