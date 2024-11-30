/* See {neuromat_eeg_func_basis.h}. */
/* Last edited on 2023-02-25 15:24:12 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <rn.h>
#include <rmxn.h>
#include <affirm.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_geom.h>
#include <neuromat_eeg_func_basis.h>

void neuromat_eeg_func_basis_eval
  ( int32_t ne, 
    double bval[],
    r3_t *p3D,
    neuromat_eeg_func_basis_eval_t mother,
    double L[],
    bool_t norm
  )
  {
    mother(ne, bval, p3D);
    if (L != NULL)
      { /* Apply specified linear map: */
        double ival[ne];
        rmxn_map_col(ne, ne, L, bval, ival); 
        rn_copy(ne, ival, bval);
      }
    if (norm)
      { double sum = rn_sum(ne, bval);
        if (fabs(sum) < 1.0e-8)
          { /* Set all basis values to zero: */
            for (uint32_t ie = 0;  ie < ne; ie++) { bval[ie] = 0.0;  }
          }
        else
          { rn_scale(ne, 1/sum, bval, bval); }
      }
    return;
  }

double *neuromat_eeg_func_basis_nearest_dists(int32_t ne, r3_t pos3D[])
  {
    double *erad = rn_alloc(ne);
    for (uint32_t ie = 0;  ie < ne; ie++)
      { /* Find the distance squared {d2min} from {pos[ie]} to its nearest neighbor: */
        r3_t *pi = &(pos3D[ie]);
        double d2min = +INF;
        for (uint32_t je = 0;  je < ne; je++)
          { if (ie != je) 
              { r3_t *pj = &(pos3D[je]);
                double d2j = r3_dist_sqr(pi, pj);
                demand(d2j > 1.0e-5, "two electrodes are nearly coincident");
                if (d2j < d2min) { d2min = d2j; }
              }
          }
        erad[ie] = sqrt(d2min);
      }
    return erad;
  }
  
double neuromat_eeg_func_basis_shepard_weight(r3_t *p, r3_t *ctr, double rho, double sigma, double ord)
  {
    double d2 = r3_dist_sqr(p, ctr);
    double t2shep = (d2 + 1.0e-38)/(rho*rho);
    double zshep = (ord == 2 ? t2shep : pow(t2shep, ord/2));
    double shep = 1/zshep;
    double z2bell = d2/(sigma*sigma);
    if (z2bell >= 90.0) { return 0.0; }
    double bell = exp(-z2bell/2);
    return shep*bell;
  }
               
double neuromat_eeg_func_basis_gauss_bell(r3_t *p, r3_t *ctr, double sigma)
  { double d2 = r3_dist_sqr(p, ctr);
    if (d2 == 0) { return 1.0; }
    double z2bell = d2/(sigma*sigma);
    if (z2bell >= 90) { return 0.0; }
    double bell = exp(-z2bell/2);
    return bell;
  }

double neuromat_eeg_func_basis_mexican_hat(r3_t *p, r3_t *ctr, double sigma, double tau)
  { double d2 = r3_dist_sqr(p, ctr);
    if (d2 == 0) { return 1.0; }
    double z2bell = d2/(sigma*sigma);
    if (z2bell >= 90) { return 0.0; }
    double bell = exp(-z2bell/2);
    double zsinc = sqrt(d2)/tau;
    double sinc = sin(M_PI*zsinc)/(M_PI*zsinc*(1 + zsinc));
    return bell*sinc;
  }

double neuromat_eeg_func_basis_voronoi_ind(r3_t *p, int32_t ie, int32_t ne, r3_t pos3D[])
  { /* Find the nearest electrode to {p3D}: */
    double d2min = +INF;
    int32_t jemin = -1;
    for (uint32_t je = 0;  je < ne; je++) 
       { double d2 = r3_dist_sqr(p, &(pos3D[je]));
         if (d2 < d2min) { d2min = d2; jemin = je; }
       }
    assert(jemin >= 0);
    return (double)(jemin == ie);
  }

double *neuromat_eeg_func_basis_lagrangian_matrix(int32_t ne, neuromat_eeg_func_basis_eval_t mother, r3_t pos3D[])
  {
    bool_t verbose = FALSE;

    /* Build colocation matrix {A} such that {A[ie,je]} is the value of element {ie} on point {pos3d[je]}: */
    double bval[ne];
    double *A = rmxn_alloc(ne, ne);
    for (uint32_t je = 0;  je < ne; je++)
      { r3_t *pj = &(pos3D[je]);
        mother(ne, bval, pj);
        /* Fill column {je} of matrix {A}:  */
        for (uint32_t ie = 0;  ie < ne; ie++)
          { double Aij = bval[ie];
            if (verbose && (Aij != 0)) { fprintf(stderr, "  A[%3d,%3d] = %+12.7f\n", ie, je, Aij); }
            A[ie*ne + je] = Aij;
          }
      }
               
    /* Invert matrix {A} to get the matrix {L} that maps from raw elem values to interp elem values: */
    assert(ne == ne);
    double *L = rmxn_alloc(ne, ne);
    (void)rmxn_inv(ne, A, L);
    /* Check inversion: */
    { double *R = rmxn_alloc(ne, ne);
      rmxn_mul(ne, ne, ne, L, A, R);
      int32_t nerr = 0;
      for (uint32_t ie = 0;  ie < ne; ie++)
        { for (uint32_t je = 0;  je < ne; je++)
            { double Rij = R[ie*ne + je];
              double Iij = (ie == je ? 1.0 : 0.0);
              double err = fabs(Rij - Iij);
              if (err > 1.0e-7) 
                { fprintf(stderr, "  inversion error R[%3d,%3d] = %+12.8f\n", ie, je, Rij);
                  nerr++;
                }
            }
        }
      demand(nerr == 0, "aborted");
      free(R);
    }
    free(A);
    return L;
  }

    
