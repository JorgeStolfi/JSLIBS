/* See {hr2_pmap_opt.h}. */
/* Last edited on 2023-10-08 11:38:21 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
 
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <r3x3.h>
#include <i2.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <affirm.h>
#include <sve_minn.h>

#include <hr2_pmap_opt.h>

void hr2_pmap_opt_quad
  ( hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    r3x3_t *R,                   /* Max adjustment for {A.dir}. */
    double rtol,                 /* Desired relative adjustment precision for {*A}. */
    hr2_pmap_t *A,               /* (IN/OUT) The affine map to adjust. */
    double *f2AP                 /* (OUT) Goal function value for the output {*A}. */
  )
  {
    int32_t maxIters = 10;
    bool_t debug = TRUE;
    
    (*f2AP) = f2(A); /* Current goal function on {A}. */

    /* Identify the elements to optimize: */
    int32_t nv;
    double *Ap[8]; /* Pointers to adjustable elements of {A.dir}. */
    double Re[8];  /* Max adjustment amount of each element. */
    hr2_pmap_opt_get_var_elems(&(A->dir), R, Ap, Re, &nv);
    fprintf(stderr, "%d adjustable elements found\n", nv);
    assert (nv <= 8);
   
    if (nv > 0)
      { 
        /* Apply nonlinear optimization. */

        /* Save initial elem values (center of search region): */ 
        double Ce[nv]; 
        for (int32_t iv = 0; iv < nv; iv++) { Ce[iv] = (*(Ap[iv])); }

        /* The optimization variable {x[iv]}, for {iv} in {0..nv-1}, is the 
          difference {((*Ap)[iv]-Ce[iv]} divided by the corresponding radius
          {Re[iv]}. */
          
        auto void extract_vars(double y[]);
          /* Converts the mutable elements of {*A} into 
            the optmizaton variables {y[0..nv-1]}. */

        auto void store_vars(double y[]);
          /* Converts the optimization variables {y[0..nv-1]} into 
            the mutable elements of {*A}. */

        auto double sve_goal(int32_t nx, double x[]);
          /* Computes the minimization goal function from the given argument {x[0..nx-1]}.
            Expects {nx == nv}. Also sets {*A} from {x[0..nx-1]}. */

        /* Convert the projective map {*A} to the optimization variables {z[0..nv-1]}: */
        double z[nv];
        extract_vars(z);
        /* Compute the initial goal function value: */
        double Fz = sve_goal(nv, z);
        sign_t dir = -1; /* Look for minimum. */
        sve_minn_iterate
          ( nv, 
            &sve_goal, NULL, 
            z, &Fz,
            dir,
            /*dMax:*/ 1.0,
            /*dBox:*/ TRUE,
            /*rIni:*/ 0.5, 
            /*rMin:*/ 0.05, 
            /*rMax:*/ 0.5, 
            /*stop:*/ 0.01,
            maxIters,
            debug
          );

        /* Convert the optimization variables {z[0..nv-1]} to projective map {*A}: */
        store_vars(z);

        /* Local implementations: */

        void extract_vars(double y[])
          { for (int32_t iv = 0; iv < nv; iv++)
              { assert (Re[iv] != 0);
                y[iv] = ((*(Ap[iv])) - Ce[iv])/Re[iv];
              }
          }

        void store_vars(double y[])
          { for (int32_t iv = 0; iv < nv; iv++)
              { (*(Ap[iv])) = Ce[iv] + y[iv]*Re[iv]; }
            r3x3_inv(&(A->dir), &(A->inv));
          }

        double sve_goal(int32_t nx, double x[])
          { assert(nx == nv);
            /* Convert variables {x[0..nx-1]} to displacements {*A}: */
            store_vars(x);
            /* Evaluate the client function: */
            double Q2 = f2(A);
            (*f2AP) = Q2;
            return Q2;
          }
      }

    /* Compue the final mismatch: */
    (*f2AP) = f2(A);
    return;
  }
    
void hr2_pmap_opt_get_var_elems
  ( r3x3_t *M, 
    r3x3_t *R, 
    double *Mp[],
    double Re[],
    int32_t *nvP
  )
  {
    int32_t nv = 0;
    for (int32_t i = 0; i < 3; i++) 
      { for (int32_t j = 0; j < 3; j++) 
          { double *Rp = &(R->c[i][j]); /* Pointer to element of {*R}. */
            assert((*Rp) >= 0.0);
            if ((*Rp) != 0)
              { demand(nv <= 8, "at least one element of {M} must be fixed}");
                Mp[nv] = &(M->c[i][j]); /* Pointer to element of {*M}. */
                Re[nv] = (*Rp);
                nv++;
              }
          }
      }
    (*nvP) = nv;
  }

