/* See {r2_aff_adjust.h}. */
/* Last edited on 2020-10-16 03:38:00 by jstolfi */

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
#include <sve_minn.h>

#include <r2_aff_adjust.h>

/* INTERNAL PROTOTYPES */

void r2_aff_adjust_get_var_elems
  ( r2_aff_map_t *A, 
    r2_aff_map_t *R, 
    double *Ap[],
    double Re[],
    int32_t *nvP
  );
  /* Identifies the elements of {*A} that are to be adjusted 
    (that is, whose corresponding element of {*R} is nonzero.
    Returns the number {nv} of such elements (a number in {0..6}) in {*nvP},
    the addresses of those elemens of {*A} in {Ap[0..nv-1]},
    and the values of the corresponding elements of {*R}
    in {Re[0..nv-1]}.  These vectors should have at leat 6 elements. */

/* IMPLEMENTATIONS */

void r2_aff_adjust_quad
  ( r2_aff_adjust_func_t *f2,  /* Goal function to minimize. */
    r2_aff_map_t *R,     /* Max adjustment for {*A}. */
    double tol,          /* Desired relative adjustment precision for {*A}. */
    r2_aff_map_t *A,     /* (IN/OUT) The affine map to adjust. */
    double *f2P          /* (OUT) Goal function value for the output {*A}. */
  )
  {
    int32_t maxIters = 10;
    bool_t debug = TRUE;
    
    (*f2P) = f2(A); /* Current goal function on {A}. */

    /* Identify the elements to optimize: */
    int32_t nv;
    double *Ap[6]; /* Pointers to adjustable elements. */
    double Re[6];  /* Max adjustment amount of each element. */
    r2_aff_adjust_get_var_elems(A, R, Ap, Re, &nv);
    fprintf(stderr, "%d adjustable elements found\n", nv);
    
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

        /* Convert the affine maps {*A} to the optimization variables {z[0..nv-1]}: */
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

        /* Convert the optimization variables {z[0..nv-1]} to affine maps {*A}: */
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
          }

        double sve_goal(int32_t nx, double x[])
          { assert(nx == nv);
            /* Convert variables {x[0..nx-1]} to displacements {*A}: */
            store_vars(x);
            /* Evaluate the client function: */
            double Q2 = f2(A);
            (*f2P) = Q2;
            return Q2;
          }
      }

    /* Compue the final mismatch: */
    (*f2P) = f2(A);
    return;
  }
    
void r2_aff_adjust_get_var_elems
  ( r2_aff_map_t *A, 
    r2_aff_map_t *R, 
    double *Ap[],
    double Re[],
    int32_t *nvP
  )
  {
    int32_t nv = 0;
    for (int32_t s = 0; s < 6; s++) 
      { double *Rp = r2_aff_map_elem_addr(R, s); /* Pointer to element of {*R}. */
        assert((*Rp) >= 0.0);
        if ((*Rp) != 0)
          { Ap[nv] = r2_aff_map_elem_addr(A, s); /* Pointer to element of {*A}. */
            Re[nv] = (*Rp);
            nv++;
          }
      }
    (*nvP) = nv;
  }


