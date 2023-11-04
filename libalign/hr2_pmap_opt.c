/* See {hr2_pmap_opt.h}. */
/* Last edited on 2023-10-20 18:27:04 by stolfi */

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
#include <hr2_pmap_opt_translation.h>
#include <hr2_pmap_opt_congruence.h>
#include <hr2_pmap_opt_similarity.h>
#include <hr2_pmap_opt_affine.h>
#include <hr2_pmap_opt_generic.h>
 
int32_t hr2_pmap_opt_degrees_from_type(hr2_pmat_type_t type);
  /* Number of degrees of freedom in a projective map of the specified {type}. */

void hr2_pmap_opt_check_type(hr2_pmap_t *M, hr2_pmat_type_t type);
  /* Checks that the map {*A} is of the given type. */

void hr2_pmap_opt_encode(hr2_pmap_t *M, hr2_pmat_type_t type, int32_t nv, double  rad[], double y[]);
  /* Converts the map {*A} into  the optmizaton variables {y[0..nv-1]}. */

void hr2_pmap_opt_decode(int32_t nv, double y[], double  rad[], hr2_pmat_type_t type, hr2_pmap_t *M);
  /* Converts the optimization variables {y[0..nv-1]} into 
    the projective map {*A}. */

void hr2_pmap_opt_aff_quadratic
  ( hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    hr2_pmat_type_t type,     /* Desired map type. */
    int32_t max_iter,         /* Max outer loop iterations. */
    double max_f2,            /* Stopping criterion. */
    double max_mod,           /* Maximum change allowed in each matrix element. */
    hr2_pmap_t *A,            /* (IN/OUT) The affine map to adjust. */
    double *f2A_P             /* (OUT) Goal function value for the output {*A}. */
  )
  {
    bool_t debug = TRUE;
    
    demand(type != hr2_pmat_type_PROJECTIVE, "cannot be used for generic proj map");
    
    /* Identify the elements to optimize: */
    int32_t nv = hr2_pmap_opt_degrees_from_type(type);
    fprintf(stderr, "%d degrees of freedom in map\n", nv);
    assert (nv <= 6);
   
    /* Ensure that the map {A} is of the proper type: */
    hr2_pmap_opt_check_type(A, type);
    double Fz = f2(A); /* Current goal function on {A}. */
    (*f2A_P) = Fz;
     
    if ((Fz > max_f2) && (nv > 0))
      { 
        /* Apply nonlinear optimization. */

        /* Initialize constant elements of the incremental map {M}: */
        hr2_pmap_t M = hr2_pmap_identity(); 
        
        /* Scaling factors from the unit cube to the mod range: */
        double rad[nv];
        for (int32_t iv = 0; iv < nv; iv++) { rad[iv] = max_mod; }
            
        /* Save the initial map: */
        hr2_pmap_t A0 = (*A);
        
        auto double sve_goal(int32_t ny, double y[]);
          /* Converts the optimization variables {y[0..ny-1]} to the map {*A}, then computes {f2(A)}.
            Expects {ny == nv}. */
        
        /* Convert given map tp optimization variables: */
        double z[nv]; 
        hr2_pmap_opt_encode(A, type, nv, rad, z);

        sign_t dir = -1; /* Look for minimum. */
        sve_minn_iterate
          ( nv, 
            &sve_goal, NULL, 
            z, &Fz,
            dir,
            /*dMax:*/ 1.0,  
            /*dBox:*/ TRUE, 
            /*rIni:*/ 0.5,  
            /*rMin:*/ 0.01, 
            /*rMax:*/ 0.5,  
            /*stop:*/ 0.0001, 
            max_iter,
            debug
          );

        /* Ensure that the optimum {z} is reflected in {A,f2A): */
        Fz = sve_goal(nv, z);
        (*f2A_P) = Fz;

        /* Local implementations: */

        double sve_goal(int32_t ny, double y[])
          { assert(ny == nv);
            /* Convert variables {y[0..ny-1]} to a modification {*A}: */
            hr2_pmap_opt_decode(nv, y, rad, type, &M);
            (*A) = hr2_pmap_compose(&A0, &M);
            /* Evaluate the client function: */
            double Q2 = f2(A);
            (*f2A_P) = Q2;
            return Q2;
          }
      }
    assert(Fz == (*f2A_P));
    return;
  }

int32_t hr2_pmap_opt_degrees_from_type(hr2_pmat_type_t type)
  {
    switch(type)
      { 
        case hr2_pmat_type_TRANSLATION:
          return 2;
        case hr2_pmat_type_CONGRUENCE:
          return 3;
        case hr2_pmat_type_SIMILARITY:
          return 4;
        case hr2_pmat_type_AFFINE:
          return 6;
        case hr2_pmat_type_PROJECTIVE:
          demand(FALSE, "generic projective map not supported");
        default:
          demand(FALSE, "unimplemented map type");
      }
  }

void hr2_pmap_opt_check_type(hr2_pmap_t *M, hr2_pmat_type_t type)
  { 
    fprintf(stderr, "!! {hr2_pmap_opt_check_type} NOT IMPLEMENTED"); 
  }

void hr2_pmap_opt_encode(hr2_pmap_t *M, hr2_pmat_type_t type, int32_t nv, double y[], double rad[])
  {
    /* Assumes that {M} is a positive matrix with {M[0,0] = 1}, near the identity. */
    switch(type)
      { 
        case hr2_pmat_type_TRANSLATION:
          { assert(nv == 2);
            hr2_pmap_opt_translation_encode(M, rad, y);
            return;
          }
        case hr2_pmat_type_CONGRUENCE:
          { assert(nv == 3);
            hr2_pmap_opt_congruence_encode(M, rad, y);
            return;
          }
        case hr2_pmat_type_SIMILARITY:
          { assert(nv == 4);
            hr2_pmap_opt_similarity_encode(M, rad, y);
            return;
          }
        case hr2_pmat_type_AFFINE:
          { assert(nv == 6);
            hr2_pmap_opt_affine_encode(M, rad, y);
            return;
          }
        case hr2_pmat_type_PROJECTIVE:
          demand(FALSE, "generic projective map not supported");
        default:
          demand(FALSE, "unimplemented map type");
      }
  }

void hr2_pmap_opt_decode(int32_t nv, double y[], double rad[], hr2_pmat_type_t type, hr2_pmap_t *M)
  { 
    /* Assumes that the non-variable elements of {M} are those of the identity matrix. */
    switch(type)
      { 
        case hr2_pmat_type_TRANSLATION:
          { assert(nv == 2);
            hr2_pmap_opt_translation_decode(y, rad, M);
            return;
          }
        case hr2_pmat_type_CONGRUENCE:
          { assert(nv == 3);
            hr2_pmap_opt_congruence_decode(y, rad, M);
            return;
          }
        case hr2_pmat_type_SIMILARITY:
          { assert(nv == 4);
            hr2_pmap_opt_similarity_decode(y, rad, M);
            break;
          }
        case hr2_pmat_type_AFFINE:
          { assert(nv == 6);
            hr2_pmap_opt_affine_decode(y, rad, M);
            break;
          }
        case hr2_pmat_type_PROJECTIVE:
          demand(FALSE, "generic projective map not supported");
        default:
          demand(FALSE, "unimplemented map type");
      }
  }
 
