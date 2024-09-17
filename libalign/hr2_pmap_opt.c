/* See {hr2_pmap_opt.h}. */
/* Last edited on 2024-09-17 13:42:01 by stolfi */

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
#include <rn.h>
#include <r3x3.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <affirm.h>
#include <minn.h>
#include <minn_plot.h>
#include <sve_minn.h>

#include <hr2_pmap_encode.h>

#include <hr2_pmap_opt.h>

double hr2_pmap_opt_eval_encoded
  ( int32_t ny,
    double y[],
    hr2_pmap_type_t type,
    sign_t sgn, 
    hr2_pmap_t *M,
    hr2_pmap_opt_func_t *f2
  );
  /* Performs {hr2_pmap_decode(ny, y, type, M, sgn)} and then
    evaluates {f2} on {M}.  If {type} is general projective, also adds 
    a small multiple of {hr2_pmap_opt_homo_scaling_bias(M)}. */

double hr2_pmap_opt_homo_scaling_bias(hr2_pmap_t *M);
  /* Returns {f2(M.dir)+f2(M.inv)}, where {f2} is a function of a 3x3
    matrix that is minimum when the sum of the squares of the elements
    is 1.
    
    This function can be used as a bias term when optimizing a general
    projective map, to remove the spurious degree of freedom due to
    homogeneous scaling of the matrices {M.dir} and {M.inv}. */
  
void hr2_pmap_opt_quadratic
  ( hr2_pmap_type_t type,     /* Desired map type. */
    sign_t sgn,               /* Handedness, when applicable. */
    hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    int32_t maxIter,          /* Max outer loop iterations. */
    double maxMod,            /* Maximum change allowed in each matrix element. */
    double f2Stop,            /* Stopping criterion. */
    hr2_pmap_t *M,            /* (IN/OUT) The affine map to adjust. */
    double *f2A_P,            /* (OUT) Goal function value for the output {*M}. */
    bool_t verbose
  )
  {
    bool_t debug = TRUE;
    
    if (verbose || debug) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }

    /* Get the number of optimization variables {nz} */
    int32_t nz = hr2_pmap_encode_num_parameters(type);
    if (verbose) { fprintf(stderr, "    %d optimization variables\n", nz); }
    assert (nz <= 9);
   
    double f2z = f2(M); /* Current goal function on {M}. */
    if (verbose) { fprintf(stderr, "    initial f2(M) = %20.14f\n", f2z); }

    /* Ensure that the map {M} is of the proper type and handedness: */
    double tol = 1.0e-14; /* Tolerance for type match. */
    demand(hr2_pmap_is_type(M, type, sgn, tol), "wrong map type or sign");
     
    if ((f2z > f2Stop) && (nz > 0))
      { 
        /* Apply nonlinear optimization. */
            
        auto double sve_goal(int32_t ny, double y[]);
          /* Converts the optimization variables {y[0..ny-1]} to the map {*M}, then computes {f2(M)}.
            Expects {ny == nz}. */
        
        auto bool_t sve_ok(int32_t ny, double y[], double f2y);
          /* Returns {TRUE} iff the matrix defined by {y[0..ny-1]} is good enough.
            Assumes that {f2y} is the value of {sve_goal(ny, y)}. */
        
        /* Convert given map tp optimization variables: */
        double z[nz];
        hr2_pmap_encode(M, type, nz, z);

        sign_t dir = -1; /* Look for minimum. */
        double ctr[nz]; rn_copy(nz, z, ctr);
        double dMax = maxMod * sqrt(nz);
        bool_t dBox = FALSE;
        double rIni = 0.5*dMax;
        double rMin = 0.01;
        double rMax = 2.0*dMax;
        double stop = 0.0001;  /* Min {z} change between iterations. */

        if (verbose) 
          { fprintf(stderr, "  max matrix element change = %13.6f\n", maxMod);
            fprintf(stderr, "  estimated distance from optimum = %13.6f\n", dMax);
            fprintf(stderr, "  probe radius = %13.6f [ %13.6f _ %13.6f ]\n", rIni, rMin, rMax);
          }

        sve_minn_iterate
          ( nz, &sve_goal, &sve_ok, z, &f2z,
            dir, ctr, dMax, dBox, rIni, rMin, rMax, stop,
            maxIter,
            debug
          );

        /* Ensure that the optimum {z} is reflected in {M,f2A): */
        f2z = sve_goal(nz, z);
        (*f2A_P) = f2z;

        /* Local implementations: */

        double sve_goal(int32_t ny, double y[])
          { assert(ny == nz);
            double Q2 = hr2_pmap_opt_eval_encoded(ny, y, type, sgn, M, f2);
            return Q2;
          }
          
        bool_t sve_ok(int32_t ny, double y[], double f2y)
          { return f2y <= f2Stop; }
      }

    if (verbose) { fprintf(stderr, "    final f2(M) =   %20.14f\n", f2z); }
    (*f2A_P) = f2z;
    
    if (verbose || debug) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
    return;
  }
  
double hr2_pmap_opt_eval_encoded
  ( int32_t ny,
    double y[],
    hr2_pmap_type_t type,
    sign_t sgn, 
    hr2_pmap_t *M,
    hr2_pmap_opt_func_t *f2
  )
  {
    hr2_pmap_decode(ny, y, type, sgn, M);
   
    /* Evaluate the client function: */
    double Q2 = f2(M);
    if (type == hr2_pmap_type_GENERIC)
      { /* Add hoogeneous bias term: */
        Q2 += 0.001 * hr2_pmap_opt_homo_scaling_bias(M);
      }
    return Q2;
  }
  
double hr2_pmap_opt_homo_scaling_bias(hr2_pmap_t *M)
  { double sumA2 = r3x3_norm_sqr(&(M->dir));
    double sumB2 = r3x3_norm_sqr(&(M->inv));
    return 0.25*(sumA2 + 1/sumA2 + sumB2 + 1/sumB2);
  }
        
void hr2_pmap_opt_1D_plot
  ( char *outPrefix,
    char *tag,
    hr2_pmap_t *M,
    hr2_pmap_type_t type,
    sign_t sgn,
    hr2_pmap_opt_func_t *f2,
    int32_t nu,
    double urad,
    int32_t ns
  )
  {
    demand(type != hr2_pmap_type_IDENTITY, "cannot vary maps of {type} indentity");

    hr2_pmap_t N = (*M);
    
    auto double sve_goal(int32_t ny, double y[]);
      /* Converts the optimization variables {y[0..ny-1]} to the map {*N}, then computes {f2(N)}.
        Expects {ny == nz}. */

    /* Get the number of optimization variables {nz} */
    int32_t nz = hr2_pmap_encode_num_parameters(type);
    assert (nz <= 9);
    double z[nz];
    hr2_pmap_encode(M, type, nz, z);
 
    char *fname = NULL;
    asprintf(&fname, "%s-%s-1D-plot.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    double step = urad/ns;
    minn_plot_1D_gnuplot(wr, nz, z, urad, step, sve_goal);
    fclose(wr);
    free(fname);

    double sve_goal(int32_t ny, double y[])
      { assert(ny == nz);
        double Q2 = hr2_pmap_opt_eval_encoded(ny, y, type, sgn, &N, f2);
        return Q2;
      }
  }
