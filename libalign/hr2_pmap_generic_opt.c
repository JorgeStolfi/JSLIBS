/* See {hr2_pmap_generic_opt.h}. */
/* Last edited on 2025-02-16 20:20:42 by stolfi */

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

#include <hr2_pmap_generic_opt.h>

void hr2_pmap_generic_opt_map_to_vars(hr2_pmap_t *M, uint32_t ny, double y[]);
  /* Copies the elements of {M.dir}, in row order, to {y[0..ny-1]}.
    Expects that {M} has the sign {sng} and {ny==9}. */

void hr2_pmap_generic_opt_vars_to_mat(uint32_t ny, double y[], sign_t sgn, hr2_pmap_t *M);
  /* Copies {y[0..ny-1]} to the the elements of {M.dir}, in row order.
    Then computes {M.inv}. Then forces the sign of {M}
    to be {sgn} by {hr2_pmap_set_sign} if necessary.
    Expects that {ny==9}. */

double hr2_pmap_generic_opt_eval(hr2_pmap_t *M, hr2_pmap_opt_func_t *f2);
  /* Evaluates {f2} on {M}, then adds a small multiple of 
    {hr2_pmap_opt_homo_scaling_bias(M)}. */
  
void hr2_pmap_generic_opt_quadratic
  ( sign_t sgn,               /* Handedness, when applicable. */
    hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    hr2_pmap_opt_pred_t *OK,  /* Client acceptance predicate. */
    int32_t maxIter,          /* Max outer loop iterations. */
    double maxMod,            /* Maximum change allowed in each matrix element. */
    hr2_pmap_t *M,            /* (IN/OUT) The projective map to adjust. */
    double *f2M_P,            /* (OUT) Goal function value for the output {*M}. */
    bool_t verbose
  )
  {
    bool_t debug = FALSE;
    
    if (verbose || debug) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    
    demand((sgn == +1) || (sgn == -1), "invalid map sign");
    demand(hr2_pmap_sign(M) == sgn, "initial guess has the wrong sign");

    /* The optimization variables are the elements of {M.dir}: */
    uint32_t nz = 9;
    double z[nz];

    hr2_pmap_generic_opt_map_to_vars(M, nz, z);
        
    double f2z = f2(M); /* Current goal function on {M}. */
    if (verbose) { fprintf(stderr, "    initial f2(M) = %20.14f\n", f2z); }

    /* Ensure that the map {M} is of the proper handedness: */
    demand(hr2_pmap_sign(M), "wrong map sign");
     
    bool_t ok = (OK != NULL) && OK(M, f2z);
    if (! ok)
      { 
        /* Apply nonlinear optimization. */
        
        bool_t sve_debug = debug;
        bool_t sve_debug_probes = FALSE; 
            
        auto double sve_goal(int32_t ny, double y[]);
          /* Sets the elements of {M.dir} to {y[0..ny-1]},
            computes {M.inv}, then computes {f2(M)}. 
            Expects {ny == nz}. */
        
        auto bool_t sve_ok(int32_t iter, int32_t ny, double y[], double f2y, double dist, double step, double radius);
          /* Assumes that {y[0..ny-1]} are 
            the elements of {M.dir}, {M.inv} is the 
            inverse of {M.dir}, and {f2y} is {f2(M}}. Returns {TRUE} iff 
            {OK(M)} is true.  Expects {OK != NULL} and {ny == nz}. */
        
        auto double sve_proj(int32_t ny, double y[], double f2y);
          /* Normalizes {y[0..n-1]} to unit norm, then performs
            {sve_goal(y)} to update {M.dir}, {M.inv} and recompute {f2(M)}, which
            should be the same as {f2y}. Expects {ny == nz}. */
        
        sign_t dir = -1; /* Look for minimum. */
        double ctr[nz]; rn_copy(nz, z, ctr);
        double dMax = maxMod * sqrt(nz);
        bool_t dBox = FALSE;
        double rIni = 0.5*dMax;
        double rMin = 0.0001;
        double rMax = 2.0*dMax;
        double minStep = 0.01*rMin;  /* Min {z} change between iterations. */
        
        if (verbose) 
          { fprintf(stderr, "  max matrix element change = %13.6f\n", maxMod);
            fprintf(stderr, "  estimated distance from optimum = %13.6f\n", dMax);
            fprintf(stderr, "  probe radius = %13.6f [ %13.6f _ %13.6f ]\n", rIni, rMin, rMax);
          }

        sve_minn_iterate
          ( nz, &sve_goal, (OK == NULL ? NULL : &sve_ok), &sve_proj,
            z, &f2z, dir, 
            ctr, dMax, dBox, rIni, rMin, rMax, 
            minStep, maxIter, sve_debug, sve_debug_probes
          );

        /* Ensure that the optimum {z} is reflected in {M,f2z): */
        f2z = sve_goal(nz, z);

        /* Local implementations: */

        double sve_goal(int32_t ny, double y[])
          { assert(ny == nz);
            hr2_pmap_generic_opt_vars_to_mat(ny, y, sgn, M);
            double Q2 = hr2_pmap_generic_opt_eval(M, f2);
            return Q2;
          }
          
        bool_t sve_ok(int32_t iter, int32_t ny, double y[], double f2y, double dist, double step, double radius)
          { assert(OK != NULL);
            assert(ny == nz);
            /* Assumes {M.inv} is the inverse of {M.dir} and {f2y} is {f2(M)}: */
            return OK(M, f2y);
          }
          
        double sve_proj(int32_t ny, double y[], double f2y)
          { assert(ny == nz);
            double mag = rn_dir(ny, y, y);
            assert(mag > 1.0e-18);
            double f2y_new = sve_goal(ny, y);
            if (fabs(f2y_new - f2y) > 1.0e-6*(fabs(f2y_new) + fabs(f2y)) + 1.0e-12) 
              { fprintf(stderr, "      !! F(M) changed by projection\n"); }
            return f2y_new;
          }
      }

    if (verbose) { fprintf(stderr, "    final f2(M) =   %20.14f\n", f2z); }
    (*f2M_P) = f2z;
    
    if (verbose || debug) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
    return;
  }

void hr2_pmap_generic_opt_map_to_vars(hr2_pmap_t *M, int32_t ny, double y[])
  { assert(ny == 9);
    int32_t ky = 0;
    for (uint32_t i = 0;  i < 3; i++)
      { for (uint32_t j = 0;  j < 3; j++)
          { y[ky] = M->dir.c[i][j]; ky++; }
      }
    assert(ky == ny);
    r3x3_inv(&(M->dir), &(M->inv));
  }

void hr2_pmap_generic_opt_vars_to_mat(int32_t ny, double y[], sign_t sgn, hr2_pmap_t *M)
  { assert(ny == 9);
    int32_t ky = 0;
    for (uint32_t i = 0;  i < 3; i++)
      { for (uint32_t j = 0;  j < 3; j++)
          { M->dir.c[i][j] = y[ky]; ky++; }
      }
    assert(ky == ny);
    hr2_pmap_set_sign(M, sgn);
  }
        
void hr2_pmap_generic_opt_1D_plot
  ( char *outPrefix,
    char *tag,
    char *stage,
    hr2_pmap_t *M,
    sign_t sgn,
    hr2_pmap_opt_func_t *f2,
    int32_t nu,
    double urad,
    int32_t ns
  )
  {
    hr2_pmap_t N = (*M);
    
    demand(nu <= 7, "too many directions");
    
    auto double sve_goal(int32_t ny, double y[]);
      /* Stores the variables {y[0..ny-1]} into the map {N.dir},
        computes {, then computes {f2(N)}.
        Expects {ny == 9}. */

    /* The optimization variables are the elements of {N.dir}: */
    int32_t nz = 9;
    double z[nz];
    hr2_pmap_generic_opt_map_to_vars(&N, nz, z);
    double zmag = rn_dir(nz, z, z);
    assert(zmag > 1.0e-18);
    
    /* Choose the plot directions and make them orthogonal to {z}: */
    double U[nu*nz]; /* The rows are the directions. */
    rmxn_throw_directions(nu, nz, U);
    for (uint32_t ku = 0;  ku < nu; ku++)
      { double *uk = &(U[ku*nz]);
        /* Project {uk} perpendicular to {z}: */
        while (fabs(rn_dot(nz, uk, z)) > 0.90)
          { /* Too close to {±z}, get another one: */
            rn_throw_dir (nz, uk); 
          }
        (void)rn_decomp(nz, uk, z, NULL, uk);
        double umag = rn_dir(nz, uk, uk);
        assert(umag > 0.0001);
      }

    /* Plot the goal function along the directions: */
    char *fname = jsprintf("%s-%s-%s-1D-plot.txt", outPrefix, tag, stage);
    FILE *wr = open_write(fname, TRUE);
    double step = urad/ns;
    for (uint32_t ku = 0;  ku < nu; ku++)
      { double *uk = &(U[ku*nz]);
        minn_plot_1D_gnuplot(wr, nz, z, uk, urad, step, sve_goal);
      }
    fclose(wr);
    free(fname);
    return;

    double sve_goal(int32_t ny, double y[])
      { assert(ny == nz);
        hr2_pmap_generic_opt_vars_to_mat(ny, y, sgn, &N);
        double Q2 = hr2_pmap_generic_opt_eval(&N, f2);
        return Q2;
      }
  }

double hr2_pmap_generic_opt_eval(hr2_pmap_t *M, hr2_pmap_opt_func_t *f2)
  {
    /* Evaluate the client function: */
    double Q2 = f2(M);
    double B2 = hr2_pmap_opt_homo_scaling_bias(M);
    return Q2 + 1.0e-6*B2;
  
  }
