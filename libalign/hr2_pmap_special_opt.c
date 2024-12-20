/* See {hr2_pmap_special_opt.h}. */
/* Last edited on 2024-12-05 10:19:51 by stolfi */

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
#include <rmxn_extra.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <affirm.h>
#include <minn.h>
#include <minn_plot.h>
#include <sve_minn.h>

#include <hr2_pmap_encode.h>
#include <hr2_pmap_opt.h>

#include <hr2_pmap_special_opt.h>

double hr2_pmap_special_opt_eval_encoded
  ( int32_t n,
    double z[],
    double yctr[],
    double yrad[],
    hr2_pmap_type_t type,
    sign_t sgn, 
    hr2_pmap_t *M,
    hr2_pmap_opt_func_t *f2
  );
  /* Converts {z[0..n-1]} to {y[0..n-1]} by the formula {y[k] = yctr[k]
    + yrad[k]*z[k]} for {k} in {0..n-1}. Then performs
    {hr2_pmap_decode(n, y, type, M, sgn)} and then evaluates {f2} on
    {M}. If {type} is general projective, also adds a small multiple of
    {hr2_pmap_opt_homo_scaling_bias(M)}. */
  
void hr2_pmap_special_opt_quadratic
  ( hr2_pmap_type_t type,     /* Desired map type. */
    sign_t sgn,               /* Handedness, when applicable. */
    hr2_pmap_opt_func_t *f2,  /* Goal function to minimize. */
    hr2_pmap_opt_pred_t *OK,  /* Client acceptance predicate. */
    double yrad[],            /* Scaling factor for each encoding element. */
    int32_t maxIter,          /* Max outer loop iterations. */
    hr2_pmap_t *M,            /* (IN/OUT) The affine map to adjust. */
    double *f2M_P,            /* (OUT) Goal function value for the output {*M}. */
    bool_t verbose
  )
  {
    bool_t debug = FALSE;
    
    if (verbose || debug) { fprintf(stderr, "    > --- %s ---\n", __FUNCTION__); }
    
    demand((sgn == +1) || (sgn == -1), "invalid map sign");
    demand(type != hr2_pmap_type_NONE, "invalid map type");
    demand(type != hr2_pmap_type_GENERIC, "not for general projective maps");

    /* Get the number of optimization variables {n} */
    int32_t n = hr2_pmap_encode_num_parameters(type);
    if (verbose) { fprintf(stderr, "    %d optimization variables\n", n); }
    assert (n <= 9);
   
    double f2z = f2(M); /* Current goal function on {M}. */
    if (verbose) { fprintf(stderr, "    initial f2(M) = %20.14f\n", f2z); }

    /* Ensure that the map {M} is of the proper type and handedness: */
    double tol = 1.0e-14; /* Tolerance for type match. */
    demand(hr2_pmap_is_type(M, type, sgn, tol), "wrong map type or sign");
     
    bool_t ok = (OK != NULL) && OK(M, f2z);
    if ((! ok) && (n > 0))
      { 
        /* Apply nonlinear optimization. */
        
        /* Get the center of the search region: */
        double yctr[n];
        hr2_pmap_encode(M, type, n, yctr);
        
        bool_t sve_debug = debug;
        bool_t sve_debug_probes = FALSE; 
            
        auto double sve_goal(int32_t n1, double z1[]);
          /* Converts the optimization variables {z1[0..n-1]} into a
            projective map, using {yctr} and {yrad}, and stores it into
            {*M}. Then computes {f2(M)}. Expects {n1 == n}. */
        
        auto bool_t sve_ok(int32_t iter, int32_t n1, double z1[], double f2z, double dist, double step, double radius);
          /* Returns {TRUE} iff the matrix defined by the optimization
            variables {z1[0..n-1]} is good enough. Specifically, converts
            the optimization variables {z1[0..n-1]} into a projective
            map, using {yctr} and {yrad}, and calls the client's {OK}
            predicate. Expects {OK != NULL} and {n1 == n}. */
        
        auto double sve_proj(int32_t n1, double z1[], double f2z);
          /* Normalizes the optimization variables {z1[0..n-1]}, in case
            there are multiple vectors {z} that produce equivalent maps
            {M}. Expects {n1 == n}.
            
            Specifically, converting them to an encoding vector
            {y[0..n-1]}, using {yctr} and {yrad}, performing an decode
            and encode on {y}, and converting the resulting {y} back to
            {z1[0..n-1]}. For encodings that include angles, for
            example, this process reduces those angles to the range
            {(-\pi _ +\pi]}.
            
            The procedure then recomputes and returns {sve_goal(z1)},
            just in case that normalization changed it (e.g. because of
            roundoff errors). */
        
        /* The optimization variables {(0...0)} mean the initial map {M}: */
        double z[n];
        for (uint32_t k = 0;  k < n; k++) { z[k] = 0.0; }

        sign_t dir = -1; /* Look for minimum. */
        double *ctr = NULL; /* Domain of {z} is centered at origin. */
        double dMax = 1.0;
        bool_t dBox = FALSE;
        double rIni = 0.5*dMax;
        double rMin = 0.0001;
        double rMax = 2.0*dMax;
        double minStep = 0.01*rMin;  /* Min {z} change between iterations. */
        
        if (verbose) 
          { fprintf(stderr, "  encoding center and radius:\n");
            for (uint32_t k = 0;  k < n; k++) { fprintf(stderr, "    %d %13.6f %13.6f\n", k, yctr[k], yrad[k]); }
            fprintf(stderr, "  estimated {z} distance from optimum = %13.6f\n", dMax);
            fprintf(stderr, "  probe radius = %13.6f [ %13.6f _ %13.6f ]\n", rIni, rMin, rMax);
          }

        sve_minn_iterate
          ( n, &sve_goal, (OK == NULL ? NULL : &sve_ok), &sve_proj,
            z, &f2z, dir, 
            ctr, dMax, dBox, rIni, rMin, rMax, 
            minStep, maxIter, sve_debug, sve_debug_probes
          );

        /* Ensure that the optimum {z} is reflected in {M,f2M): */
        f2z = sve_goal(n, z);
        (*f2M_P) = f2z;

        /* Local implementations: */

        double sve_goal(int32_t n1, double z1[])
          { assert(n1 == n);
            double Q2 = hr2_pmap_special_opt_eval_encoded(n, z1, yctr, yrad, type, sgn, M, f2);
            return Q2;
          }
          
        bool_t sve_ok(int32_t iter, int32_t n1, double z1[], double f2z, double dist, double step, double radius)
          { assert(OK != NULL);
            assert(n1 == n);
            double y[n];
            for (uint32_t k = 0;  k < n; k++) { y[k] = yctr[k] + z1[k]*yrad[k]; }
            hr2_pmap_t MM;
            hr2_pmap_decode(n, y, type, sgn, &MM);
            return OK(&MM, f2z);
          }
          
        double sve_proj(int32_t n1, double z1[], double f2z)
          { assert(n1 == n);
            hr2_pmap_t MM;
            double y[n];
            for (uint32_t k = 0;  k < n; k++) { y[k] = yctr[k] + z1[k]*yrad[k]; }
            hr2_pmap_decode(n, y, type, sgn, &MM);
            hr2_pmap_encode(&MM, type, n, y);
            for (uint32_t k = 0;  k < n; k++) { z1[k] = (y[k] - yctr[k])/yrad[k]; }
            double f2z_new = sve_goal(n, z1);
            if (fabs(f2z_new - f2z) > 1.0e-6*(fabs(f2z_new) + fabs(f2z)) + 1.0e-12) 
              { fprintf(stderr, "      !! F(M) changed by projection\n"); }
            return f2z_new;
          }
      }

    if (verbose) { fprintf(stderr, "    final f2(M) =   %20.14f\n", f2z); }
    (*f2M_P) = f2z;
    
    if (verbose || debug) { fprintf(stderr, "    < --- %s ---\n", __FUNCTION__); }
    return;
  }
  
double hr2_pmap_special_opt_eval_encoded
  ( int32_t n,
    double z[],
    double yctr[],
    double yrad[],
    hr2_pmap_type_t type,
    sign_t sgn, 
    hr2_pmap_t *M,
    hr2_pmap_opt_func_t *f2
  )
  {
    double y[n];
    for (uint32_t k = 0;  k < n; k++) { y[k] = yctr[k] + z[k]*yrad[k]; }
    hr2_pmap_decode(n, y, type, sgn, M);
    /* Evaluate the client function: */
    double Q2 = f2(M);
    return Q2;
  }
        
void hr2_pmap_special_opt_1D_plot
  ( char *outPrefix,
    char *tag,
    char *stage,
    hr2_pmap_t *M,
    hr2_pmap_type_t type,
    sign_t sgn,
    hr2_pmap_opt_func_t *f2,
    double yrad[],
    double urad,
    int32_t nu,
    int32_t ns
  )
  {
    demand(type != hr2_pmap_type_IDENTITY, "cannot vary maps of {type} indentity");
    
    /* Get the number of optimization variables {n} */
    int32_t n = hr2_pmap_encode_num_parameters(type);
    assert (n <= 9);
    
    double yctr[n];
    hr2_pmap_encode(M, type, n, yctr);

    hr2_pmap_t N = (*M);
    
    auto double sve_goal(int32_t n1, double z[]);
      /* Converts the optimization variables {z[0..n-1]} to the map
        {*N}, then computes {f2(N)}. Expects {n1 == n}. */

    
    double z[n];
    rn_zero(n, z);
    
    /* Choose the plot directions: */
    double U[nu*n]; /* The rows are the directions. */
    rmxn_throw_directions(nu, n, U);

    /* Plot the goal function along the directions: */
    char *fname = jsprintf("%s-%s-%s-1D-plot.txt", outPrefix, tag, stage);
    FILE *wr = open_write(fname, TRUE);
    double step = urad/ns;
    for (uint32_t ku = 0;  ku < nu; ku++)
      { double *uk = &(U[ku*n]);
        minn_plot_1D_gnuplot(wr, n, z, uk, urad, step, sve_goal);
      }
    fclose(wr);
    free(fname);

    double sve_goal(int32_t n1, double z1[])
      { assert(n1 == n);
        double Q2 = hr2_pmap_special_opt_eval_encoded(n1, z1, yctr, yrad, type, sgn, &N, f2);
        return Q2;
      }
  }
