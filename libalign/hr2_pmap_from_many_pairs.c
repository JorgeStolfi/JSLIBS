/* See {hr2_pmap_from_many_pairs.h}. */
/* Last edited on 2023-10-09 19:36:03 by stolfi */

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

#include <hr2_pmap_from_many_pairs.h>

hr2_pmap_t hr2_pmap_from_many_pairs
  ( int32_t np,
    r2_t p1[], 
    r2_t h1[],
    r2_t p2[], 
    r2_t h2[],
    hr2_pmat_type_t type,
    int32_t maxIter,
    double maxErr,
    hr2_pmap_opt_report_proc_t report
  )
  {
    bool_t debug = TRUE;
    if (debug) fprintf(stderr, "  > begin {hr2_pmap_from_many_pairs}\n");
    
    hr2_pmap_t M0 = hr2_pmap_from_many_pairs_initial(type, np, p1, p2);
                                        
    if (np < 4)
      { fprintf(stderr, "too few point pairs (%d) for projective map, using affine map\n", np); 
        return (*M0);
      }
    
    int32_t nx = 8; /* Number of packed parameters. */
    
    auto hr2_pmap_t compute_S_map(void);
      /* Computes the {S} map used by {unpack_parameters}. */
         
    hr2_pmap_t S = compute_S_map();  /* Pmap from image1's domain to canonical square. */
    
    auto void pack_parameters(hr2_pmap_t *M, int32_t nx, double x[]);
      /* Packs the homogeneous projective map {M} to the parameter
        vector {x[0..nx-1]}, as the images of the four corners 
        of image 1. Requires {nx==8}. */
         
    auto void unpack_parameters(int32_t nx, double x[], hr2_pmap_t *M);
      /* Unpacks the parameter vector {x[0..nx-1]} to the 
        homogeneous map matrix {M}.  The matrix
        will have {M[0,0] == 1}. Requires {nx==8}. */
         
    auto double goalf(int32_t nx, double x[]);
      /* Computes the mean squared distance between the positions of the
         mapped points {p1*M.dir} and {p2*M.inv}, given the packed
         parameters {x[0..nx-1]}. */
    
    auto bool_t is_ok(int32_t nx, double x[], double Fx);
      /* Returns true if the squared error {Fx} is small enough. */
    
    double esq = image_stitch_mean_err_sqr(p1, &(M0->dir), p2, &(M0->inv));
    double dMax = 20*sqrt(esq);      /* Search radius around initial guess. */
    bool_t dBox = FALSE;             /* Search in ball, not in box. ??? */
    double rIni = 0.250*dMax;        /* Initial probe radius. */
    double rMin = 0.5;               /* Minimum probe radius. */
    double rMax = 0.500*dMax;        /* Maximum probe radius. */
    double stop = 0.1*maxErr;        /* Stopping criterion. */
    sign_t dir = -1;                 /* Look for minimum. */
    
    if (verbose) 
      { fprintf(stderr, "initial mean error squared = %22.16e\n", esq);
        fprintf(stderr, "estimated distance from optimum = %13.6f\n", dMax);
        fprintf(stderr, "probe radius = %13.6f [ %13.6f _ %13.6f ]\n", rIni, rMin, rMax);
      }
      
    double x[nx]; /* Initial guess and final optimum parameters. */ 
    pack_parameters(M0, nx, x);
    if (verbose) 
      { char *fname = NULL;
        asprintf(&fname, "%s-f2-plot.txt", outPrefix);
        FILE *wr = open_write(fname, TRUE);
        minn_plot_1D_gnuplot(wr, nx, goalf, x, 20, rIni);
        fclose(wr);
        free(fname);
      }
    if (verbose) { fprintf(stderr, "optimizing\n"); }
    double Fx = goalf(nx, x);
    if (verbose) { fprintf(stderr, "initial rms error = %13.6f\n", Fx); }
    sve_minn_iterate(nx, &goalf, &is_ok, x, &Fx, dir, dMax, dBox, rIni, rMin, rMax, stop, maxIter, verbose);
    if (verbose) { fprintf(stderr, "final rms error = %13.6f\n", Fx); }
    
    /* Unpack projective map: */
    hr2_pmap_t M; /* Optimized map. */
    unpack_parameters(nx, x, &M);
    if (verbose) { fprintf(stderr, "--- exit image_stitch_optimize_pmap ---\n"); }
    return M;
    
    /* --- internal procs ----------------------------------------------------- */
    
    hr2_pmap_t compute_S_map(void)
      { hr2_point_t ph[4];
        int32_t ix, iy;
        for (int32_t ix = 0; ix < 2; ix++)
          { double px = (ix == 0 ? L1->c[0] : H1->c[0]);
            for (int32_t iy = 0; iy < 2; iy++)
              { double py = (iy == 0 ? L1->c[1] : H1->c[1]);
                ph[2*iy+ix] = (hr2_point_t){{{ 1, px, py }}};
              }
          }
        hr2_pmap_t S = hr2_pmap_from_four_points(&(ph[0]), &(ph[1]), &(ph[2]), &(ph[3]));
        return hr2_pmap_inv(&S);
      }
    
    void pack_parameters(hr2_pmap_t *M, int32_t nx, double x[])
      { assert(nx == 8);
        int32_t ix, iy;
        int32_t k = 0;
        for (int32_t ix = 0; ix < 2; ix++)
          { double px = (ix == 0 ? L1->c[0] : H1->c[0]);
            for (int32_t iy = 0; iy < 2; iy++)
              { double py = (iy == 0 ? L1->c[1] : H1->c[1]);
                r2_t p = (r2_t){{ px, py }};
                r2_t q; image_stitch_map_point(&p, &(M->dir), &q);
                x[k] = q.c[0]; k++;
                x[k] = q.c[1]; k++;
              }
          }
      }
    
    void unpack_parameters(int32_t nx, double x[], hr2_pmap_t *M)
      { assert(nx == 8);
        int32_t ix, iy;
        int32_t k = 0;
        hr2_point_t ph[4];
        for (int32_t ix = 0; ix < 2; ix++)
          { for (int32_t iy = 0; iy < 2; iy++)
              { double qx = x[k]; k++;
                double qy = x[k]; k++;
                ph[2*iy + ix] = (hr2_point_t){{{ 1, qx, qy }}};
              }
          }
        hr2_pmap_t R = hr2_pmap_from_four_points(&(ph[0]), &(ph[1]), &(ph[2]), &(ph[3]));
        (*M) = hr2_pmap_compose(&S, &R); 
      }
    
    double goalf(int32_t nx, double x[])
      {
        if (debug) { rn_gen_print(stderr, nx, x, "%8.1f", "\n    [ ", " ", " ]\n"); }
        hr2_pmap_t M;
        unpack_parameters(nx, x, &M);
        /* Mean square error on point lists: */
        double esq = hr2_pmap_mismatch_sqr(&M, np, p1, p2);
        /* Deformation penalty: */
        double dsq1 = image_stitch_deform_sqr(L1, H1, &(M.dir));
        double dsq2 = image_stitch_deform_sqr(L2, H2, &(M.inv));
        double F = esq + 0.0001*(dsq1 + dsq2);
        if (debug) 
          { fprintf(stderr, "    mean squared error = %13.6e", esq); 
            fprintf(stderr, " squared deform = %13.6e %13.6e", dsq1, dsq2);
            fprintf(stderr, " function = %13.6e\n", F);
          }
        return F;
        
      }
      
    bool_t is_ok(int32_t nx, double x[], double Fx)
      {
        return Fx < maxErr*maxErr;
      }
    
    return hr2_pmap_identity();
  }

void hr2_pmap_from_many_pairs_quad
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
    hr2_pmap_from_many_pairs_get_var_elems(&(A->dir), R, Ap, Re, &nv);
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
