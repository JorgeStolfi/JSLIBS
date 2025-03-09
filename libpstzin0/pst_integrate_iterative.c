/* See pst_integrate_iterative.h */
/* Last edited on 2025-03-02 14:16:10 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include <float_image.h>
#include <float_image_mscale.h>
#include <jsfile.h>
#include <r2.h>
#include <jswsize.h>
#include <rn.h>

#include <pst_basic.h>
#include <pst_imgsys.h>
#include <pst_imgsys_solve.h>
#include <pst_interpolate.h>
#include <pst_slope_map.h>
#include <pst_weight_map.h>
#include <pst_height_map.h>
#include <pst_integrate.h>

#include <pst_integrate_iterative.h>

void pst_integrate_iterative
  ( float_image_t *G, 
    float_image_t *H,
    double hintsWeight,
    float_image_t *Z,
    float_image_t *R,
    uint32_t maxIter,
    double convTol,
    bool_t sortSys,
    bool_t verbose,
    int32_t level,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  )
  {
    int32_t indent = (level < -1 ? 0 : 2*level+2);
    /* Get input image size and check compatibility: */
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 3, "slope map {G} must have 3 channels");
    if (verbose) { fprintf(stderr, "%*sinput slope map size = %d×%d ...\n", indent, "", NX_G, NY_G); }
  
    /* Compute output image size: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    if (verbose) { fprintf(stderr, "%*sheight map size = %d×%d ...\n", indent, "", NX_Z, NY_Z); }
    demand(Z != NULL, "height map must not be null");
    float_image_check_size(Z, 2, NX_Z, NY_Z, "bad height map");

    /* Build the system: */
    if (verbose) { fprintf(stderr, "%*sbuilding the system ...\n", indent, ""); }
    pst_imgsys_t *S = pst_integrate_build_system(G, H, hintsWeight, indent, verbose);
    if (verbose) { fprintf(stderr, "%*sfinal system has %d equations and %d variables\n", indent, "", S->N, S->N); }
 
    if (reportSys != NULL) { reportSys(level, S); }

    /* Get the initial guess {z} for the heights: */
    double *z = talloc(S->N, double);
    
    auto void get_initial_guess(void);
      /* Sets {z[k]} to the initial guess of variable {k}, derived from {Z}. */
      
    auto void reportSol(int32_t level, int32_t iter, double change, bool_t final, uint32_t N, double z[]);

    get_initial_guess();
    if (verbose) { fprintf(stderr, "%*ssolving the system ...\n", indent, ""); }
    bool_t para = FALSE; /* Should be parameter. 1 means parallel execution, 0 sequential. */
    bool_t szero = TRUE; /* Should be parameter. 1 means adjust sum to zero, 0 let it float. */
    uint32_t *ord = NULL;
    if (sortSys) { ord = pst_imgsys_sort_equations(S); }
    pst_imgsys_solve_iterative
      ( S, z, ord, maxIter, convTol, 
        para, szero, verbose, level, 
        reportStep,
        (reportHeights == NULL ? NULL : &reportSol)
      );
    pst_imgsys_copy_sol_vec_to_image(S, z, Z, 0.0);
        
    free(z);
    pst_imgsys_free(S);
    if (ord != NULL) { free(ord); }
    return;

    void reportSol(int32_t level, int32_t iter, double change, bool_t final, uint32_t N, double z[])
      { assert(reportHeights != NULL);
        pst_imgsys_copy_sol_vec_to_image(S, z, Z, 0.0);
        reportHeights(level, iter, change, final, Z, R);
      }

    void get_initial_guess(void)
      { if (verbose) { fprintf(stderr, "%*sobtaining initial guess of heights from {Z} ...\n", indent, ""); }
        for (int32_t k = 0; k < S->N; k++)
          { /* Index of equation's main variable: */
            int32_t uidk = k;
            /* Compute indices of pixel corresponding to the variable {Z[k]}: */
            int32_t xk = S->col[uidk]; int32_t yk = S->row[uidk];
            assert((xk >= 0) && (xk < NX_Z) && (yk >= 0) && (yk < NY_Z));
            double wzk = float_image_get_sample(Z, 1, xk, yk);
            z[uidk] = (wzk == 0.0 ? 0.0 : float_image_get_sample(Z, 0, xk, yk));
          }
      }
  }
