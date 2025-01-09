/* See pst_integrate_iterative.h */
/* Last edited on 2025-01-08 05:34:06 by stolfi */

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
  ( float_image_t *IG, 
    float_image_t *IW, 
    bool_t keepNull,
    float_image_t *IZ, 
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    float_image_t **OZP,
    float_image_t **OWP,
    bool_t verbose,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  )
  {
    /* Get input image size and check compatibility: */
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(IG, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 2, "wrong {IG} channels");
    if (IW != NULL) { float_image_check_size(IW, 1, NX_G, NY_G); }  
    if (verbose) { fprintf(stderr, "input slope map size = %d×%d ...\n", NX_G, NY_G); }
  
    /* Compute output image size: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    if (IZ != NULL) { float_image_check_size(IZ, 1, NX_Z, NY_Z); }
    if (verbose) { fprintf(stderr, "height map size = %d×%d ...\n", NX_Z, NY_Z); }
    float_image_t *OZ = float_image_new(1, (int32_t)NX_Z, (int32_t)NY_Z); 
    (*OZP) = OZ;

    float_image_t *OW = ((IW == NULL) || (OWP == NULL) ? NULL : float_image_new(1, (int32_t)NX_Z, (int32_t)NY_Z));
    if (OWP != NULL) { (*OWP) = OW; }

    /* Build the system: */
    pst_imgsys_t *S = pst_integrate_build_system(IG, IW, verbose);
    if (keepNull)
      { pst_imgsys_fill_holes(S); }
    else
      { pst_imgsys_remove_holes(S); }
    if (verbose) { fprintf(stderr, "system has %d equations and %d unknowns\n", S->N, S->N); }
    if (OW != NULL) { pst_imgsys_extract_system_weight_image(S, OW); }

    /* Get the initial guess {OZ} for the heights: */
    double *VZ = talloc(S->N, double);
    if (IZ != NULL)
      { pst_imgsys_copy_image_to_sol_vec(S, IZ, VZ, 0.0); }
    else
      { for(int32_t uid = 0; uid < S->N; uid++) { VZ[uid] = 0.0; } }
      
    auto void reportSol(uint32_t indent, uint32_t iter, double change, bool_t final, uint32_t N, double Z[]);

    if (verbose) { fprintf(stderr, "solving the system ...\n"); }
    bool_t para = FALSE; /* Should be parameter. 1 means parallel execution, 0 sequential. */
    bool_t szero = TRUE; /* Should be parameter. 1 means adjust sum to zero, 0 let it float. */
    uint32_t *ord = NULL;
    if (topoSort) { ord = pst_imgsys_sort_equations(S); }
    pst_imgsys_solve_iterative
      ( S, VZ, ord, maxIter, convTol, 
        para, szero, verbose, 0, 
        reportStep,
        (reportHeights == NULL ? NULL : &reportSol)
      );
    pst_imgsys_copy_sol_vec_to_image(S, VZ, OZ, 0.0);
        
    free(VZ);
    pst_imgsys_free(S);
    if (ord != NULL) { free(ord); }
    return;

    void reportSol(uint32_t indent, uint32_t iter, double change, bool_t final, uint32_t N, double Z[])
      { assert(reportHeights != NULL);
        pst_imgsys_copy_sol_vec_to_image(S, VZ, OZ, 0.0);
        reportHeights(0, iter, change, final, OZ);
      }
  }
