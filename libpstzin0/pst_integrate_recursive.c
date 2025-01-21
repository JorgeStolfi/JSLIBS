/* See pst_integrate_recursive.h */
/* Last edited on 2025-01-18 12:32:59 by stolfi */

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

#include <pst_integrate_recursive.h>

void pst_integrate_recursive
  ( float_image_t *G, 
    bool_t keepNull,
    float_image_t *Z, 
    int32_t level,
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    bool_t verbose,
    pst_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  )
  {
    int32_t indent = (level < -1 ? 0 : 2*level+2); /* Indentation for messages. */

    /* Get input image size and check compatibility: */
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 3, "wrong {G} channels");
    /* Compute output image size: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    if (verbose) { fprintf(stderr, "%*sheight map size = %d×%d ...\n", indent, "", NX_Z, NY_Z); }
    demand(Z != NULL, "height map must not be null");
    float_image_check_size(Z, 2, NX_Z, NY_Z, "bad height map");

    if (verbose)
      { fprintf(stderr, "%*sEntering level %d with G size %d×%d ...\n", indent, "", level, NX_G, NY_G); }

    if (reportData != NULL) { reportData(level, G, Z); }

    /* Decide whether to use multi-scale: */
    bool_t trivial = ((NX_G == 1) && (NY_G == 1));
    
    /* Get the initial guess {OZ} for the heights: */
    if (trivial)
      { /* Solution for the 1x1 system: */
        /* Get slopes: */
        float dZdX = float_image_get_sample(G, 0, 0, 0); 
        float dZdY = float_image_get_sample(G, 1, 0, 0); 
        /* Assume an affine function that is 0 at center, with the given gradient: */
        double z10 = 0.5*(dZdX - dZdY);
        double z11 = 0.5*(dZdX + dZdY);
        float_image_set_sample(Z, 0, 0, 0, (float)-z11);
        float_image_set_sample(Z, 0, 0, 1, (float)-z10);
        float_image_set_sample(Z, 0, 1, 0, (float)+z10);
        float_image_set_sample(Z, 0, 1, 1, (float)+z11);
        float_image_fill_channel(Z, 1, 1.0);
        /* The initial guess is the solution: */
        if (reportSys != NULL) { reportSys(level); }
        if (reportHeights != NULL) { reportHeights(level, 0, 0.0, TRUE, Z); }
      }
    else
      { /* Shrink slope maps and weight map by half: */
        if (verbose) { fprintf(stderr, "%*sShrinking slope and weight maps ...\n", indent, ""); }
        float_image_t *RG = pst_slope_map_shrink(G);
        int32_t NC_RG, NX_RG, NY_RG;
        float_image_get_size(RG, &NC_RG, &NX_RG, &NY_RG);
        assert(NC_RG == 2);
        int32_t NX_RZ = NX_RG+1;
        int32_t NY_RZ = NY_RG+1;
        
        /* Compute the half-scale height map: */
        float_image_t *RZ = float_image_new(2, NX_RZ, NY_RZ);
        pst_integrate_recursive
          ( RG, keepNull, RZ,
            level+1, 2*maxIter, convTol/2, topoSort,
            verbose, reportData, reportSys, reportStep, reportHeights
          );
        
        /* Expand the computed height map to double size: */
        if (verbose) { fprintf(stderr, "%*sExpanding height map to %d×%d ...\n", indent, "", NX_Z, NY_Z); }
        Z = pst_height_map_expand(RZ);
          
        /* Free the working storage: */
        float_image_free(RG);
        float_image_free(RZ);
        
        /* Build the linear system: */
        uint32_t NXY_Z = (uint32_t)(NX_Z*NY_Z);
        if (verbose) { fprintf(stderr, "%*sBuilding linear system for %d pixels ...\n", indent, "", NXY_Z); }
        pst_imgsys_t *S = pst_integrate_build_system(G, W, verbose);
        assert(S->N == NXY_Z);
        if (keepNull)
          { pst_integrate_fill_holes(S, NX_Z, NY_Z); }
        else
          { pst_imgsys_remove_holes(S, NULL); }
        if (verbose) { fprintf(stderr, "%*sSystem has %d equations and %d variables\n", indent, "", S->N, S->N); }
        if (reportSys != NULL) { reportSys(level, S, U); }

        if (verbose) { fprintf(stderr, "%*scopying var/eq weights to {U} map ...\n", indent, ""); }
        pst_imgsys_extract_system_eq_tot_weight_image(S, U, 0.0);

        /* Get the initial guess {Z} for the heights: */
        double *h = talloc(S->N, double);
        if (verbose) { fprintf(stderr, "%*staking initial guess of heights ...\n", indent, ""); }
        pst_imgsys_copy_image_to_sol_vec(S, Z, h, 0.0); 

        auto void reportSol(int32_t level, int32_t iter, double change, bool_t final, uint32_t N, double h[]);
        
        if (verbose) { fprintf(stderr, "%*sSolving the system ...\n", indent, ""); }
        bool_t para = FALSE; /* Should be parameter. 1 means parallel execution, 0 sequential. */
        bool_t szero = TRUE; /* Should be parameter. 1 means adjust sum to zero, 0 let it float. */
        uint32_t *ord = NULL;
        if (topoSort) { ord = pst_imgsys_sort_equations(S); }
        pst_imgsys_solve_iterative
          ( S, h, ord, maxIter, convTol, 
            para, szero, verbose, level, 
            reportStep,
            (reportHeights == NULL ? NULL : &reportSol)
          );
        pst_imgsys_copy_sol_vec_to_image(S, h, Z, 0.0);

        free(h);
        pst_imgsys_free(S);
        if (ord != NULL) { free(ord); }
        
        void reportSol(int32_t level, int32_t iter, double change, bool_t final, uint32_t N, double h[])
          { assert(reportHeights != NULL);
            pst_imgsys_copy_sol_vec_to_image(S, h, Z, 0.0);
            reportHeights(level, iter, change, final, Z);
          }
      }

   if (verbose)
      { fprintf
          ( stderr, "%*sleaving level %d (Z = %d×%d) ...\n", 
            indent, "", level, NX_Z, NY_Z
          );
      }
  }
