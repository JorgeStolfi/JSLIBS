/* See pst_integrate_recursive.h */
/* Last edited on 2025-01-08 00:35:32 by stolfi */

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
  ( float_image_t *IG, 
    float_image_t *IW, 
    uint32_t level,
    uint32_t maxIter,
    double convTol,
    bool_t topoSort,
    float_image_t **OZP,
    float_image_t **OWP,
    bool_t verbose,
    uint32_t reportStep,
    pst_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    pst_integrate_report_heights_proc_t *reportHeights
  )
  {
    /* Get input image size and check compatibility: */
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(IG, &NC_G, &NX_G, &NY_G);
    demand(NC_G == 2, "wrong {IG} channels");
    if (IW != NULL) { float_image_check_size(IW, 1, NX_G, NY_G); }
    /* Compute output image size: */
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    
    float_image_t *OZ = float_image_new(1, NX_Z, NY_Z);  
    (*OZP) = OZ;

    float_image_t *OW = ((IW == NULL) || (OWP == NULL) ? NULL : float_image_new(1, NX_Z, NY_Z));
    if (OWP != NULL) { (*OWP) = OW; }

    uint32_t indent = 2*level+2; /* Indentation for messages. */

    if (verbose)
      { fprintf(stderr, "%*sEntering level %d with G size %d×%d ...\n", indent, "", level, NX_G, NY_G); }

    if (reportData != NULL) { reportData(level, IG, IW); }

    /* Decide whether to use multi-scale: */
    bool_t trivial = ((NX_G == 1) && (NY_G == 1));
    
    /* Get the initial guess {OZ} for the heights: */
    if (trivial)
      { /* Solution for the 1x1 system: */
        /* Get slopes: */
        float dZdX = float_image_get_sample(IG, 0, 0, 0); 
        float dZdY = float_image_get_sample(IG, 1, 0, 0); 
        /* Assume an affine function that is 0 at center, with the given gradient: */
        double z10 = 0.5*(dZdX - dZdY);
        double z11 = 0.5*(dZdX + dZdY);
        float_image_set_sample(OZ, 0, 0, 0, (float)-z11);
        float_image_set_sample(OZ, 0, 0, 1, (float)-z10);
        float_image_set_sample(OZ, 0, 1, 0, (float)+z10);
        float_image_set_sample(OZ, 0, 1, 1, (float)+z11);
        if (OW != NULL) { float_image_fill_channel(OW, 0, 1.0); }
        /* The initial guess is the solution: */
        if (reportSys != NULL) { reportSys(level, NULL); }
        if (reportHeights != NULL) { reportHeights(level, 0, 0.0, TRUE, OZ); }
      }
    else
      { /* Initialize {OZ} with the solution to a coarser problem. */
      
        /* Shrink slope maps and weight map by half: */
        if (verbose) { fprintf(stderr, "%*sShrinking slope and weight maps ...\n", indent, ""); }
        float_image_t *SIG; 
        float_image_t *SIW;
        pst_slope_and_weight_map_shrink(IG,IW, &SIG,&SIW);
        
        /* Compute the half-scale height map: */
        float_image_t *SOZ, *SOW;
        pst_integrate_recursive
          ( SIG, SIW, level+1, 
            2*maxIter, convTol/2, topoSort,
            &SOZ, &SOW,
            verbose, reportStep,
            reportData, reportSys, reportHeights
          );
        
        /* Expand the computed height map to double size: */
        if (verbose) { fprintf(stderr, "%*sExpanding height map to %d×%d ...\n", indent, "", NX_Z, NY_Z); }
        pst_height_map_expand(SOZ, SOW, OZ, OW);
          
        /* Free the working storage: */
        float_image_free(SIG);
        if (SIW != NULL) { float_image_free(SIW); }
        float_image_free(SOZ);
        if (SOW != NULL) { float_image_free(SOW); }
        
        /* Build the linear system: */
        uint32_t NP_Z = (uint32_t)(NX_Z*NY_Z);
        if (verbose) { fprintf(stderr, "%*sBuilding linear system for %d pixels ...\n", indent, "", NP_Z); }
        bool_t full = FALSE; /* Should be parameter. FALSE means exclude indeterminate pixels. */
        pst_imgsys_t *S = pst_integrate_build_system(IG, IW, verbose);
        if (full)
          { pst_imgsys_fill_holes(S); }
        else
          { pst_imgsys_remove_holes(S); }
        if (OW != NULL) { pst_imgsys_extract_system_weight_image(S, OW); }
        if (verbose) { fprintf(stderr, "%*sSystem has %d equations and %d unknowns\n", indent, "", S->N, S->N); }
        if (reportSys != NULL) { reportSys(level, S); }

        /* Solve the system iteratively for the corner heights: */
        
        auto void reportSol(uint32_t indent, uint32_t iter, double change, bool_t final, uint32_t N, double h[]);
        
        double *VZ = talloc(S->N, double);
        pst_imgsys_copy_image_to_sol_vec(S, OZ, VZ, 0.0);
        if (verbose) { fprintf(stderr, "%*sSolving the system ...\n", indent, ""); }
        bool_t para = FALSE; /* Should be parameter. 1 means parallel execution, 0 sequential. */
        bool_t szero = TRUE; /* Should be parameter. 1 means adjust sum to zero, 0 let it float. */
        uint32_t *ord = NULL;
        if (topoSort) { ord = pst_imgsys_sort_equations(S); }
        pst_imgsys_solve_iterative
          ( S, VZ, ord, maxIter, convTol, 
            para, szero, verbose, level, 
            reportStep,
            (reportHeights == NULL ? NULL : &reportSol)
          );
        pst_imgsys_copy_sol_vec_to_image(S, VZ, OZ, 0.0);

        free(VZ);
        pst_imgsys_free(S);
        if (ord != NULL) { free(ord); }
        
        void reportSol(uint32_t indent, uint32_t iter, double change, bool_t final, uint32_t N, double h[])
          { if (reportHeights != NULL)
              { pst_imgsys_copy_sol_vec_to_image(S, VZ, OZ, 0.0);
                reportHeights(level, iter, change, final, OZ);
              }
          }
      }

   if (verbose)
      { fprintf
          ( stderr, ("%*sLeaving level %d (OZ = %" int64_d_fmt "×%" int64_d_fmt ") ...\n"), 
            indent, "", level, OZ->sz[1], OZ->sz[2]
          );
      }
  }
