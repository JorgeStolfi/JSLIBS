/* See pst_integrate_recursive.h */
/* Last edited on 2025-01-25 09:31:51 by stolfi */

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

#include <pst_integrate_recursive.h>

#define MIN_SLOPE_MAP_SIZE 3
  /* Stop the recursion when the slope map has no more than 
    this number of rows and columns. */

void pst_integrate_recursive
  ( float_image_t *G, 
    float_image_t *H,
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
    float_image_check_size(Z, 2, NX_Z, NY_Z, "bad unknown height map {Z}");
    if (H != NULL) { float_image_check_size(H, 2, NX_Z, NY_Z, "bad external height map {H}"); }

    if (verbose)
      { fprintf(stderr, "%*sentering level %d with G size %d×%d ...\n", indent, "", level, NX_G, NY_G); }

    /* Decide whether to use multi-scale: */
    bool_t trivial = ((NX_G <= MIN_SLOPE_MAP_SIZE) && (NY_G <= MIN_SLOPE_MAP_SIZE));
    if (! trivial)
      { /* Obtain the initial guess {Z} recursively: */
        /* Shrink all maps by half: */
        if (verbose) { fprintf(stderr, "%*sshrinking slope and weight maps ...\n", indent, ""); }
        float_image_t *RG = pst_slope_map_shrink(G);
        int32_t NC_RG, NX_RG, NY_RG;
        float_image_get_size(RG, &NC_RG, &NX_RG, &NY_RG);
        assert(NC_RG == 2);
        
        int32_t NX_RZ = NX_RG+1;
        int32_t NY_RZ = NY_RG+1;
        float_image_t *RZ = pst_height_map_shrink(Z);
        float_image_check_size(RZ, 2, NX_RZ, NY_RZ, "shrunk height map {RZ} has wrong size");
        float_image_t *RH = (H == NULL ? NULL : pst_height_map_shrink(H));
        if (H != NULL) { float_image_check_size(RH,2, NX_RZ, NY_RZ, "shrunk height map {RH} has wrong size"); }
        
        /* Compute the half-scale height maps: */
        pst_integrate_recursive
          ( RG, RH, RZ,
            level+1, 2*maxIter, convTol/2, topoSort,
            verbose, reportData, reportSys, reportStep, reportHeights
          );
        
        /* Expand the computed height map to double size: */
        if (verbose) { fprintf(stderr, "%*sexpanding height map to %d×%d ...\n", indent, "", NX_Z, NY_Z); }
        Z = pst_height_map_expand(RZ, NX_Z,NY_Z);

        /* Free the working storage: */
        float_image_free(RG);
        if (RH != NULL) { float_image_free(RH); }
        float_image_free(RZ);
      }

    if (reportData != NULL) { reportData(level, G, H, Z); }

    pst_integrate_iterative
      ( G , H, Z, maxIter, convTol, topoSort, verbose, level,
        reportSys, reportStep, reportHeights
      );

   if (verbose)
      { fprintf
          ( stderr, "%*sleaving level %d (Z = %d×%d) ...\n", 
            indent, "", level, NX_Z, NY_Z
          );
      }
    
    return;
  }
