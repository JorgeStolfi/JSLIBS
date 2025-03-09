/* See pst_integrate_recursive.h */
/* Last edited on 2025-03-03 20:46:01 by stolfi */

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
  ( int32_t level,
    float_image_t *G, 
    float_image_t *H,
    double hintsWeight,
    float_image_t *Z, 
    float_image_t *R, 
    uint32_t maxLevel,
    uint32_t maxIter,
    double convTol,
    bool_t sortSys,
    bool_t verbose,
    pst_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  )
  {
    int32_t indent = (level < -1 ? 0 : 2*level+2); /* Indentation for messages. */
    if (verbose) { fprintf(stderr, "%*s--- entering {%s} level = %d ---\n", indent, "", __FUNCTION__, level); }
    
    /* Check slope map size: */
    int32_t NC_G, NX_G, NY_G;
    float_image_get_size(G, &NC_G, &NX_G, &NY_G);
    if (verbose) { fprintf(stderr, "%*sslope map has %d channels size = %d×%d ...\n", indent, "", NC_G, NX_G, NY_G); }
    demand((NC_G == 2) || (NC_G == 3), "slope map {G} must have 2 or 3 channels");

    if (verbose)
      { fprintf(stderr, "%*sentering level %d with G size %d×%d ...\n", indent, "", level, NX_G, NY_G); }

    /* Check height map size: */
    demand(Z != NULL, "height map must not be null");
    int32_t NX_Z = NX_G + 1;
    int32_t NY_Z = NY_G + 1;
    if (verbose) { fprintf(stderr, "%*sheight map size = %d×%d ...\n", indent, "", NX_Z, NY_Z); }
    float_image_check_size(Z, 2, NX_Z, NY_Z, "bad unknown height map {Z}");

    /* Check hints map size: */
    int32_t NC_H = 0;
    if (H != NULL) 
      { int32_t NX_H, NY_H;
        float_image_get_size(H, &NC_H, &NX_H, &NY_H);
        if (verbose) { fprintf(stderr, "%*shints map {H} has %d channels size = %d×%d ...\n", indent, "", NC_H, NX_H, NY_H); }
        demand((NC_H == 1) || (NC_H == 2), "hints map {H} must have 1 or 2 channels");
        demand((NX_H == NX_Z) && (NY_H == NY_Z), "hints map {H} has wrong size");
      }

    /* Check reference map size: */
    int32_t NC_R = 0;
    if (R != NULL) 
      { int32_t NX_R, NY_R;
        float_image_get_size(R, &NC_R, &NX_R, &NY_R);
        if (verbose) { fprintf(stderr, "%*sreference map {R} has %d channels size = %d×%d ...\n", indent, "", NC_R, NX_R, NY_R); }
        demand((NC_R == 1) || (NC_R == 2), "reference map {R} must have 1 or 2 channels");
        demand((NX_R == NX_Z) && (NY_R == NY_Z), "reference map {R} has wrong size");
      }

    if (reportData != NULL) { reportData(level, G, H, R); }

    /* Decide whether to use multi-scale: */
    bool_t trivial = ((NX_G <= MIN_SLOPE_MAP_SIZE) && (NY_G <= MIN_SLOPE_MAP_SIZE));
    if ((! trivial) && (level < maxLevel))
      { /* Obtain the initial guess {Z} recursively: */
        
        /* Shrink slope map {G} by half: */
        if (verbose) { fprintf(stderr, "%*sshrinking slope and weight maps ...\n", indent, ""); }
        float_image_t *G_red = pst_slope_map_shrink(G, 1.0);
        int32_t NC_G_red, NX_G_red, NY_G_red;
        float_image_get_size(G_red, &NC_G_red, &NX_G_red, &NY_G_red);
        assert(NC_G_red == NC_G);
        
        /* Shrink height map {Z} by one more row and col than {G}: */
        int32_t NX_Z_red = NX_G_red+1;
        int32_t NY_Z_red = NY_G_red+1;
        float_image_t *Z_red = pst_height_map_shrink(Z, 0.5);
        float_image_check_size(Z_red, 2, NX_Z_red, NY_Z_red, "shrunk height map {Z_red} has wrong size");
        
        /* Shrink hints map {H} by half: */
        float_image_t *H_red = (H == NULL ? NULL : pst_height_map_shrink(H, 0.5));
        if (H != NULL) { float_image_check_size(H_red, NC_H, NX_Z_red, NY_Z_red, "shrunk height map {H_red} has wrong size"); }
        
        /* Shrink reference map {R} by half: */
        float_image_t *R_red = (R == NULL ? NULL : pst_height_map_shrink(R, 0.5));
        if (R != NULL) { float_image_check_size(R_red, NC_R, NX_Z_red, NY_Z_red, "shrunk height map {R_red} has wrong size"); }
        
        /* Compute the half-scale height maps: */
        pst_integrate_recursive
          ( level+1, G_red, H_red, hintsWeight, Z_red, R_red,
            maxLevel, 2*maxIter, convTol/2, sortSys,
            verbose, reportData, reportSys, reportStep, reportHeights
          );
        
        /* Expand the computed height map to double size: */
        if (verbose) { fprintf(stderr, "%*sexpanding computed height map from %dx%d to %d×%d ...\n", indent, "", NX_Z_red, NY_Z_red, NX_Z, NY_Z); }
        float_image_t *Z_exp = pst_height_map_expand(Z_red, NX_Z,NY_Z, 2.0);
        float_image_assign(Z, Z_exp);

        /* Free the working storage: */
        float_image_free(G_red);
        float_image_free(Z_red);
        if (H_red != NULL) { float_image_free(H_red); }
        if (R_red != NULL) { float_image_free(R_red); }
      }

    pst_integrate_iterative
      ( G, H, hintsWeight, Z, R, maxIter, convTol, sortSys, verbose, level,
        reportSys, reportStep, reportHeights
      );

    if (verbose) { fprintf(stderr, "%*s--- leaving {%s} level = %d ---\n", indent, "", __FUNCTION__, level); }
    
    return;
  }
