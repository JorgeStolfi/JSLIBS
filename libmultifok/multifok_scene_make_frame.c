/* See {multifok_scene.h}. */
/* Last edited on 2024-10-28 15:14:10 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <r3.h>
#include <frgb.h>
#include <r2.h>
#include <i2.h>
#include <frgb_ops.h>
#include <jsrandom.h>
#include <wt_table.h>
#include <wt_table_hann.h>

#include <multifok_scene.h>
#include <multifok_sampling.h>
#include <multifok_raytrace.h>
#include <multifok_scene.h>
#include <multifok_scene_tree.h>

#include <multifok_scene_make_frame.h>
  
#define DASHES "------------------------------------------------------------"

multifok_frame_t *multifok_scene_make_frame
  ( int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene, 
    multifok_scene_tree_t *tree,
    multifok_scene_raytrace_pattern_t *pattern,
    double zFoc, 
    double zDep, 
    int32_t HS,
    int32_t KR,
    bool_t verbose,
    multifok_raytrace_debug_pred_t *debug_pix,
    multifok_raytrace_report_ray_proc_t *report_ray
  )
  {
    demand((NX >= 10) && (NY >= 10), "image is too small");
    demand(zDep > 0.0, "invalid {zDep}");
    
    double zMin = scene->dom[2].end[0]; /* Lower bound for {Z} in scene. */
    
    /* Pixel sub-sampling points: */
    int32_t NW = 2*HS+1; /* Pixel subsampling window width and height. */
    int32_t NS = NW*NW;  /* Number of pixel sub-sampling points. */
    r2_t *uSmp;
    double *wSmp;
    multifok_sampling_choose_pixel_sampling_points_and_weights(HS, &NS, &uSmp, &wSmp, verbose);
    
    /* Choose the relative ray tilts and weights: */
    int32_t NR; /* Actual number of aperture rays. */
    r2_t *tRay;
    double *wRay;
    multifok_sampling_choose_ray_tilts_and_weights(KR, NS, &NR, &tRay, &wRay, verbose);
    if (KR == 0)
      { assert(NR == 1); }
    else
      { assert(NR >= NS);
        assert((NR % NS) == 0);
        assert(NR/NS >= KR);
      }
     
    auto void trace_scene(r3_t *pR, r3_t *dR, bool_t debug_loc, r3_t *pHit_P, int32_t NC_loc, float colr[]);
      /* A ray-tracing function suitable for the {trace_ray} argument of
        {multifok_raytrace_make_frame}. It shoots a ray {R} through the
        point {pR} with direction {dR}, which must be downward
        ({dR.c[2] < 0}), both in the scene coordinate system.
        
        The procedure finds the first point {pHit(R)} on the scene's
        surface hit by the ray. It returns that point in {*pHit_P}.
        It also returns in {colr[0..NC-1]} the
        color of the scene's surface at that point. 
        
        If the ray does not hit any object, returns in {*pHit_P} the 
        intersection of the ray with the plane {Z=zMin}, and 
        sets {colr[0..NC-1]} to {0.500}. */
        
    auto void map_point(r2_t *p2_img, r3_t *p3_scene);
      /* An image-to-scene coordinate conversion fucntion, suitable for
        the {map_point} argument of {multifok_raytrace_make_frame}. The point
        {p2_img} is assumed to be in the 2D image coordinate system,
        inside the image domain {[0 _ NX] × [0 _ NY]}. The mapping
        assumes that this rectangle is horizontal in the scene
        coordinate system, fits snugly and centered in the {XY} domain
        of the scene, namely {scene->dom[0] × scene->dom[1]}, and has
        {Z} scene coordinate equal to {zFoc}. The procedure returns the
        corresponding scene coordinates in {*p3_scene}. */

    /* Create the frame images: */
    int32_t NC = 3;
    r3_t dRef = (r3_t){{ 0, 0, -1 }};
    
    multifok_frame_t *frame = multifok_raytrace_make_frame
      ( NC, NX, NY, &trace_scene, &map_point, 
        &dRef, zFoc, zDep, NS, uSmp, wSmp, NR, tRay, wRay, 
        verbose, debug_pix, report_ray
      );
    
    free(tRay);
    free(wRay);

    return frame;
    
    void trace_scene(r3_t *pR, r3_t *dR, bool_t debug_loc, r3_t *pHit_P, int32_t NC_loc, float colr[])
      { 
        assert(NC_loc == NC);
          
        if (debug_loc) { fprintf(stderr, "        %s\n", DASHES); }
        
        /* Ray-trace the scene: */
        multifok_scene_object_t *oHit_ray; /* Object hit, or {NULL} if floor. */
        r3_t pHit_ray; /* Hit point. */
        multifok_scene_raytrace(scene, tree, pR, dR, debug_loc, &oHit_ray, &pHit_ray);
        if (oHit_ray == NULL)
          { /* Ray missed scene entirely, provide a gray plane at {zMin}: */
            for (uint32_t ic = 0;  ic < NC; ic++) { colr[ic] = 0.500; }
            double tHit = (zMin - pR->c[2])/dR->c[2];
            r3_mix(1.0, pR, tHit, dR, &pHit_ray);
          }
        else
          { /* Ray hit object {oHit_ray} at {pHit_ray}: */
            double hHit_ray = pHit_ray.c[2];
            assert(! isnan(hHit_ray));
            interval_t zRange = scene->dom[2];
            /* The scene surface must be entirely contained in the scene's box: */
            if ((hHit_ray < zRange.end[0]) || (hHit_ray > zRange.end[1]))
              { fprintf(stderr, "        ray hit at Z = %16.12f, outside range {%16.12f _ %16.12f]\n", hHit_ray, zRange.end[0], zRange.end[1]);
                assert(FALSE);
              }
            /* Compute the hit point's color: */
            r3_t light_dir = (r3_t){{ 2.0, 3.0, 4.0 }};
            (void)r3_dir(&light_dir, &light_dir);
            frgb_t colr_ray = multifok_scene_raytrace_compute_hit_color(oHit_ray, &(pHit_ray), pattern, &light_dir);
            for (uint32_t ic = 0;  ic < NC; ic++) { colr[ic] = colr_ray.c[ic]; }
          }
        (*pHit_P) = pHit_ray;
        if (debug_loc) { fprintf(stderr, "        %s\n", DASHES); }
      }

    void map_point(r2_t *p2_img, r3_t *p3_scene)
      {
        r3_t p3_img = (r3_t){{ p2_img->c[0], p2_img->c[1], 0.0 }};
        (*p3_scene) = multifok_scene_coords_from_image_coords(&p3_img, NX, NY, scene->dom, zFoc);
      }
  }

#define multifok_scene_make_frame_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"
