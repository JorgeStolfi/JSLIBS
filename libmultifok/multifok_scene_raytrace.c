/* See {multifok_scene_raytrace.h}. */
/* Last edited on 2025-02-06 21:22:38 by stolfi */

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
#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsqroots.h>

#include <multifok_scene.h>
#include <multifok_scene_raytrace.h>
#include <multifok_scene_object_raytrace.h>
#include <multifok_scene_object_normal.h>
#include <multifok_scene_tree.h>

#define FUDGE (1.0e-6)
  /* Fudge amount to expand bounding boxes to account for roundoff. */

double multifok_scene_raytrace_tree
  ( multifok_scene_tree_t *tr, 
    r3_t *p, 
    r3_t *d, 
    double tMin,
    double tMax,
    interval_t ray_bbox[],
    bool_t verbose, 
    int32_t level,
    multifok_scene_object_t **oHit_P
  );
  /* Same as {multifok_scene_raytrace} but for a subtree {tr} of the whole
    scene tree, at depth {level}, with the ray clipped to the parameter
    range {[tMin _ tMax]}. */
 
void multifok_scene_raytrace
  ( multifok_scene_t *scene,
    multifok_scene_tree_t *tree,
    r3_t *p, 
    r3_t *d, 
    bool_t verbose,
    multifok_scene_object_t **oHit_P, 
    r3_t *pHit_P
  )
  {
    if (verbose || isnan(p->c[2]) || isnan(d->c[2])) 
      { r3_gen_print(stderr, p, "%12.8f", "        tracing ray p = ( ", " ", " )");
        r3_gen_print(stderr, d, "%12.8f", "  d = ( ", " ", " )\n");
      }
      
    demand(d->c[2] < 0, "invalid direction vector");
    demand(p->c[2] >= 0, "invalid focus plane {Z}");
    
    double zMin = scene->dom[2].end[0];
    double zMax = scene->dom[2].end[1];
   
    /* Assumes that the {Z} of the scene's visible surface is in {[zMin _ zMax]}: */
    double tMin = (zMax + 2*FUDGE - p->c[2])/d->c[2];  /* Note: {tMin} is at {zMax}. */
    double tMax = (zMin + 2*FUDGE - p->c[2])/d->c[2];  /* Note: {tMax} is at {zMin}. */
    assert(tMin < tMax);
    int32_t level = 0;
    interval_t ray_bbox[3];
    multifok_scene_raytrace_get_ray_bbox(p, d, tMin, tMax, ray_bbox);
    double tHit = multifok_scene_raytrace_tree(tree, p, d, tMin, tMax, ray_bbox, verbose, level, oHit_P);
    if (tHit != +INF)
      { assert((tHit >= tMin) && (tHit <= tMax));
        assert((*oHit_P) != NULL);
        double zHit = p->c[2] + tHit*d->c[2];
        double xHit = p->c[0] + tHit*d->c[0];
        double yHit = p->c[1] + tHit*d->c[1];
        (*pHit_P) = (r3_t) {{ xHit, yHit, zHit }};
      }
    else
      { (*pHit_P) = (r3_t) {{ NAN, NAN, NAN }}; }
  }

double multifok_scene_raytrace_tree
  ( multifok_scene_tree_t *tr, 
    r3_t *p, 
    r3_t *d, 
    double tMin,
    double tMax,
    interval_t ray_bbox[],
    bool_t verbose, 
    int32_t level,
    multifok_scene_object_t **oHit_P
  )
  { 
    multifok_scene_object_t *oHit = NULL;
    double tHit = +INF;
    
    if (tr != NULL)
      { if (verbose) { fprintf(stderr, "        %*stracing tree with ray tMin = %12.8f tMax = %12.8f\n", 2*level, "", tMin, tMax); }
        assert(tMin <= tMax);

        /* Check if it it intersects the tree's bbox: */
        bool_t ok = TRUE; /* Set to false if ray's bbox is disjoint from tree's bbox. */
        for (uint32_t j = 0;  (j < 3) & ok; j++) 
          { if (ray_bbox[j].end[1] < tr->bbox[j].end[0]) { ok = FALSE; }
            if (ray_bbox[j].end[0] > tr->bbox[j].end[1]) { ok = FALSE; }
          }
        if (! ok)
          { if (verbose){ fprintf(stderr, "        %*smissed bounding box\n", 2*level, ""); } }
        else
          { /* Bounding boxes intersect, we have a chance. */
            /* Ray-trace the root object: */
            multifok_scene_object_t *oRoot = tr->obj;
            if (verbose) 
              { char *typeX = multifok_scene_object_type_to_string(oRoot->type);
                fprintf(stderr, "        %*strying root obj %d type %s\n", 2*level, "", oRoot->ID, typeX);
              }
            double tHit_root = multifok_scene_object_raytrace(oRoot, p, d, tMin, tMax, verbose);
            if (tHit_root != +INF)
              { /* Root object hit replaces previous hit: */
                assert((tHit_root >= tMin) && (tHit_root <= tMax));
                if (verbose) { fprintf(stderr, "        %*shit root obj %d at t = %12.8f tMax = %12.8f\n", 2*level, "", oRoot->ID, tHit_root, tMax); }
                /* Save object, bring down {tMax}: */
                oHit = oRoot;
                tHit = tHit_root;
                tMax = tHit_root;
                /* Update the ray's bounding box {ray_bbox[0..2]}: */
                multifok_scene_raytrace_get_ray_bbox(p, d, tMin, tMax, ray_bbox);                
              }
            else
              { if (verbose) { fprintf(stderr, "        %*smissed/rejected root obj %d\n", 2*level, "", oRoot->ID); } 
              }

            /* Ray trace the subtrees. We start with the subtree
              {sub[ic]} that seems closest to the ray in direction
              {axis}. Then we ray-trace {sub[ic]} and {sub[1-ic]},
              in that order, remembering the highest hit as we
              go. */
            uint8_t axis = tr->axis;
            int32_t ic; /* First subtree to try. */
            if ((tr->sub[0] != NULL) && (tr->sub[1] != NULL))
              { double obc = interval_mid(&(tr->obj->bbox[axis]));
                double chd[2]; /* Distance to coord {axis} of child box centers. */
                for (uint32_t jc = 0;  jc < 2; jc++) 
                  { chd[jc] = fabs(obc - interval_mid(&(tr->sub[jc]->bbox[axis]))); }
                ic = (chd[0] < chd[1] ? 0 : 1);
              }
            else
              { /* No subtrees, {ic} is irrelevant. */
                ic = 0;
              }

            for (uint32_t kc = 0;  kc < 2; kc++) 
              { /* Ray-trace {sub[ic]} */
                if (tr->sub[ic] != NULL)
                  { if (verbose) { fprintf(stderr, "        %*strying subtree %d\n", 2*level, "", ic); }
                    multifok_scene_object_t *oHit_sub;
                    double tHit_sub = multifok_scene_raytrace_tree(tr->sub[ic], p, d, tMin, tMax, ray_bbox, verbose, level+1, &oHit_sub);
                    if (tHit_sub != +INF)
                      { /* We got a hit: */
                        assert((tHit_sub >= tMin) && (tHit_sub <= tMax));
                        assert(oHit_sub != NULL);
                        if (verbose) 
                          { fprintf(stderr, "        %*sgot hit with subtree %d obj %d at t = %12.8f tMax = %12.8f\n", 2*level, "", ic, oHit_sub->ID, tHit_sub, tMax); }
                        oHit = oHit_sub;
                        tHit = tHit_sub;
                        tMax = tHit_sub;
                        /* Update the ray's bounding box {ray_bbox[0..2]}: */
                        multifok_scene_raytrace_get_ray_bbox(p, d, tMin, tMax, ray_bbox);
                      }
                    else
                      { if (verbose) { fprintf(stderr, "        %*smissed hit with subtree %d\n", 2*level, "", ic); } 
                      }
                  }
                /* Try the other subtree: */
                ic = 1 - ic;
              }
          }
      }
    /* Return to caller: */
    if (verbose) 
      { if (oHit != NULL) 
          { fprintf(stderr, "        %*sreturning hit with obj %d at t = %12.8f\n", 2*level, "", oHit->ID, tHit);
          }
        else
          { fprintf(stderr, "        %*sreturning with no hit\n", 2*level, ""); }
      }
    (*oHit_P) = oHit;
    return tHit;
  }
        
void multifok_scene_raytrace_get_ray_bbox
  ( r3_t *p,
    r3_t *d,
    double tMin,
    double tMax,
    interval_t bbox[]
  )
  { 
    for (uint32_t j = 0;  j < 3; j++) 
      { double cLo = p->c[j] + d->c[j]*tMin; /* Coordinate at bottom of ray. */
        double cHi = p->c[j] + d->c[j]*tMax; /* Coordinate at top of ray. */
        bbox[j] = (interval_t){{ fmin(cLo, cHi) - 0.5*FUDGE, fmax(cLo, cHi) + 0.5*FUDGE }};
      }
  }

/* SCENE TO IMAGE COORDINATE MAPPING */

r3_t multifok_scene_coords_to_image_coords
  ( r3_t *p_scene,
    double zFoc,
    interval_t scene_dom[],
    int32_t NX,
    int32_t NY
  )
  {
    double pixSize = multifok_scene_pixel_size(NX, NY, scene_dom);
    
    double x_scene_ctr = interval_mid(&(scene_dom[0]));
    double y_scene_ctr = interval_mid(&(scene_dom[1]));
    double x_img_ctr = 0.5*NX;
    double y_img_ctr = 0.5*NY;

    double x_img = y_img_ctr + (p_scene->c[0] - x_scene_ctr)/pixSize;
    double y_img = x_img_ctr + (p_scene->c[1] - y_scene_ctr)/pixSize;
    double z_img = (p_scene->c[2] - zFoc)/pixSize;
    return (r3_t){{ x_img, y_img, z_img }};    
  }
  
r3_t multifok_scene_coords_from_image_coords
  ( r3_t *p3_img,
    int32_t NX,
    int32_t NY,
    interval_t scene_dom[],
    double zFoc
  )
  {
    double pixSize = multifok_scene_pixel_size(NX, NY, scene_dom);

    double x_scene_ctr = interval_mid(&(scene_dom[0]));
    double y_scene_ctr = interval_mid(&(scene_dom[1]));
    double x_img_ctr = 0.5*NX;
    double y_img_ctr = 0.5*NY;

    double x_scene = x_scene_ctr + (p3_img->c[0] - x_img_ctr)*pixSize;
    double y_scene = y_scene_ctr + (p3_img->c[1] - y_img_ctr)*pixSize;
    double z_scene = zFoc + p3_img->c[2]*pixSize;
    return (r3_t){{ x_scene, y_scene, z_scene }};
  }
    
double multifok_scene_pixel_size(int32_t NX, int32_t NY, interval_t scene_dom[])
  { demand((NX > 0) && (NY > 0), "invalid {NX,NY}");
    double WX = 2*interval_rad(&(scene_dom[0]));
    double WY = 2*interval_rad(&(scene_dom[1]));
    demand((WX > 0) && (WY > 0), "invalid scene domain");
    double pixSize = fmin(WX/((double)NX), WY/((double)NY));
    return pixSize;
  }


#define multifok_scene_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

