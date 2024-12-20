/* See {multifok_scene_tree.h}. */
/* Last edited on 2024-12-06 07:09:46 by stolfi */

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

#include <multifok_scene.h>
#include <multifok_scene_object.h>
#include <multifok_scene_tree.h>

void multifok_scene_tree_sort_objects(uint32_t NO, multifok_scene_object_t objs[], uint8_t axis);
  /* Rearranges {objs[0..NO-1]} in increasing order of coordinate {axis} (0 or 1). */
    
void multifok_scene_tree_check_object_order(uint32_t NO, multifok_scene_object_t objs[], uint8_t axis);
  /* Checks whether the centers of objects {objs[0..NO-1]} are sorted by increasing {axis}
    coordinate. */

multifok_scene_tree_t *multifok_scene_tree_build
  ( uint32_t NO,
    multifok_scene_object_t objs[],
    int32_t debug_level
  )
  {
    if (debug_level > 0) { fprintf(stderr, "  %*s> %s NO = %d\n", 2*debug_level, "", __FUNCTION__, NO); }
    multifok_scene_tree_t *tr = NULL;
    if (NO != 0)
      { tr = talloc(1, multifok_scene_tree_t);
        tr->sub[0] = NULL;
        tr->sub[1] = NULL;
        /* Choose {axis}, sort objects: */
        if (NO == 1)
          { tr->axis = UINT8_MAX; }
        else
          { /* We must choose {axis}, sort the objects, build {sub[0],sub[1]}. */
         
            /* Find the largest bounding dimension, {X} or {Y}, and split perp to that axis: */
            interval_t bboxXY[2]; /* Bounding box of all objects in {X} and {Y} only. */
            double wd[2]; /* Width of {bboxXY[j]}. */
            for (uint32_t jxy = 0;  jxy < 2; jxy++)
              { bboxXY[jxy] = (interval_t){{ +INF, -INF }}; /* Empty interval. */
                for (uint32_t ko = 0;  ko < NO; ko++)
                  { bboxXY[jxy] = interval_join(&(bboxXY[jxy]), &(objs[ko].bbox[jxy])); }
                wd[jxy] = interval_width(&(bboxXY[jxy]));
                assert(wd[jxy] > 0);
              }
            tr->axis = (wd[0] > wd[1] ? 0 : 1);
            
            /* Sort objects along that axis: */
            multifok_scene_tree_sort_objects(NO, objs, tr->axis);
          }
        /* Pick root object {objs[mo]}: */
        uint32_t mo = NO/2; /* Index of middle object */
        tr->obj = &(objs[mo]);
        
        /* Set {tr->bbox[0..2]} to the bounding box of the root object: */
        for (uint32_t j = 0;  j < 3; j++) { tr->bbox[j] = tr->obj->bbox[j]; }
            
        /* Recurse on the two halves: */
        for (uint32_t ic = 0;  ic < 2; ic++)
          { /* Define the objects {objs[ko_min..ko_max]} that go into {sub[ic]}: */
            uint32_t ko_min = (ic == 0 ? 0 : mo + 1);
            uint32_t ko_max = (uint32_t)(ic == 0 ? mo-1 : NO-1);
            uint32_t NO_sub = (uint32_t)(ko_max + 1 - ko_min);
            assert((NO_sub >= 0) && (NO_sub <= NO));
            if (NO_sub == 0)
              { tr->sub[ic] = NULL; }
            else
              { /* Build {sub[ic]}: */
                int32_t debug_sub = (debug_level < 0 ? -1 : debug_level+1);
                multifok_scene_tree_t *sub = multifok_scene_tree_build(NO_sub, &(objs[ko_min]), debug_sub);
                assert(sub != NULL);
                /* Merge {sub[ic].bbox} into {bbox}: */
                for (uint32_t j = 0;  j < 3; j++) { tr->bbox[j] = interval_join(&(tr->bbox[j]), &(sub->bbox[j])); }
                tr->sub[ic] = sub;
              }
          }
        /* Check that the objects are still OK: */
        /* multifok_scene_check_object_IDs(NO, objs, -1); */
      }
    if (debug_level > 0) { fprintf(stderr, "  %*s< %s\n", 2*debug_level, "", __FUNCTION__); }
    return tr;
  }
            
void multifok_scene_tree_sort_objects
  ( uint32_t NO,
    multifok_scene_object_t objs[],
    uint8_t axis
  )  
  { demand((axis == 0) || (axis == 1), "invalid {axis}");
    auto int32_t objcomp(const void*a, const void*b);
      /* Assumes {*a,*b} are {}. Compares their centers along {axis}, returns {-1,0,+1}. */
            
    /* multifok_scene_check_object_IDs(NO, objs); */
    if (NO > 1) { qsort((void*)objs, NO, sizeof(multifok_scene_object_t), objcomp); }
    /* multifok_scene_check_object_IDs(NO, objs); */
    /* multifok_scene_check_object_order(NO, objs, axis); */

    return;
                
    int32_t objcomp(const void*a, const void*b)
      { multifok_scene_object_t *obja = (multifok_scene_object_t *)a;
        multifok_scene_object_t *objb = (multifok_scene_object_t *)b;
        double ca = interval_mid(&(obja->bbox[axis]));
        double cb = interval_mid(&(objb->bbox[axis]));
        if (ca < cb)
          { return -1; }
        else if (ca > cb)
          { return +1; }
        else
          { return 0; }
      }
  }

void multifok_scene_tree_check_object_order(uint32_t NO, multifok_scene_object_t objs[], uint8_t axis)
  { demand((axis == 0) || (axis == 1), "invalid {axis}");
    fprintf(stderr, "checking center ordering by axis %d...\n", axis);
    double ctr_prev = -INF;
    for (uint32_t ko = 0;  ko < NO; ko++) 
      { double ctr = interval_mid(&(objs[ko].bbox[axis]));
        assert(ctr >= ctr_prev);
        ctr_prev = ctr;
      }
  }

#define multifok_scene_tree_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

