/* See {multifok_scene.h}. */
/* Last edited on 2024-12-06 05:56:27 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <interval_io.h>
#include <r3.h>
#include <frgb.h>
#include <r2.h>
#include <jsrandom.h>

#include <multifok_scene.h>
#include <multifok_scene_object.h>

uint32_t multifok_scene_choose_object_count(bool_t floorOnly, interval_t dom[], double rMin, double rMax, double minSep);
  /* Chooses the ideal number of objects (disks or balls) that {multifok_scene_throw_busy} 
    should try to put in a scene.  The parameter {minSep} has the meaning described 
    under {multifok_scene_throw_busy}. */

/* IMPLEMENTATIONS */

multifok_scene_t *multifok_scene_new(interval_t dom[], bool_t verbose)
  {
    /* Allocate record: */
    multifok_scene_t *scene = talloc(1, multifok_scene_t);
    
    /* Store data in scene: */
    for (uint32_t j = 0;  j < 3; j++) { scene->dom[j] = dom[j]; }
    if (verbose) { fprintf(stderr, "creating empty scene,\n"); }
    scene->NO = 0;
    scene->objs = NULL;
    return scene;
  }

void multifok_scene_throw_objects
  ( multifok_scene_t *scene,
    bool_t floorOnly,
    bool_t flatFloor,
    double rMin, 
    double rMax,
    double minSep, 
    bool_t verbose
  )
  {
    bool_t debug = FALSE;
    
    demand((scene->NO == 0) && (scene->objs == NULL), "scene already has objects");
    
    /* Determine the max number of objects: */
    uint32_t NO_max = multifok_scene_choose_object_count(floorOnly, scene->dom, rMin, rMax, minSep);
    if (verbose) { fprintf(stderr, "trying to generate %d objects\n", NO_max); }

    uint32_t NO = 0; /* Number of objects already generated. */
    multifok_scene_object_t *objs = talloc(NO_max, multifok_scene_object_t);

    interval_t *dom = scene->dom;

    /* Generate the objects {objs[0..NO-1]}. */
    /* If overlapping, every try is valid, otherwise we need many more tries: */
    uint32_t NT = (minSep >= 0 ? 50 : 1)*NO_max; /* Number of tries. */
    
    for (uint32_t kt = 0;  (kt < NT) && (NO < NO_max); kt++)
      { multifok_scene_object_t obj;
        if (NO == 0)
          { /* First object is a background flat or ramp floor: */
            obj = multifok_scene_object_background_make(dom, flatFloor);
          }
        else
          { /* Subsequent objects are foreground ones: */
            assert(! floorOnly);
            /* Generate a random object {obj}: */
            obj = multifok_scene_object_foreground_throw(dom, minSep, rMin, rMax, debug);
          }
        if (debug) { multifok_scene_object_print(stderr, "  ", &obj, ""); }
        
        /* Check containment in {Z}: */
        assert(dom[2].end[0] <= obj.bbox[2].end[0]);
        assert(obj.bbox[2].end[0] <= obj.bbox[2].end[1]);
        assert(obj.bbox[2].end[1] <= dom[2].end[1]);
        
        int32_t ko_overlap = -1; /* Index of object overlapped by {obj}, or {-1}. */
        if (minSep >= 0 && NO > 0)
          { demand(flatFloor, "inconsistent parameters");
            /* Reject {obj} if it overlaps previous foreground objects in {X} and {Y}: */
            /* Check for {XY} overlaps: */
            for (int32_t ko = 0;  (ko < NO) && (ko_overlap < 0); ko++)
              { multifok_scene_object_t *objk = &(objs[ko]);
                bool_t overlap = multifok_scene_object_XY_overlap(objk, &obj, minSep);
                if (overlap) { ko_overlap = ko; }
              }
            assert(multifok_scene_object_XY_is_inside(&obj, scene->dom, minSep));
          }
        if (ko_overlap < 0)
          { if (debug) { fprintf(stderr, " (accepted, ID = %d)\n", NO); }
            obj.ID = (multifok_scene_object_ID_t)NO;
            if (verbose && (! debug)) {  multifok_scene_object_print(stderr, "  ", &obj, "\n"); }
            objs[NO] = obj;
            NO++;
          }
        else
          { if (debug) { fprintf(stderr, " (overlaps %d, rejected)\n", ko_overlap); } }
      }
    
    assert(NO <= NO_max);
    if (NO < NO_max)
      { if (verbose) { fprintf(stderr, "generated only %d objects\n", NO); }
        /* Trim array: */
        objs = retalloc(objs, NO, multifok_scene_object_t);
      }
     
    /* Store data in scene: */
    scene->NO = NO;
    scene->objs = objs;
    return;
  }

uint32_t multifok_scene_choose_object_count
  ( bool_t floorOnly,
    interval_t dom[],
    double rMin,
    double rMax,
    double minSep
  )
  {
    if (floorOnly) { return 1; }
    
    /* Compute the average area {aObj} of an object with radius in {[rMin _ rMax]}, accounting for min sep: */
    double rfat = (minSep >= 0 ? 0.5*minSep : 0);
    double rr0 = rMin + rfat; /* Min radius object occup. */
    double rr1 = rMax + rfat; /* Max radius object occup. */
    double aObj = M_PI*(rr1*rr1*rr1 - rr0*rr0*rr0)/(rr1 - rr0)/3; /* Average disk area. */
    double rObj = sqrt(aObj/M_PI); /* RMS average object radius. */
    
    /* Compute the {XY} area {aBox} available for those disks: */
    double wd[2];
    for (uint32_t j = 0;  j < 2; j++)
      { wd[j] = dom[j].end[1] - dom[j].end[0]; 
        if (minSep >= 0)
          { wd[j] -= 2*minSep; }
        else
          { wd[j] += 2*rObj; }
        demand (wd[j] >= 0, "dom too tight");
      }
    double aBox = wd[0]*wd[1]; /* Total useful area of dom. */

    /* Decide max number of foreground objects {NO_max}. */
    /* If non-overlapping, limited by area ratio, else more than that: */
    uint32_t NO_max = (minSep >= 0 ? 1 : 2)*(uint32_t)ceil(aBox/aObj);
    /* Ensure at least one foregreound object: */
    if (NO_max <= 0) { NO_max = 1; }
    
    return NO_max + 1; /* Including the background object. */
  }
      
r3_t multifok_scene_box_center(interval_t box[])
  { r3_t ctr;
    for (uint32_t j = 0;  j < 3; j++) { ctr.c[j] = interval_mid(&(box[j])); }
    return ctr;
  }
  
r3_t multifok_scene_box_radius(interval_t box[])
  { r3_t rad;
    for (uint32_t j = 0;  j < 3; j++) { rad.c[j] = interval_rad(&(box[j])); }
    return rad;
  }

void multifok_scene_print_box(FILE *wr, char *pref, interval_t box[], char *suff)
  { if (pref != NULL) { fputs(pref, wr); }
    for (uint32_t j = 0;  j < 3; j++)
      { if (j > 0) { fputs(" × ", wr); }
        interval_gen_print(wr, &(box[j]), "%12.6f", "[ ", " _ ", " ]");
      }
    if (suff != NULL) { fputs(suff, wr); }
  }
    

void multifok_scene_check_object_IDs(uint32_t NO, multifok_scene_object_t objs[])
  { fprintf(stderr, "checking object IDs...\n");
    bool_t seen[NO];
    for (uint32_t ko = 0;  ko < NO; ko++) { seen[ko] = FALSE; }
    for (uint32_t ko = 0;  ko < NO; ko++) 
      { multifok_scene_object_ID_t ID = objs[ko].ID; 
        assert((0 <= ID) && (ID < NO));
        assert(! seen[ID]);
        seen[ID] = TRUE;
      }
  }

#define multifok_scene_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

