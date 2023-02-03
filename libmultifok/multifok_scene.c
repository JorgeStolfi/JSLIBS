/* See {multifok_scene.h}. */
/* Last edited on 2023-01-31 19:17:10 by stolfi */

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
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsqroots.h>

#include <multifok_scene.h>

#define ZMAX multifok_scene_ZMAX
  /* Shorter name. */
  
#define FUDGE (1.0e-6)
  /* Fudge amount to expand bounding boxes to account for roundoff. */
  
#define DASHES "------------------------------------------------------------"

int32_t multifok_scene_choose_object_count(interval_t dom[], double rMin, double rMax, double minSep);
  /* Chooses the ideal number of objects (disks or balls) that {multifok_scene_throw_busy} 
    should try to put in a scene.  The parameter {minSep} has the meaning described 
    under {multifok_scene_throw_busy}. */

void multifok_scene_throw_colors(frgb_t *tone, frgb_t *bg, frgb_t *fg);
  /* Generates two random contrasting colors, with general hue depending on {warm}. */

typedef struct multifok_scene_tree_t
  { multifok_scene_object_t *obj; 
    int8_t axis;                            /* Axis perpendicular to split plane, 0 or 1. */
    interval_t bbox[3];                     /* Bounding box of all objects. */
    struct multifok_scene_tree_t *child[2]; /* Children nodes, or NULL. */
  } multifok_scene_tree_t;
  /* Node in an acceleration tree. 
  
    The node stores a pointer to one object {obj}, and to two subtrees
    {child[0,1]}. The two sets of objects are disjoint and {obj} is in
    neither of them. The {bbox}s of the children may overlap, and
    overlap with {obj}; but the objects in {child[0]} tend to have lower
    coordinates along the specified {axis} than those in {child[1]}.
    
    If both children are {NULL} then {axis} is irrelevant. */

void multifok_scene_sort_objects(int32_t NO, multifok_scene_object_t objs[], int8_t axis);
  /* Rearranges {objs[0..NO-1]} in increasing order of coordinate {axis} (0 or 1). */

void multifok_scene_check_objs(int32_t NO, multifok_scene_object_t objs[], int8_t axis);
  /* Checks whether the {ID} fields of {objs[0..NO-1]} are OK. If {axis} is not 
    negative, also checks that the object centers are sorted by increasing {axis}
    coordinate. */

/* RAY-TRACING SCENES */
    
multifok_scene_tree_t *multifok_scene_build_tree(int32_t NO, multifok_scene_object_t objs[]);
  /* Builds an acceleration tree for the objects {objs[0..NO-1]}.  Will rearrange the
    objects in the process. */
    
void multifok_scene_ray_trace_tree
  ( multifok_scene_tree_t *tr, 
    r3_t *p, 
    r3_t *d, 
    double zMin,
    bool_t debug, 
    int32_t level,
    multifok_scene_object_t **hob_P, 
    r3_t *hpt_P
  );
  /* Finds the first intersection of the ray defined by {p} and {d} and the objects in the tree {tr}.
    Only considers the part of the ray between {Z=zMin} and {Z=ZMAX}.
    
    Assumes that the ray is mostly vertical ({d.c[2] >> 0}) and directed towards {-d}. If the ray hits at least one object,
    returns in {*hpt_P} the hit point with largest {Z}, and in {*hob_P} the corresponding object.
    If the ray misses all objects, returns {NULL} in {*hob_P}, and {*hpt_P} will be undefined.
    
    The {level} is used for indenting printouts. */

bool_t multifok_scene_ray_trace_object(multifok_scene_object_t *obj, r3_t *p, r3_t *d, double zMin, bool_t debug, r3_t *hpt_P);
  /* Traces one ray {R} that goes through the point {p} with the direction parallel to {d}, assumed to be 
    not horizontal, and limited in {Z} to the range {zMin,ZMAX}.
    If the ray hits the object {obj}, returns {TRUE} and stores the hit point in {*hpt_P}.
    Otherwise returns {FALSE} and leaves {*hpt_P} unchanged. */
    
void multifok_scene_ray_trace_floor(interval_t dom[], bool_t flatFloor, r3_t *p, r3_t *d, bool_t debug, r3_t *hpt_P);   
  /* Traces one ray {R} that goes through the point {p} with the direction parallel to {d}, assumed to be 
    not horizontal.  Computes the intersection of the ray with the floor and returns that point in {*hpt_P}.
    The floor is defined by the scene's domain {scene.dom[0..2]} and the flag {scene.flatFloor}, as
    described under {multifok_scene_t}. */
    
void multifok_scene_get_ray_bbox(r3_t *p, r3_t *d, double zMin, interval_t bbox[]);
  /* Sets {bbox[0..2]} to the bounding box of the part of the ray between {Z=zMin} and {Z=ZMAX}. */
  
/* IMPLEMENTATIONS */

multifok_scene_t *multifok_scene_new(interval_t dom[], bool_t flatFloor, bool_t verbose)
  {
    demand((dom[2].end[0] >= 0.0) && (dom[2].end[1] <= ZMAX), "bad {dom[2]}");
     
    /* Allocate record: */
    multifok_scene_t *scene = (multifok_scene_t *)notnull(malloc(sizeof(multifok_scene_t)), "no mem");
    
    /* Store data in scene: */
    for (int32_t j = 0; j < 3; j++) { scene->dom[j] = dom[j]; }
    if (verbose) { fprintf(stderr, "creating empty scene, {flatFloor = %c}\n", "FT"[flatFloor]); }
    scene->flatFloor = flatFloor;
    scene->NO = 0;
    scene->objs = NULL;
    scene->tree = NULL;
    return scene;
  }

void multifok_scene_throw_objects
  ( multifok_scene_t *scene,
    double rMin, 
    double rMax,
    double minSep, 
    bool_t verbose
  )
  {
    demand((scene->NO == 0) && (scene->objs == NULL), "scene already has objects");

    int32_t NR = 128; /* Size of {random} state array in bytes */
    char newstate[128];
    char *oldstate = initstate(4615*417, newstate, NR);

    interval_t dom[3]; /* Box for throwing the object in. */
    for (int32_t j = 0; j < 3; j++) { dom[j] = scene->dom[j]; }
    
    /* Determine the max number of objects: */
    int32_t NO_max = multifok_scene_choose_object_count(dom, rMin, rMax, minSep);
    if (verbose) { fprintf(stderr, "trying to generate %d objects\n", NO_max); }
    multifok_scene_object_t *objs = (multifok_scene_object_t*)notnull(malloc(NO_max*sizeof(multifok_scene_object_t)), "no mem");

    /* Generate the objects {objs[0..NO-1]} If overlapping, every try is valid, otherwise we need many more tries: */
    int32_t NT = (minSep >= 0 ? 50 : 1)*NO_max; /* Number of tries. */
    int32_t NO = 0; /* Number of objects actually generated. */
    for (int32_t kt = 0; (kt < NT) && (NO < NO_max); kt++)
      { /* Generate a random object {obj}: */
        multifok_scene_object_t obj = multifok_scene_object_throw(NO, dom, minSep, rMin, rMax);
        bool_t ok = TRUE; /* False if object overlaps the previous ones in {XY}. */
        if (minSep >= 0)
          { /* Check for {XY} overlaps: */
            for (int32_t jd = 0; (jd < NO) && ok; jd++)
              { multifok_scene_object_t *objj = &(objs[jd]);
                /* Check for overlap between {objj} and {obj} plus min separation: */
                double dxy = hypot(objj->ctr.c[0] - obj.ctr.c[0], objj->ctr.c[1] - obj.ctr.c[1]);
                if (dxy < obj.rad + minSep + objj->rad) { ok = FALSE; }
              }
          }
        if (ok) 
          { objs[NO] = obj;
            if (verbose) 
              { r3_gen_print(stderr, &(obj.ctr), "%12.8f", "  ctr = ( ", " ", " )");
                fprintf(stderr, " rad = %12.8f", obj.rad);
                fprintf(stderr, " %s ID = %d\n", (obj.flat ? "disk" : "ball"), obj.ID);
              }
            /* Sanity check: object's containment relative to {scene->dom}: */
            for (int32_t j = 0; j < 3; j++)
              { /* Determine bounds {cmin,cmax} of object along axis {j}: */
                double cmin = obj.bbox[j].end[0];
                double cmax = obj.bbox[j].end[1];
                /* Determine allowed extent of object along axis {j}:: */
                double dmin = scene->dom[j].end[0];
                double dmax = scene->dom[j].end[1];
                if (j != 2)
                  { if (minSep >= 0)
                      { /* Should avoid {minSep} margins: */
                        dmin += minSep;
                        dmax -= minSep; 
                      }
                    else
                      { /* Object center only must be inside {dom[j]}: */
                        cmin = obj.ctr.c[j]; 
                        cmax = cmin;
                      }
                  }
                /* Check if coordinate {j} of object is contained in {scene->dom[j]}: */
                affirm(cmin > dmin, "object outside scene dom lo");
                affirm(cmax < dmax, "object outside scene dom hi");
              }
            NO++;
          }
      }
    assert(NO <= NO_max);
    if (NO < NO_max)
      { if (verbose) { fprintf(stderr, "generate only %d objects\n", NO); }
        /* Trim array: */
        objs = (multifok_scene_object_t*)notnull(realloc(objs, NO*sizeof(multifok_scene_object_t)), "no mem");
      }
    
    /* Store data in scene: */
    scene->NO = NO;
    scene->objs = objs;

    scene->tree = multifok_scene_build_tree(NO, objs);
    /* Check that the objects are still OK: */
    multifok_scene_check_objs(NO, objs, -1);
    
    /* Restore {random} state: */
    char *prestate = setstate(oldstate);
    assert(prestate == newstate);
  }

int32_t multifok_scene_choose_object_count(interval_t dom[], double rMin, double rMax, double minSep)
  {
    /* Compute the average area {aObj} of an object with radius in {[rMin _ rMax]}, accounting for min sep: */
    double rfat = (minSep >= 0 ? 0.5*minSep : 0);
    double rr0 = rMin + rfat; /* Min radius object occup. */
    double rr1 = rMax + rfat; /* Max radius object occup. */
    double aObj = M_PI*(rr1*rr1*rr1 - rr0*rr0*rr0)/(rr1 - rr0)/3; /* Average disk area. */
    double rObj = sqrt(aObj/M_PI); /* RMS average object radius. */
    
    /* Compute the {XY} area {aBox} available for those disks: */
    double wd[2];
    for (int32_t j = 0; j < 2; j++)
      { wd[j] = dom[j].end[1] - dom[j].end[0]; 
        if (minSep >= 0)
          { wd[j] -= 2*minSep; }
        else
          { wd[j] += 2*rObj; }
        demand (wd[j] >= 0, "dom too tight");
      }
    double aBox = wd[0]*wd[1]; /* Total useful area of dom. */

    /* Decide max number of objects. If non-overlapping, limited by area ratio, else more than that: */
    int32_t NO_max = (minSep >= 0 ? 1 : 2)*(int32_t)ceil(aBox/aObj); /* Max number of objects. */
    return NO_max;
  }

multifok_scene_object_t multifok_scene_object_throw(int32_t ID, interval_t dom[], double minSep, double rMin, double rMax)
  {
    demand(rMin <= rMax, "invalid radius interval");

    multifok_scene_object_t obj;
    obj.ID = ID;
    
    /* Choose disk or ball: */
    obj.flat = (drandom() < 0.5);
    
    /* Get the preliminary range of positions {[blo[0..2] _ bho[0..2]} for the object: */
    double blo[3], bhi[3];
    for (int32_t j = 0; j < 3; j++) 
      { blo[j] = dom[j].end[0]; 
        bhi[j] = dom[j].end[1];
        /* Shrink the box by a small amount: */
        double eps = 0.0001*(bhi[j] - blo[j]);
        blo[j] += eps;
        bhi[j] -= eps;
        /* If {minSep} is positive, shrink it by {minSep} in {X} and {Y}: */
        if ((minSep > 0) && ((j == 0) || (j == 1)))
          { blo[j] += minSep;
            bhi[j] -= minSep;
          }
      }

    /* Choose the object's radius: */
    obj.rad = rMin + drandom()*drandom()*(rMax - rMin);

    /* Reduce {rad} as needed so that the ball can be placed at all: */
    for (int32_t j = 0; j < 3; j++) 
      { double bw = bhi[j] - blo[j];
        if (((minSep >= 0) || (j == 2)) && (2*obj.rad > bw))
          { obj.rad = 0.4999*bw; }
      }

    /* Choose the center {ctr} and define the bounding box {bbox}: */
    for (int32_t j = 0; j < 3; j++)
      { double xlo = blo[j];
        double xhi = bhi[j];
        /* Decide the object's size along axis {j}: */
        double radj = ((j == 2) && obj.flat ? 0.0 : obj.rad);
        /* Decide if the object's {j}-projection must be fully inside {[blo[j] _ bhi[j]]}: */
        bool_t fully_inside = (j == 2) || (minSep >= 0);
        /* If it has to be fully inside, shrink the range for {ctr.c[j]}: */
        if (fully_inside) { xlo += radj; xhi -= radj; }
        demand (xlo < xhi, "dom is too tight");
        /* Choose the center's coordinate: */
        double cctr = xlo + drandom()*(xhi - xlo);
        obj.ctr.c[j] = cctr;
        /* Define the bounding box: */
        double cmin = cctr - radj;
        double cmax = cctr + radj;
        obj.bbox[j] = (interval_t) {{ cmin - FUDGE, cmax + FUDGE }};
      }
      
    /* Choose the background colors: */
    frgb_t tone_disk = (frgb_t){{ 1.0f, 0.3f, 0.0f }};
    frgb_t tone_ball = (frgb_t){{ 0.0f, 0.3f, 1.0f }};
    multifok_scene_throw_colors((obj.flat ? &tone_disk : &tone_ball), &(obj.bg), &(obj.fg));

    return obj;
  }
  
void multifok_scene_throw_colors(frgb_t *tone, frgb_t *bg, frgb_t *fg)
  {
    double wmax = +INF;
    for (int32_t ic = 0; ic < 3; ic++) 
      { double vi = tone->c[ic];
        bg->c[ic] = (float)(0.25*vi);
        wmax = fmin(wmax, 1.0 - bg->c[ic]);
      }
    for (int32_t ic = 0; ic < 3; ic++) 
      { fg->c[ic] = (float)(bg->c[ic] + wmax);
        assert(fg->c[ic] <= 1.00001);
        fg->c[ic] = (float)fmin(fg->c[ic], 1.0);
      }
  }

multifok_scene_tree_t *multifok_scene_build_tree(int32_t NO, multifok_scene_object_t objs[])
  {
    if (NO == 0)
      { return NULL; }
    else 
      { multifok_scene_tree_t *tr = (multifok_scene_tree_t *)notnull(malloc(sizeof(multifok_scene_tree_t)), "no mem");
        tr->child[0] = NULL;
        tr->child[1] = NULL;
        /* Choose {axis}, sort objects: */
        if (NO == 1)
          { tr->axis = -1; }
        else
          { /* We must choose {axis}, sort the objects, build {child[0],child[1]}. */
         
            /* Find the largest bounding dimension, {X} or {Y}, and split perp to that axis: */
            interval_t bboxXY[2]; /* Bounding box of all objects in {X} and {Y} only. */
            double wd[2]; /* Width of {bboxXY[j]}. */
            for (int32_t jxy = 0; jxy < 2; jxy++)
              { bboxXY[jxy] = (interval_t){{ +INF, -INF }}; /* Empty interval. */
                for (int32_t ko = 0; ko < NO; ko++)
                  { bboxXY[jxy] = interval_join(&(bboxXY[jxy]), &(objs[ko].bbox[jxy])); }
                wd[jxy] = interval_width(&(bboxXY[jxy]));
                assert(wd[jxy] > 0);
              }
            tr->axis = (wd[0] > wd[1] ? 0 : 1);
            
            /* Sort objects along that axis: */
            multifok_scene_sort_objects(NO, objs, tr->axis);
          }
        /* Pick root object {objs[mo]}: */
        int32_t mo = NO/2; /* Index of middle object */
        tr->obj = &(objs[mo]);
        
        /* Set {tr->bbox[0..2]} to the bounding box of the root object: */
        for (int32_t j = 0; j < 3; j++) { tr->bbox[j] = tr->obj->bbox[j]; }
            
        /* Recurse on the two halves: */
        for (int32_t ic = 0; ic < 2; ic++)
          { /* Define the objects {objs[ko_min..ko_max]} that go into {child[ic]}: */
            int32_t ko_min = (ic == 0 ? 0 : mo + 1);
            int32_t ko_max = (ic == 0 ? mo-1 : NO-1);
            int32_t NO_ch = ko_max - ko_min + 1;
            assert(NO_ch >= 0);
            if (NO_ch == 0)
              { tr->child[ic] = NULL; }
            else
              { /* Build {child[ic]}: */
                multifok_scene_tree_t *ch = multifok_scene_build_tree(NO_ch, &(objs[ko_min]));
                assert(ch != NULL);
                /* Merge {child[ic].bbox} into {bbox}: */
                for (int32_t j = 0; j < 3; j++) { tr->bbox[j] = interval_join(&(tr->bbox[j]), &(ch->bbox[j])); }
                tr->child[ic] = ch;
              }
          }
        /* Check that the objects are still OK: */
        /* multifok_scene_check_objs(NO, objs, -1); */
        return tr;
      }
  }
            
void multifok_scene_sort_objects(int32_t NO, multifok_scene_object_t objs[], int8_t axis)  
  { auto int32_t objcomp(const void*a, const void*b);
      /* Assumes {*a,*b} are {}. Compares their centers along {axis}, returns {-1,0,+1}. */
            
    /* multifok_scene_check_objs(NO, objs, -1); */
    if (NO > 1) { qsort((void*)objs, NO, sizeof(multifok_scene_object_t), objcomp); }
    /* multifok_scene_check_objs(NO, objs, axis); */

    return;
            
    int32_t objcomp(const void*a, const void*b)
      { multifok_scene_object_t *oba = (multifok_scene_object_t *)a;
        multifok_scene_object_t *obb = (multifok_scene_object_t *)b;
        double ca = oba->ctr.c[axis];
        double cb = obb->ctr.c[axis];
        if (ca < cb)
          { return -1; }
        else if (ca > cb)
          { return +1; }
        else
          { return 0; }
      }
  }

void multifok_scene_check_objs(int32_t NO, multifok_scene_object_t objs[], int8_t axis)
  { 
    /* Check id IDs are correct: */
    fprintf(stderr, "checking object IDs...\n");
    int32_t NO_max = 10000; /* Hope it is enough. */
    bool_t seen[NO_max];
    for (int32_t ko = 0; ko < NO_max; ko++) { seen[ko] = FALSE; }
    for (int32_t ko = 0; ko < NO; ko++) 
      { int32_t ID = objs[ko].ID; 
        assert((0 <= ID) && (ID < NO_max));
        assert(! seen[ID]);
        seen[ID] = TRUE;
      }
    if (axis >= 0)
      { /* Check if centers are sorted by increasing {axis} coordinate: */
        fprintf(stderr, "checking center ordering by axis %d...\n", axis);
        double cMax = -INF;
        for (int32_t ko = 0; ko < NO; ko++) 
          { double c = objs[ko].ctr.c[axis];
            assert(c >= cMax - FUDGE);
            cMax = c;
          }
      }
  }

void multifok_scene_ray_trace
  ( multifok_scene_t *scene,
    r3_t *p, 
    r3_t *d, 
    bool_t debug,
    multifok_scene_object_t **hob_P, 
    r3_t *hpt_P
  )
  {
    bool_t verbose = debug;
    multifok_scene_tree_t *tr = scene->tree;
    
    if (verbose) 
      { fprintf(stderr, "    " DASHES "\n");
        r3_gen_print(stderr, p, "%12.8f", "    tracing ray p = ( ", " ", " )");
        r3_gen_print(stderr, d, "%12.8f", "  d = ( ", " ", " )\n");
      }
      
    demand(d->c[2] == 1.0, "invalid direction vector");
    demand(p->c[2] >= 0, "invalid focus plane {Z}");
    
    /* Get highest ray-object hit: */
    multifok_scene_object_t *obj_hob; /* Highest object that was hit, or {NULL}. */
    r3_t obj_hpt; /* Hit point on {obj_hob}, or {(NAN,NAN,NAN)} */
    double zMin = 0.0 - FUDGE; /* Assume that the ray does not go below this {Z}. */
    int32_t level = 0;
    multifok_scene_ray_trace_tree(tr, p, d, zMin, debug, level, &obj_hob, &obj_hpt);
    if (verbose) 
      { if (obj_hob == NULL) 
          { fprintf(stderr, "    hit no object\n"); }
        else
          { fprintf(stderr, "    hit object %d", obj_hob->ID);
            r3_gen_print(stderr, p, "%12.8f", " at ( ", " ", " )\n");
          }
      }
    /* Compute ray hit with floor: */
    r3_t flo_hpt;
    multifok_scene_ray_trace_floor(scene->dom, scene->flatFloor, p, d, debug, &flo_hpt);
 
    /* Decide if ray hit floor or objects first. */
    multifok_scene_object_t *hob; /* Object that was hit, or {NULL}. */
    r3_t hpt;  /* Point where {hob} was hit. */
    if ((obj_hob == NULL) || (flo_hpt.c[2] > obj_hpt.c[2]))
      { /* Floor hit: */
        hob = NULL; hpt = flo_hpt;
      }
    else
      { /* Object hit: */
        hob = obj_hob;  hpt = obj_hpt;
      }

    (*hob_P) = hob;
    (*hpt_P) = hpt;
    if (verbose) { fprintf(stderr, "    " DASHES "\n"); }
  }

void multifok_scene_ray_trace_tree
  ( multifok_scene_tree_t *tr, 
    r3_t *p, 
    r3_t *d, 
    double zMin,
    bool_t debug, 
    int32_t level,
    multifok_scene_object_t **hob_P, 
    r3_t *hpt_P
  )
  { 
    bool_t verbose = debug;
    
    multifok_scene_object_t *hob = NULL;
    r3_t hpt = (r3_t){{ NAN, NAN, NAN }};
    
    if (tr != NULL) 
      { if (verbose) { fprintf(stderr, "    %*sray zMin = %12.8f tree zMax = %12.8f\n", 2*level, "", zMin, tr->bbox[2].end[1]); }
        if (zMin <= tr->bbox[2].end[1])
          { /* Get the ray's bounding box {ray_bbox[0..2]}: */
            interval_t ray_bbox[3];
            multifok_scene_get_ray_bbox(p, d, zMin, ray_bbox);

            /* Check if it it intersects the tree's bbox: */
            bool_t ok = TRUE; /* Set to false if ray's bbox is disjoint from tree's bbox. */
            for (int32_t j = 0; (j < 3) & ok; j++) 
              { if (ray_bbox[j].end[1] < tr->bbox[j].end[0]) { ok = FALSE; }
                if (ray_bbox[j].end[0] > tr->bbox[j].end[1]) { ok = FALSE; }
              }
            if (ok)
              { /* Bounding boxes intersect, we have a chance. */
                /* Ray-trace the root object: */
                if (verbose) { fprintf(stderr, "    %*strying root obj %d\n", 2*level, "", tr->obj->ID); }
                multifok_scene_object_t *obj = tr->obj;
                r3_t obj_hpt;
                bool_t obj_hit = multifok_scene_ray_trace_object(obj, p, d, zMin, debug, &obj_hpt);
                if (obj_hit)
                  { /* Root object hit replaces previous hit: */
                    hob = obj;
                    hpt = obj_hpt;
                    if (verbose || (hpt.c[2] < zMin-FUDGE))
                      { fprintf(stderr, "    %*shit root obj %d at Z = %12.8f zMin = %12.8f\n", 2*level, "", hob->ID, hpt.c[2], zMin); }
                    /* Adjust ray {zMin}: */
                    assert(hpt.c[2] >= zMin-FUDGE);
                    zMin = hpt.c[2];
                  }
                else
                  { if (verbose) { fprintf(stderr, "    %*smissed/rejected root obj %d\n", 2*level, "", obj->ID); } }

                /* We start with the child {ic} that seems closest to the ray in direction {axis}.
                  Then we ray-trace {child[ic]} and {child[1-ic]}, in that order, remembering the
                  highest hit as we go. */
                int8_t axis = tr->axis;
                int32_t ic = 0;
                if ((tr->child[0] != NULL) && (tr->child[1] != NULL))
                  { double obc = tr->obj->ctr.c[axis];
                    double chd[2]; /* Distance to coord {axis} of child box centers. */
                    for (int32_t jc = 0; jc < 2; jc++) 
                      { chd[jc] = fabs(obc - interval_mid(&(tr->child[jc]->bbox[axis]))); }
                    ic = (chd[0] < chd[1] ? 0 : 1);
                  }
                
                for (int32_t kc = 0; kc < 2; kc++) 
                  { /* Ray-trace {child[ic]} */
                    if (verbose) { fprintf(stderr, "    %*strying child %d\n", 2*level, "", ic); }
                    multifok_scene_object_t *ch_hob;
                    r3_t ch_hpt;
                    multifok_scene_ray_trace_tree(tr->child[ic], p, d, zMin, debug, level+1, &ch_hob, &ch_hpt);
                    if (ch_hob != NULL)
                      { /* We got a hit: */
                       if ((hob == NULL) || (ch_hpt.c[2] >= hpt.c[2]))
                          { hob = ch_hob;
                            hpt = ch_hpt;
                            if (verbose || (hpt.c[2] < zMin-FUDGE)) 
                              { fprintf(stderr, "    %*skept hit with child %d obj %d at Z = %12.8f zMin = %12.8f\n", 2*level, "", ic, hob->ID, hpt.c[2], zMin); }
                            /* Adjust ray {zMin}: */
                            assert(hpt.c[2] >= zMin-FUDGE);
                            zMin = hpt.c[2];
                          }
                        else
                          { if (verbose) { fprintf(stderr, "    %*srehected hit with child %d\n", 2*level, "", ic); } }
                      }
                    else
                      { if (verbose) { fprintf(stderr, "    %*smissed child %d\n", 2*level, "", ic); } }
                    /* Try the other child: */
                    ic = 1 - ic;
                  }
              }
          }
      }
    /* Return to caller: */
    if (verbose) 
      { if (hob != NULL) 
          { fprintf(stderr, "    %*sreturning hit with obj %d at Z = %12.8f\n", 2*level, "", hob->ID, hpt.c[2]); }
        else
          { fprintf(stderr, "    %*sreturning with no hit\n", 2*level, ""); }
      }
    (*hob_P) = hob;
    (*hpt_P) = hpt;
  }
        
void multifok_scene_get_ray_bbox(r3_t *p, r3_t *d, double zMin, interval_t bbox[])
  { 
    for (int32_t j = 0; (j < 3); j++) 
      { if (j == 2) 
          { bbox[j] = (interval_t){{ zMin, ZMAX }} ; }
        else
          { double cLo = p->c[j] + d->c[j]*(zMin - p->c[2]); /* Coordinate at bottom of ray. */
            double cHi = p->c[j] + d->c[j]*(ZMAX - p->c[2]); /* Coordinate at top of ray. */
            bbox[j] = (interval_t){{ fmin(cLo, cHi) - FUDGE, fmax(cLo, cHi) + FUDGE }};
          }
      }
  }

bool_t multifok_scene_ray_trace_object(multifok_scene_object_t *obj, r3_t *p, r3_t *d, double zMin, bool_t debug, r3_t *hpt_P)
  { 
    bool_t verbose = debug;
    
    /* Get the ray's bounding box */
    interval_t ray_bbox[3];
    multifok_scene_get_ray_bbox(p, d, zMin, ray_bbox);

    
    /* Check bbox intersection: */
    bool_t hit = TRUE; /* Set to false if bboxes are disjoint. */
    double eZ_hpt = NAN;  /* Value of {Z - p.c[2]} at hit point, if {hit} is true. */
    for (int32_t j = 0; j < 3; j++)
      { if (ray_bbox[j].end[1] < obj->bbox[j].end[0]) { hit = FALSE; }
        if (ray_bbox[j].end[0] > obj->bbox[j].end[1]) { hit = FALSE; }
      }
    if (hit)
      { /* Bounding boxes intersect, we have a chance. */

        /* Grab the relevant coordinates: */
        double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
        double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
        demand(dZ = 1.0, "invalid direction vector");

        r3_t *c = &(obj->ctr);
        double cX = c->c[0], cY = c->c[1], cZ = c->c[2];
        double r_obj = obj->rad, r2_obj = r_obj*r_obj;

        /* Let {ray(Z) = p + d*(Z - pZ)} be the point on the ray at height {Z}. */
        /* For generic {Z}, let {eZ} be {Z - pZ}, then {ray(Z) = p + d*eZ}. */
        /* For generic {Z}, let {s} be {ray(Z) - c}, that is {(p - c) + d*eZ}. */

        if (obj->flat)
          { /* Compute the point {h = ray(cZ)} where the ray hits the object's plane: */
            eZ_hpt = cZ - pZ; /* {Z} distance from object plane to {p}. */
            /* Compute displacement {sX,sY} from object's center to {h}: */ 
            double sX = pX + dX*eZ_hpt - cX; /* {X} position of ray hit rel to object ctr */
            double sY = pY + dY*eZ_hpt - cY; /* {Y} position of ray hit rel to object ctr */
            /* Check if ray hits object: */
            double r2_ray = sX*sX + sY*sY;
            hit = (r2_ray <= r2_obj);
          }
        else
          { /* Form the coefficients {A*eZ*eZ + B*eZ + C} of */
            /* {r2_obj = sX*sX + sY*sY + sZ*sZ - r2_obj} as function of {eZ}: */

            double LX = pX - cX;
            double LY = pY - cY;
            double LZ = pZ - cZ;

            /* sX = pX + dX*eZ - cX  = eZ*dX + LX; */
            /* sY = pY + dY*eZ - cY  = eZ*dY + LY; */
            /* sZ = Z - cZ           = eZ + LZ;    */

            double A = dX*dX + dY*dY + 1.0;
            double B = 2*(LX*dX + LY*dY + LZ);
            double C = LX*LX + LY*LY + LZ*LZ - r2_obj;

            double Delta = B*B - 4*A*C;
            hit = (Delta > 0);
            if (hit) { eZ_hpt = (sqrt(Delta) - B)/(2*A); }
          }
      }

    /* If we still have a hit, check its {Z} against {zMIn}: */
    if (hit) 
      { /* Ray hits object,  but maybe below {zMin}: */
        double hZ = eZ_hpt + p->c[2];
        char *what = (obj->flat ? "disk" : "ball");
        if (verbose) { fprintf(stderr, "    hit %s %d at Z = %+12.8f", what, obj->ID, hZ); } 
        if (hZ < zMin)
          { if (verbose) { fprintf(stderr, " (rejected, zMin = %12.8f)\n", zMin); }
            hit = FALSE;
          }
        else
          { if (verbose) { fprintf(stderr, " (accepted)\n"); }
            double hX = p->c[0] + d->c[0]*eZ_hpt;
            double hY = p->c[1] + d->c[1]*eZ_hpt;
            (*hpt_P) = (r3_t) {{ hX, hY, hZ }};
          }
       }
    return hit;
  }
  
void multifok_scene_ray_trace_floor(interval_t dom[], bool_t flatFloor, r3_t *p, r3_t *d, bool_t debug, r3_t *back_hpt_P)
  {     
    bool_t verbose = debug;
    
    /* Let {ray(Z) = p + d*(Z - pZ)} be the point on the ray at height {Z}. */

    /* Grab the relevant coordinates: */
    double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];

    /* Grab the relevant box coordinates: */
    double bXlo = dom[0].end[0];
    double bXhi = dom[0].end[1];
    double bZlo = dom[2].end[0];
    double bZhi = dom[2].end[1];
    
    double Z_hpt; /* {Z} of point where ray hits floor. */
    if (! flatFloor)
      { demand(dZ = 1.0, "invalid direction vector");
        /* The floor is a slanted plane that contains the */
        /* left bottom and top right edges of the {dom} box. */
        /* Determine the floor equation {(Z - bZ) = A*(X - bX)}: */
        /* But we must fudge it a bit because rays may be cast a {skosh} oustde the box. */
        double skosh = 1.0; /* Limit to box overflow due to pixel subsampling. */
        assert(pX >= bXlo - skosh);
        assert(pX <= bXhi + skosh);
        double bX = bXlo - skosh;
        double dom_dx = bXhi - bXlo + 2*skosh;
        double bZ = bZlo;
        double dom_dz = bZhi - bZlo;
        double A = dom_dz/dom_dx; /* Floor tilt, {dZ/dX}. */
        /* From the ray equation, {(X - pX) = dX*(Z - pZ)}. */
        /* Hence at the hit point, {Z - bZ + A*bX = A*X = A*(dX*Z - dX*pZ + pX)}. */
        /* That is, {Z*(1 - A*dX) = bZ + A*(pX - dX*pZ) - A*bX: */
        /* That is, {Z = (bZ + A*(pX - bX - dX*pZ))/(1 - A*dX): */
        double den = 1 - A*dX;
        demand(den > 1.0e-4, "floor too tilted for aperture");
        Z_hpt = (bZ + A*(pX - bX - dX*pZ))/den;
        /* We must clip the floor to the dom {Z} range anyway because of tilted rays. */
        if (Z_hpt < bZlo) { Z_hpt = bZlo; }
        if (Z_hpt > bZhi) { Z_hpt = bZhi; }
      }
    else
      { Z_hpt = bZlo; }
    
    if (verbose) { fprintf(stderr, "    hit floor Z = %+12.8f\n", Z_hpt); }

    /* For generic {Z}, let {eZ} be {Z - pZ}, then {ray(Z) = p + d*eZ}. */
    double eZ = Z_hpt - pZ; /* {Z} distance from hit point to {p}. */
    (*back_hpt_P) = (r3_t){{ pX + dX*eZ, pY + dY*eZ, Z_hpt }};
  }
   
frgb_t multifok_scene_compute_hit_color(multifok_scene_object_t *obj, r3_t *q, multifok_scene_pattern_t *pattern)
  { int32_t ID = (obj == NULL ? -1 : obj->ID);
    double r = pattern(q->c[0],q->c[1],q->c[1],ID);
    frgb_t clr;
    if (obj == NULL)
      { for (int32_t j = 0; j < 3; j++) { clr.c[j] = (float)r; } }
    else
      { clr = frgb_mix((1-r), &(obj->bg), r, &(obj->fg)); }
    return clr;
  }
 
#define multifok_scene_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

