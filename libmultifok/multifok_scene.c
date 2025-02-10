/* See {multifok_scene.h}. */
/* Last edited on 2025-02-10 03:36:19 by stolfi */

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

uint32_t multifok_scene_choose_throw_object_count(interval_t dom[], double wMinXY, double wMaxXY, double minSep);
  /* Chooses the ideal number of foreground objects that
    {multifok_scene_throw_foreground_objects} should try to put in a
    scene. The parameter {minSep} has the meaning described under that
    procedure. */

/* Shorter names for object types: */
#define ot_FLAT multifok_scene_object_type_FLAT
#define ot_RAMP multifok_scene_object_type_RAMP
#define ot_DISK multifok_scene_object_type_DISK
#define ot_BALL multifok_scene_object_type_BALL
#define ot_CONE multifok_scene_object_type_CONE
#define ot_PYRA multifok_scene_object_type_PYRA

#define DASHES "----------------------------------------------------------------------"
#define TILDES "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

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

void multifok_scene_add_floor
  ( multifok_scene_t *scene,
    multifok_scene_object_type_t type, 
    bool_t verbose
  )
  {
    demand((scene->NO == 0) && (scene->objs == NULL), "scene already has objects");
    
    if (verbose) 
      { char *typeX = multifok_scene_object_type_to_string(type);
        fprintf(stderr, "trying to add a floor object of type %s\n", typeX); 
      }
      
    interval_t *dom = scene->dom;
    multifok_scene_object_t obj = multifok_scene_object_background_make(type, dom, verbose);
    
    /* Add to scene: */
    uint32_t NO = scene->NO; /* Number of objects already generated. */
    scene->objs = retalloc(scene->objs, NO+1, multifok_scene_object_t);
    obj.ID = (multifok_scene_object_ID_t)NO; 
    scene->objs[NO] = obj; NO++;
    scene->NO = NO;
  }

void multifok_scene_add_foreground_object
  ( multifok_scene_t *scene,
    multifok_scene_object_type_t type,
    interval_t bbox[],
    frgb_t *fgGlo,
    frgb_t *bgGlo,
    frgb_t *fgLam,
    frgb_t *bgLam,
    bool_t verbose
  )
  {
    if (verbose)  { fprintf(stderr, "  entering {multifok_scene_add_foreground_object}\n"); }

    demand((scene->NO >= 1) && (scene->objs != NULL), "scene must have at least one object");
    demand((scene->objs[0].type == ot_FLAT) || (scene->objs[0].type == ot_RAMP), "the first object must be a floor");
    
    demand((type != ot_FLAT) && (type != ot_RAMP), "invalid object type");
    
    uint32_t NO = scene->NO; /* Number of objects already generated. */
    scene->objs = retalloc(scene->objs, NO+1, multifok_scene_object_t);
      
    multifok_scene_object_t obj = multifok_scene_object_foreground_make(type, bbox, fgGlo, bgGlo, fgLam, bgLam, verbose);
    
    /* Check object's containment in {scene.dom}: */
    interval_t *dom = scene->dom;
    for (uint32_t j = 0; j < 3; j++) 
      { assert(dom[j].end[0] <= dom[j].end[1]);
        assert(obj.bbox[j].end[0] <= obj.bbox[j].end[1]);
        if (j == 2)
          { bool_t inside = (LO(dom[j]) <= LO(obj.bbox[j])) && (HI(obj.bbox[j]) <= HI(dom[j]));
            demand(inside, "object's extends outside scene's {Z} range"); }
        else
          { double ctrj = interval_mid(&(obj.bbox[j]));
            bool_t inside = (LO(dom[j]) <= ctrj) && (ctrj <= HI(dom[j]));
            demand(inside, "object's center is outside scene's {XY} range"); 
          }
      }

    /* Add to scene: */
    obj.ID = (multifok_scene_object_ID_t)NO; 
    scene->objs[NO] = obj; NO++;
    scene->NO = NO;

    if (verbose)  { fprintf(stderr, "  exiting {multifok_scene_add_foreground_object}, ID = %d\n", obj.ID); }
  }

void multifok_scene_throw_foreground_objects
  ( multifok_scene_t *scene,
    double wMinXY, 
    double wMaxXY,
    frgb_t *fgGlo,
    double minSep, 
    bool_t verbose
  )
  {
    demand((scene->NO == 1) && (scene->objs != NULL), "scene is empty pr has 2+ objects");
    multifok_scene_object_t *ground = &(scene->objs[0]);
    if (verbose) 
      { fprintf(stderr, "floor object:\n");
        multifok_scene_object_print(stderr, 2, ground);
      }
    demand(ground->type == ot_FLAT, "invalid floor object type");
    
    /* Determine the max number of objects to generate: */
    uint32_t NO = scene->NO; /* Number of objects already generated. */
    uint32_t NO_max = NO + multifok_scene_choose_throw_object_count(scene->dom, wMinXY, wMaxXY, minSep);
    if (verbose) { fprintf(stderr, "trying to generate %d foreground objects\n", NO_max-NO); }

    scene->objs = retalloc(scene->objs, NO_max, multifok_scene_object_t);

    interval_t *dom = scene->dom;

    /* Generate the new objects {objs[1..]}. */
    /* If overlapping, every try is valid, otherwise we need many more tries: */
    uint32_t NT = (minSep >= 0 ? 50 : 1)*NO_max; /* Number of tries. */
    
    for (uint32_t kt = 0;  (kt < NT) && (NO < NO_max); kt++)
      { if (verbose) { fprintf(stderr, "" TILDES "\n"); }
        /* Generate a random foreground object {obj}: */
        multifok_scene_object_t obj = multifok_scene_object_foreground_throw(dom, minSep, wMinXY, wMaxXY, fgGlo, verbose);
        if (verbose) { fprintf(stderr, "trying to add object ...\n"); }

        demand((obj.type != ot_FLAT) && (obj.type != ot_RAMP), "invalid object type");
        
        /* Check containment in {Z}: */
        assert(dom[2].end[0] <= obj.bbox[2].end[0]);
        assert(obj.bbox[2].end[0] <= obj.bbox[2].end[1]);
        assert(obj.bbox[2].end[1] <= dom[2].end[1]);
        
        int32_t ko_overlap = -1; /* Index of object overlapped by {obj}, or {-1}. */
        if (minSep >= 0)
          { /* Reject {obj} if it overlaps previous foreground objects in {X} and {Y}: */
            /* Check for {XY} overlaps, skipping the floor: */
            for (int32_t ko = 0;  (ko < NO) && (ko_overlap < 0); ko++)
              { multifok_scene_object_t *objk = &(scene->objs[ko]);
                if ((objk->type != ot_FLAT) && (objk->type != ot_RAMP))
                  { bool_t overlap = multifok_scene_object_XY_overlap(&obj, objk, minSep);
                    if (overlap) { ko_overlap = ko; }
                  }
              }
            assert(multifok_scene_object_XY_is_inside(&obj, scene->dom, minSep));
          }
        if (ko_overlap < 0)
          { if (verbose) { fprintf(stderr, "accepted, ID = %d\n", NO); }
            assert(NO < NO_max);
            obj.ID = (multifok_scene_object_ID_t)NO;
            scene->objs[NO] = obj;  NO++;
          }
        else
          { if (verbose) { fprintf(stderr, "overlaps %d, rejected\n", ko_overlap); } }
      }
    if (verbose) { fprintf(stderr, "" TILDES "\n"); }
        
    assert(NO <= NO_max);
    if (NO < NO_max)
      { if (verbose) { fprintf(stderr, "generated only %d objects\n", NO); }
        /* Trim array: */
        scene->objs = retalloc(scene->objs, NO, multifok_scene_object_t);
      }
    scene->NO = NO;
  }

uint32_t multifok_scene_choose_throw_object_count
  ( interval_t dom[],
    double wMinXY,
    double wMaxXY,
    double minSep
  )
  {
    /* Compute the average area {aObj} of an object with width in {[wMinXY _ wMaxXY]}, accounting for min sep: */
    double rfat = (minSep >= 0 ? 0.5*minSep : 0);
    double rr0 = 0.5*wMinXY + rfat; /* Min radius object occup. */
    double rr1 = 0.5*wMaxXY + rfat; /* Max radius object occup. */
    double avgPI = (3*M_PI + 4)/4; /* Average value of {PI} assuming equal chances. */
    double areaObj = avgPI*(rr1*rr1*rr1 - rr0*rr0*rr0)/(rr1 - rr0)/3; /* Average {XY} proj area. */
    double rObj = sqrt(areaObj/avgPI); /* RMS average object radius. */
    
    /* Compute the {XY} area {areaBox} available for those objects: */
    double wd[2];
    for (uint32_t j = 0;  j < 2; j++)  /* Note 0 and 1 only. */
      { wd[j] = dom[j].end[1] - dom[j].end[0]; 
        if (minSep >= 0)
          { /* Objects wholly inside: */ wd[j] -= 2*minSep; }
        else
          { /* Objects may go outside: */ wd[j] += 2*rObj; }
        demand (wd[j] >= 0, "dom too tight");
      }
    double areaBox = wd[0]*wd[1]; /* Total useful area of dom. */

    /* Decide max number of foreground objects {NO_gen}. */
    /* If non-overlapping, limited by area ratio, else more than that: */
    uint32_t NO_gen = (minSep >= 0 ? 1 : 2)*(uint32_t)ceil(areaBox/areaObj);
    /* Ensure at least one foreground object: */
    if (NO_gen <= 0) { NO_gen = 1; }
    
    return NO_gen;
  }
  
void multifok_scene_print(FILE *wr, multifok_scene_t *scene)
  { multifok_scene_print_box(wr, "dom = ", scene->dom, "\n");
    fprintf(stderr, "scene has %d objects:\n", scene->NO);
    fputs("  " DASHES "\n", wr);
    for (int32_t i = 0; i < scene->NO; i++)
      { multifok_scene_object_print(wr, 2, &(scene->objs[i]));
        fputs("  " DASHES "\n", wr);
      }
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

