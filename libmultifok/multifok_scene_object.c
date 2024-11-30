/* See {multifok_scene_object.h}. */
/* Last edited on 2024-10-29 13:25:27 by stolfi */

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
#include <interval_io.h>
#include <r3.h>
#include <frgb.h>
#include <r2.h>
#include <frgb_ops.h>
#include <jsrandom.h>

#include <multifok_scene.h>
#include <multifok_scene_object.h>
  
#define FUDGE (1.0e-6)
  /* Fudge amount to expand bounding boxes to account for roundoff. */

/* Shorter names for object types: */
#define ot_FLAT multifok_scene_object_type_FLAT
#define ot_RAMP multifok_scene_object_type_RAMP
#define ot_DISK multifok_scene_object_type_DISK
#define ot_BALL multifok_scene_object_type_BALL
#define ot_CONE multifok_scene_object_type_CONE

void multifok_scene_object_throw_colors(frgb_t *tone, frgb_t *bg, frgb_t *fg);
  /* Generates two random contrasting colors, with general hue depending on {warm}. */

/* RADIUS CHOSERS

  Each of the following procedures chooses randomly the half-size
  {rad[0..2]} of the bounding box of an object of the given type so
  that the box will fit in a box with radius {srad[0..2]}.
  
  If {strict}, the {XY} projection of the object will fit entirely in
  a box of size {srad[0] × srad[1]}, otherwise only the center will need to fit, and
  {rad[0]} and {rad[1]} may be as big as .
  
  The {XY} radii {rad[0]} and {rad[1]} will be in {[rMin _ rMax]}. The
  procedures will fail if {srad} is too small given {rMin},
  {strict}, and the object type. */
  
void multifok_scene_object_throw_radius_DISK
  ( double rMin,
    double rMax,
    double srad[],
    bool_t strict,
    double rad[]
  );

void multifok_scene_object_throw_radius_BALL
  ( double rMin,
    double rMax,
    double srad[],
    bool_t strict,
    double rad[]
  );

void multifok_scene_object_throw_radius_CONE
  ( double rMin,
    double rMax,
    double srad[],
    bool_t strict,
    double rad[]
  );

void   multifok_scene_object_throw_radius_round
  ( double rMin,
    double rMax,
    double hRel,
    double srad[],
    bool_t strict,
    double rad[]
  );
  /* Chooses the radius of an object with round {XY} projection
    ({DIS}, {BALL}, or {CONE}) whose height is {hRel}
    times the radius of that projection. */
  
/* IMPLEMENTATIONS */
             
void multifok_scene_object_print(FILE *wr, char *pref, multifok_scene_object_t *obj, char *suff)
  { if (pref != NULL) { fputs(pref, wr); }
    char *typeX = multifok_scene_object_type_to_string(obj->type);
    if (obj->ID >= 0) { fprintf(wr, "ID = %3d ", obj->ID); }
    fprintf(wr, "type = %s", typeX);
    r3_t ctr = multifok_scene_box_center(obj->bbox);
    r3_gen_print(wr, &(ctr), "%12.8f", "  ctr = ( ", " ", " )");
    r3_t rad = multifok_scene_box_radius(obj->bbox);
    r3_gen_print(wr, &(rad), "%12.8f", "  rad = ( ", " ", " )");
    frgb_print(wr, " bg = ( ", &(obj->bg), 3, "%5.3f", " )");
    frgb_print(wr, " fg = ( ", &(obj->fg), 3, "%5.3f", " )");
    if (suff != NULL) { fputs(suff, wr); }
  }

bool_t multifok_scene_object_XY_is_inside
  ( multifok_scene_object_t *obj,
    interval_t dom[],
    double margin
  )
  { for (uint32_t j = 0;  j < 2; j++)
      { /* Determine bounds {cmin,cmax} of object along axis {j}: */
        double omin, omax;
        /* Determine allowed extent of object along axis {j}:: */
        double dmin = dom[j].end[0];
        double dmax = dom[j].end[1];
        if (margin < 0)
          { /* Object center only must be inside {dom[j]}: */
            omin = interval_mid(&(obj->bbox[j])); 
            omax = omin;
          }
        else
          { /* Whole range {obj->bbox[j]} must be inside {dom[j]}: */
            omin = obj->bbox[j].end[0];
            omax = obj->bbox[j].end[1];
            if (margin > 0)
              { /* Should avoid {margin}: */
                dmin += margin;
                dmax -= margin; 
              }
          }
        if ((omin < dmin) || (omax > dmax)) { return FALSE; }
      }
    return TRUE;
  }
                   
bool_t multifok_scene_object_XY_overlap
  ( multifok_scene_object_t *obja,
    multifok_scene_object_t *objb,
    double minSep
  )
  { demand(minSep >=0, "invalid {minSep}");
  
    /* The {RAMP}  object is assumed to overlap with anything: */
    if ((obja->type == ot_RAMP) || (objb->type == ot_RAMP)) { return TRUE; }
    
    /* The {FLAT} object only overlaps if it cuts the other in {Z}: */
    if ((obja->type == ot_FLAT) || (objb->type == ot_FLAT))
      { if (obja->bbox[2].end[1] <= objb->bbox[2].end[0]) { return FALSE; }
        if (objb->bbox[2].end[0] >= obja->bbox[2].end[1]) { return FALSE; }
        return TRUE;
      }
    
    /* Neither {FLAT} nor {RAMP}. Check if boxes overlap in {X} and {Y}: */
    for (uint32_t j = 0;  j < 2; j++) 
      { /* Check if boxes have at least {minSep} separation along axis {j}: */
        if (obja->bbox[j].end[0] - objb->bbox[j].end[1] >= minSep) { return FALSE; } 
        if (objb->bbox[j].end[0] - obja->bbox[j].end[1] >= minSep) { return FALSE; } 
      }
      
    /* Boxes overlap, but maybe the objects don't: */
    bool_t around = (obja->type == ot_DISK) || (obja->type == ot_BALL) || (obja->type == ot_CONE);
    bool_t bround = (objb->type == ot_DISK) || (objb->type == ot_BALL) || (objb->type == ot_CONE);
    if (around && bround)
      { /* More detailed check using radii: */
        r3_t ctra = multifok_scene_box_center(obja->bbox);
        r3_t rada = multifok_scene_box_radius(obja->bbox);
        double rxya = (rada.c[0] + rada.c[1])/2;
        assert(fabs(rada.c[0] - rada.c[1]) < 1.0e-8*rxya); 
        
        r3_t ctrb = multifok_scene_box_center(objb->bbox);
        r3_t radb = multifok_scene_box_radius(objb->bbox);
        double rxyb = (radb.c[0] + radb.c[1])/2;
        assert(fabs(radb.c[0] - radb.c[1]) < 1.0e-8*rxyb); 
        
        double dxy = hypot(ctra.c[0] - ctrb.c[0], ctra.c[1] - ctrb.c[1]);
        return (dxy < rxya + minSep + rxyb);
      }
    else
      { /* That is it then: */
        return TRUE;
      }
  }

multifok_scene_object_t multifok_scene_object_background_make(interval_t dom[], bool_t flatFloor)
  {
    multifok_scene_object_t obj;
    obj.ID = -1;
    /* The bounding box is 3x the {dom} width in {X} and {Y}: */
    for (uint32_t j = 0;  j < 2; j++)
      { double rd = interval_rad(&(dom[j]));
        obj.bbox[j] = dom[j];
        interval_widen(&(obj.bbox[j]), rd);
      }
    if (flatFloor)
      { obj.type = ot_FLAT;
        /* The object's {Z} range is a singleton just above the scene's bottom: */
        double Zlo = dom[2].end[0] + 2*FUDGE; 
        obj.bbox[2] = (interval_t){{ Zlo, Zlo }};
      }
    else
      { obj.type = ot_RAMP;
        /* The object's {Z} range is just a tad less than the scene's {Z} range: */
        obj.bbox[2] = dom[2];
        interval_widen(&(obj.bbox[2]), -FUDGE);
      }
    obj.fg = (frgb_t){{ 0.000, 0.000, 0.000 }};
    obj.bg = (frgb_t){{ 1.000, 1.000, 1.000 }};
    return obj;
  }

multifok_scene_object_t multifok_scene_object_foreground_throw
  ( interval_t dom[],
    double margin,
    double rMin,
    double rMax,
    bool_t verbose
  )
  {
    if (verbose)
      { fprintf(stderr, "  throwing foreground object in");
        for (uint32_t j = 0;  j < 3; j++) 
          { if (j > 0) { fprintf(stderr, " ×"); }
            interval_gen_print(stderr, &(dom[j]), "%+8.3f", " [ ", " _ ", " ]");
          }
        fprintf(stderr, " margin = %+10.6f rMin = %12.5f rMax = %12.5f\n", margin, rMin, rMax);
      }
    demand(rMin <= rMax, "invalid radius interval");

    multifok_scene_object_t obj;
    obj.ID = -1;
   
    /* Choose disk or ball: */
    
    obj.type =  (drandom() < 0.33 ? ot_DISK : (drandom() < 0.50 ? ot_BALL : ot_CONE));
    bool_t strict = (margin >= 0); /* Whole {XY} projection must be inside {dom}. */
    
    /* Get the preliminary range of positions {sdom[0..2]} for the object: */
    interval_t sdom[3];
    for (uint32_t j = 0;  j < 3; j++) 
      { sdom[j] = dom[j]; 
        /* Shrink the range {sdom[j]} by a small amount: */
        double eps = 0.0001*interval_rad(&(sdom[j])) + FUDGE;
        sdom[j].end[0] += eps;
        sdom[j].end[1] -= eps;
        if ((j != 2) && strict)
          { /* Shrink range {sdom[j]} by {margin}: */
            sdom[j].end[0] += margin;
            sdom[j].end[1] -= margin;
          }
      }

    /* Define the object's size: */
    r3_t srad = multifok_scene_box_radius(sdom);
    double rad[3]; /* Half-size of box in each axis. */
    switch(obj.type)
      {
        case ot_DISK:
          multifok_scene_object_throw_radius_DISK(rMin, rMax, srad.c, strict, rad);
          break;

        case ot_BALL:
          multifok_scene_object_throw_radius_BALL(rMin, rMax, srad.c, strict, rad);
          break;

        case ot_CONE:
          multifok_scene_object_throw_radius_CONE(rMin, rMax, srad.c, strict, rad);
          break;
          
        case ot_RAMP:
        case ot_FLAT:
        default:
          affirm(FALSE, "invalid object type");
      }
      
    /* Define the object's center: */
    for (uint32_t j = 0;  j < 3; j++) 
      { if ((j == 2) || strict)
          { /* Reduce the the center range so that the object is all inside: */
            sdom[j].end[0] += (rad[j] + FUDGE);
            sdom[j].end[1] -= (rad[j] + FUDGE);
          }
        double rdomj = interval_rad(&(sdom[j]));
        if (verbose) 
          { fprintf(stderr, "    axis %d rad[j] = %10.6f available center range = ", j, rad[j]); 
            interval_gen_print(stderr, &(sdom[j]), "%+8.3f", " [ ", " _ ", " ]");
            fprintf(stderr, " radius = %10.6f\n", rdomj); 
          }
        demand(rdomj >= 0.0, "scene domain is too small"); 
        /* Pick a center {ctrj} in {sdom[j]}: */
        double ctrj = interval_mid(&(sdom[j])) + (2*drandom()-1)*rdomj;
        /* Define {obj.box[j]} from {ctrj} and {rad[j]}: */
        obj.bbox[j] = interval_from_mid_rad (ctrj, rad[j]);
      }
      
    /* Choose the object colors: */
    obj.fg = (frgb_t){{ (float)drandom(), 1.000, 0.600f }}; frgb_from_HTY(&(obj.fg));
    obj.bg = (frgb_t){{ (float)drandom(), 1.000, 0.200f }}; frgb_from_HTY(&(obj.bg));
     
    return obj;
  }
      
void multifok_scene_object_throw_radius_DISK
  ( double rMin,
    double rMax,
    double srad[],
    bool_t strict,
    double rad[]
  )
  { 
    multifok_scene_object_throw_radius_round(rMin, rMax, 0.0, srad, strict, rad);
  }
      
void multifok_scene_object_throw_radius_BALL
  ( double rMin,
    double rMax,
    double srad[],
    bool_t strict,
    double rad[]
  )
  { 
    multifok_scene_object_throw_radius_round(rMin, rMax, 2.0, srad, strict, rad);
  }
       
void multifok_scene_object_throw_radius_CONE
  ( double rMin,
    double rMax,
    double srad[],
    bool_t strict,
    double rad[]
  )
  {
    multifok_scene_object_throw_radius_round(rMin, rMax, 2.0, srad, strict, rad);
  }
 
void   multifok_scene_object_throw_radius_round
  ( double rMin,
    double rMax,
    double hRel,
    double srad[],
    bool_t strict,
    double rad[]
  )
  {
    /* Reduce {rMax} so that the height is at most {srad[2]}: */
    rMax = fmin(rMax, 2*srad[2]/hRel);
    
    /* Reduce {rMax} according to {srad[0..1],strict}: */
    for (uint32_t j = 0;  j < 2; j++) 
      { if (strict)
          { /* Radius of {XY} projection must be at most {srad[j]}: */
            rMax = fmin(rMax, srad[j]);
          }
        else
          { /* Radius of {XY} projection must be at most {2*srad[j]}: */
            rMax = fmin(rMax, 2*srad[j]);
          }
      }
    demand(rMax > 0.1*rMin, "domain too small for {rMin,strict}");
    double rObj;
    if (rMin >= rMax) 
      { rObj = rMax; }
    else
      { rObj = rMin + drandom()*(rMax - rMin); }
    rad[0] = rObj; rad[1] = rObj; rad[2] = hRel*rObj/2;
  }

char *multifok_scene_object_type_to_string(multifok_scene_object_type_t type)
  { switch (type)
      { case ot_FLAT: return "FLAT";
        case ot_RAMP: return "RAMP";
        case ot_DISK: return "DISK";
        case ot_BALL: return "BALL";
        case ot_CONE: return "CONE";
        default: demand(FALSE, "unrecognized object type");
      }
  }

#define multifok_scene_object_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

