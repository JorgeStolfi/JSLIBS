/* See {multifok_scene_object.h}. */
/* Last edited on 2025-02-09 00:06:32 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <box.h>
#include <interval_io.h>
#include <r3.h>
#include <rn.h>
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
#define ot_PYRA multifok_scene_object_type_PYRA

void multifok_scene_object_throw_colors(frgb_t *tone, frgb_t *bg, frgb_t *fg);
  /* Generates two random contrasting colors, with general hue depending on {warm}. */
   
void multifok_scene_object_foreground_choose_bbox
  ( multifok_scene_object_type_t type, 
    interval_t dom[],
    double wMinXY,
    double wMaxXY,
    double margin,
    interval_t bbox[]
  );
  /* Chooses a bounding box {bbox[0..2]} for an object of the given {type}
    and scene domain {dom[0..2]}.  The box {X} and {Y] width will be 
    randomly chosen in {[wMinXY _ wMaxXY]}.   
    
    If {margin} is positive or zero, the chosen box will be entirely
    contained in {dom}, and its {XY} projection will be at least
    {margin} away from the sides of {dom}. If {margin} is negative, the
    {XY} projection of the box may extend up to {0.5*wMaxXY} outside {dom},
    but its center will lie inside {dom}.
    
    The shape of the box will
    depend on {type}. See {multifok_scene_object_adjust_bbox_size_and_shape}
    below.
    
    The procedure will fail if {sdom} is too small given {wMinXY} and
    {margin}. */

/* OBJECT SIZE ADJUSTERS

  Each of the following procedures will reduce {bbox[0..2]} as little as
  needed to make its aspect ratio conform to the specific object type. They
  fail if the reduced {X} or {Y} width of the box is less than {wdMIn}. */
 
void multifok_scene_object_adjust_bbox_size_and_shape(multifok_scene_object_type_t type, double wXY, interval_t bbox[]);
  /* Shrinks {bbox[0..2]}if needed so that its {X} and {Y} sizes are at
    most {wXY}, and its shape (the ratios between then three
    dimensions) is appropriate for an object of the given {type}. For
    instance, of {type} is {ot_BALL} the box will be a cube; if {type}
    is {ot_DISK}, the box will be a square in {XY} with zero height.
    
    If {bbox} needs to be shrunk along any axis {j}, its new range {bbox[j]}
    will be randomly placed inside the original {bbox[j]} */

void multifok_scene_object_adjust_bbox_size_and_shape_round_or_square(double width[], double wXY, double hRel);
  /* Reduces {width[0..2]} if needed so that {width[0]} and {width[1]} are equal and at most {wXY},
    and {width[2]} is {hRel} times that value. */

#define TILDES "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

/* IMPLEMENTATIONS */
             
void multifok_scene_object_print(FILE *wr, uint32_t indent, multifok_scene_object_t *obj)
  { char *typeX = multifok_scene_object_type_to_string(obj->type);
    if (obj->ID >= 0) { fprintf(wr, "%*sID = %3d ", indent, "", obj->ID); }
    fprintf(wr, "type = %s\n", typeX);
    fprintf(wr, "%*s", indent, "");
    box_gen_print(wr, 3, obj->bbox, "%12.8f", "bbox = ", " × ", "\n");
    fprintf(wr, "%*s", indent, "");
    r3_t ctr, rad; box_center_and_radii(3, obj->bbox, ctr.c, rad.c);
    r3_gen_print(wr, &ctr, "%12.8f", "ctr = ( ", " ", " )");
    r3_gen_print(wr, &rad, "%12.8f", "  rad = ( ", " ", " )\n");
    fprintf(wr, "%*s", indent, "");
    frgb_print(wr, "gloss =   ( ", &(obj->bgGlo), 3, "%5.3f", " )");
    frgb_print(wr, " -- ( ", &(obj->fgGlo), 3, "%5.3f", " )\n");
    fprintf(wr, "%*s", indent, "");
    frgb_print(wr, "lambedo = ( ", &(obj->bgLam), 3, "%5.3f", " )");
    frgb_print(wr, " -- ( ", &(obj->fgLam), 3, "%5.3f", " )\n");
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
  
    /* The {ot_RAMP}  object is assumed to overlap with anything: */
    if ((obja->type == ot_RAMP) || (objb->type == ot_RAMP)) { return TRUE; }
    
    /* The {ot_FLAT} object only overlaps if it cuts the other in {Z}: */
    if ((obja->type == ot_FLAT) || (objb->type == ot_FLAT))
      { if (obja->bbox[2].end[1] <= objb->bbox[2].end[0]) { return FALSE; }
        if (objb->bbox[2].end[0] >= obja->bbox[2].end[1]) { return FALSE; }
        return TRUE;
      }
    
    /* Neither {ot_FLAT} nor {ot_RAMP}. Check if boxes overlap in {X} and {Y}: */
    for (uint32_t j = 0;  j < 2; j++) 
      { /* Check if boxes have at least {minSep} separation along axis {j}: */
        if (obja->bbox[j].end[0] - objb->bbox[j].end[1] >= minSep) { return FALSE; } 
        if (objb->bbox[j].end[0] - obja->bbox[j].end[1] >= minSep) { return FALSE; } 
      }
      
    /* Boxes overlap in {XY}, but maybe the objects don't: */
    bool_t around = (obja->type == ot_DISK) || (obja->type == ot_BALL) || (obja->type == ot_CONE);
    bool_t bround = (objb->type == ot_DISK) || (objb->type == ot_BALL) || (objb->type == ot_CONE);
    if (around && bround)
      { /* More detailed check using radii: */
        /* multifok_scene_print_box(stderr, "    obja.bbox = ", obja->bbox, "\n"); */
        /* multifok_scene_print_box(stderr, "    objb.bbox = ", objb->bbox, "\n"); */
        
        r3_t ctra, rada; box_center_and_radii(3, obja->bbox, ctra.c, rada.c);
        double rxya = (rada.c[0] + rada.c[1])/2;
        assert(isfinite(rxya));
        assert(fabs(rada.c[0] - rada.c[1]) < 1.0e-8*rxya); 
        
        r3_t ctrb, radb; box_center_and_radii(3, objb->bbox, ctrb.c, radb.c);
        double rxyb = (radb.c[0] + radb.c[1])/2;
        assert(isfinite(rxyb));
        assert(fabs(radb.c[0] - radb.c[1]) < 1.0e-8*rxyb); 
        
        double dxy = hypot(ctra.c[0] - ctrb.c[0], ctra.c[1] - ctrb.c[1]);
        return (dxy < rxya + minSep + rxyb);
      }
    else
      { /* That is it then: */
        return TRUE;
      }
  }

multifok_scene_object_t multifok_scene_object_background_make
  ( multifok_scene_object_type_t type,
    interval_t dom[],
    bool_t verbose
  )
  {
    multifok_scene_object_t obj;
    obj.ID = multifok_scene_object_ID_NONE;
    /* The bounding box is 3x the {dom} width in {X} and {Y}: */
    for (uint32_t j = 0;  j < 2; j++)
      { double rd = interval_rad(&(dom[j]));
        obj.bbox[j] = dom[j];
        interval_widen(&(obj.bbox[j]), rd);
      }
    if (type == ot_FLAT)
      { /* The object's {Z} range is a singleton just above the scene's bottom: */
        double Zlo = dom[2].end[0] + 2*FUDGE; 
        obj.bbox[2] = (interval_t){{ Zlo, Zlo }};
      }
    else if (type == ot_RAMP)
      { /* The object's {Z} range is just a tad less than the scene's {Z} range: */
        obj.bbox[2] = dom[2];
        interval_widen(&(obj.bbox[2]), -FUDGE);
      }
    else
      { demand(FALSE, "invalid background object type"); }
      
    obj.type = type;
    /* No gloss: */
    obj.fgGlo = (frgb_t){{ 0.000, 0.000, 0.000 }};
    obj.bgGlo = (frgb_t){{ 0.000, 0.000, 0.000 }};
    /* Black-to-white lambedo: */
    obj.fgLam = (frgb_t){{ 0.000, 0.000, 0.000 }};
    obj.bgLam = (frgb_t){{ 1.000, 1.000, 1.000 }};
    if (verbose) { multifok_scene_object_print(stderr, 2, &obj); }
    return obj;
  }
      
multifok_scene_object_t multifok_scene_object_foreground_make
  ( multifok_scene_object_type_t type,
    interval_t bbox[],
    frgb_t *fgGlo,
    frgb_t *bgGlo,
    frgb_t *fgLam,
    frgb_t *bgLam,
    bool_t verbose
  )
  { 
    multifok_scene_object_t obj;
    obj.ID = multifok_scene_object_ID_NONE;
    obj.type = type;
    box_copy(3, bbox, obj.bbox);
    multifok_scene_object_adjust_bbox_size_and_shape(type, +INF, obj.bbox);
    obj.fgGlo = (fgGlo == NULL ? (frgb_t){{ 0,0,0 }} : (*fgGlo));
    obj.bgGlo = (bgGlo == NULL ? (frgb_t){{ 0,0,0 }} : (*bgGlo));
    obj.fgLam = (fgLam == NULL ? (frgb_t){{ 0,0,0 }} : (*fgLam));
    obj.bgLam = (bgLam == NULL ? (frgb_t){{ 0,0,0 }} : (*bgLam));
    if (verbose) { multifok_scene_object_print(stderr, 2, &obj); }
    return obj;
  }

multifok_scene_object_t multifok_scene_object_foreground_throw
  ( interval_t dom[],
    double margin,
    double wMinXY,
    double wMaxXY,
    frgb_t *fgGlo,
    bool_t verbose
  )
  {
    if (verbose)
      { fprintf(stderr, "throwing foreground object in");
        multifok_scene_print_box(stderr, "", dom, "\n");
        fprintf(stderr, "margin = %+10.6f width range = [ %12.6f _ %12.6f ]\n", margin, wMinXY, wMaxXY);
      }
   
    /* Choose object type: */
    
    multifok_scene_object_type_t type;
    switch (int32_abrandom(0,3))
      { case 0: type = ot_DISK; break;
        case 1: type = ot_BALL; break;
        case 2: type = ot_CONE; break;
        case 3: type = ot_PYRA; break;
        default: assert(FALSE);
      }
     
    /* Choose the object position and size: */
    interval_t bbox[3];
    multifok_scene_object_foreground_choose_bbox(type, dom, wMinXY, wMaxXY, margin, bbox);

    /* Choose the object colors: */
    frgb_t bgGlo = (frgb_t){{ 0.000f, 0.000f, 0.000f }};
    frgb_t fgLam = (frgb_t){{ (float)drandom(), 1.000, 0.600f }}; frgb_from_HTY(&(fgLam));
    frgb_t bgLam = (frgb_t){{ (float)drandom(), 1.000, 0.200f }}; frgb_from_HTY(&(bgLam));

    multifok_scene_object_t obj = multifok_scene_object_foreground_make(type, bbox, fgGlo, &bgGlo, &fgLam, &bgLam, verbose);
    return obj;
  }

void multifok_scene_object_foreground_choose_bbox
  ( multifok_scene_object_type_t type, 
    interval_t dom[],
    double wMinXY,
    double wMaxXY,
    double margin,
    interval_t bbox[]
  )
  {
    bool_t debug = FALSE;
    demand(wMinXY <= wMaxXY, "invalid width interval");
    
    box_copy(3, dom, bbox);
    /* Adjust {X} and {Y} ranges according to {margin}: */
    for (uint32_t j = 0; j < 2; j++)  /* Note, only 0 and 1. */
      { interval_t *bboxj = &(bbox[j]);
        if (margin >= 0)
          { /* Whole object must be inside {dom} with margin: */ 
            interval_widen(bboxj, -margin);
          }
        else
          { /* Object may extend outside {dom}: */
            interval_widen(bboxj, 0.5*wMaxXY);
          }
        demand(! interval_is_empty(bboxj), "failed to allocate the object in domain");
      }
    
    /* Shrink the box by a small amount, just in case: */
    for (uint32_t j = 0; j < 3; j++) 
      { double eps = 0.0001*interval_rad(&(dom[j])) + FUDGE;
        interval_widen(&(bbox[j]), -eps);
      }
      
    /* Adjust the {bbox} size and aspect ratios: */
    if (debug) { box_gen_print(stderr, 3, bbox, "%12.8f", "bbox befor adjstment = ", " × ", "\n"); }
    double wXY = dabrandom(wMinXY, wMaxXY);
    multifok_scene_object_adjust_bbox_size_and_shape(type, wXY, bbox);
    if (debug) { box_gen_print(stderr, 3, bbox, "%12.8f", "bbox after adjstment = ", " × ", "\n"); }
    
    /* Ensure that the center is inside {dom}, in case {margin} was negative: */
    for (uint32_t j = 0; j < 3; j++) /* Note, only 0 and 1. */
      { interval_t *bboxj = &(bbox[j]);
        double ctrj, radj;
        interval_mid_rad(bboxj, &ctrj, &radj);
        if (j != 2) { demand(radj > 0.5*wMinXY, "domain too small to place the object"); }
        ctrj = interval_project(&(dom[j]), ctrj);
        (*bboxj) = interval_from_mid_rad(ctrj, radj);
      }
    if (debug) { multifok_scene_print_box(stderr, "    chosen box = ", bbox, "\n");  }
  }
    
void multifok_scene_object_adjust_bbox_size_and_shape(multifok_scene_object_type_t type, double wXY, interval_t bbox[])
  {
    double width[3];
    box_widths(3, bbox, width); 
    switch(type)
      {
        case ot_DISK:
          multifok_scene_object_adjust_bbox_size_and_shape_round_or_square(width, wXY, 0.0);
          break;

        case ot_BALL:
          multifok_scene_object_adjust_bbox_size_and_shape_round_or_square(width, wXY, 1.0);
          break;

        case ot_CONE:
          multifok_scene_object_adjust_bbox_size_and_shape_round_or_square(width, wXY, 1.0);
          break;
          
        case ot_PYRA:
          multifok_scene_object_adjust_bbox_size_and_shape_round_or_square(width, wXY, 0.5);
          break;
          
        case ot_RAMP:
        case ot_FLAT:
        default:
          affirm(FALSE, "invalid object type");
      }

    /* Reposition shrunk box inside original box: */
    for (int32_t j = 0; j < 3; j++)
      { interval_t bbj = bbox[j];
        double slack = interval_width(&bbj) - width[j];
        assert(slack >= -1.0e-10);
        if (slack > 0)
          { /* Pick new range: */
            bbox[j].end[0] = LO(bbj) + slack*drandom();
            bbox[j].end[1] = fmin(HI(bbj), bbox[j].end[0] + width[j]);
          }
      }
  }
 
void multifok_scene_object_adjust_bbox_size_and_shape_round_or_square(double width[], double wXY, double hRel)
  {
    wXY = fmin(wXY, fmin(width[0], width[1]));
    width[0] = width[1] = wXY;
    if (hRel*wXY < width[2]) 
      { width[2] = hRel*wXY; } 
    else if (hRel*wXY > width[2])
      { wXY = width[2]/hRel; }
  }

char *multifok_scene_object_type_to_string(multifok_scene_object_type_t type)
  { switch (type)
      { case ot_FLAT: return "FLAT";
        case ot_RAMP: return "RAMP";
        case ot_DISK: return "DISK";
        case ot_BALL: return "BALL";
        case ot_CONE: return "CONE";
        case ot_PYRA: return "PYRA";
        default: demand(FALSE, "unrecognized object type");
      }
  }

multifok_scene_object_type_t multifok_scene_object_type_from_string(char *str)
  { 
    multifok_scene_object_type_t type;
    if      (strcmp(str, "FLAT") == 0) { type = ot_FLAT; }
    else if (strcmp(str, "RAMP") == 0) { type = ot_RAMP; }
    else if (strcmp(str, "DISK") == 0) { type = ot_DISK; }
    else if (strcmp(str, "BALL") == 0) { type = ot_BALL; }
    else if (strcmp(str, "CONE") == 0) { type = ot_CONE; }
    else if (strcmp(str, "PYRA") == 0) { type = ot_PYRA; }
    else { demand(FALSE, "no such object type"); }
    return type;
  }

void multifok_scene_object_finish
  ( multifok_scene_object_t *obj,
    r3_t *q,
    multifok_pattern_double_proc_t *patternGlo,
    frgb_t *sGlo_P,
    multifok_pattern_double_proc_t *patternLam,
    frgb_t *sLam_P
  )
  { 
    demand(sGlo_P != NULL, "null {sGlo_P}");
    demand(sLam_P != NULL, "null {sLam_P}");
    /* bool_t debug = FALSE; */
    frgb_t sGlo, sLam;
    if (obj == NULL)
      { sGlo = (frgb_t){{ 0,0,0 }};
        sLam = (frgb_t){{ 0.500, 0.500, 0.500 }};
      }
    else
      { interval_t *bbox = obj->bbox; 
        /* Express {q} relative to object: */
        r3_t u; 
        for (uint32_t j = 0;  j < 3; j++) 
          { u.c[j] = q->c[j] - interval_mid(&(bbox[j])); }
        int32_t ID = obj->ID;
        assert(ID != multifok_scene_object_ID_NONE);
        
        /* Rotate {u} about {Z} by an angle that is a function of {ID}: */
        if (patternGlo == NULL)
          { sGlo = obj->fgGlo; }
        else
          { r3_t uGlo; r3_rot_axis(&u, 0, 1, 3.0*ID + 1.0, &uGlo);
            r3_scale(0.75, &uGlo, &uGlo);
            double rGlo = patternGlo(&uGlo);
            sGlo = frgb_mix((1-rGlo), &(obj->bgGlo), rGlo, &(obj->fgGlo));
          }
          
        if (patternLam == NULL)
          { sLam = obj->fgLam; }
        else
          { r3_t uLam; r3_rot_axis(&u, 0, 1, 5.5*ID + 1.0, &uLam);
            r3_scale(1.5, &uLam, &uLam);
            double rLam = patternLam(&uLam);
            sLam = frgb_mix((1-rLam), &(obj->bgLam), rLam, &(obj->fgLam));
          }
      }
    (*sGlo_P) = sGlo;
    (*sLam_P) = sLam;
  }

#define multifok_scene_object_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

