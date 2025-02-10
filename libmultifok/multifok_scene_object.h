/* Objects for a {multifok_scene_t}. */
/* Last edited on 2025-02-08 17:26:26 by stolfi */

#ifndef multifok_scene_object_H
#define multifok_scene_object_H

#include <stdint.h>

#include <r3.h>
#include <bool.h>
#include <interval.h>
#include <frgb.h>

#include <multifok_pattern.h>

typedef enum 
  { multifok_scene_object_type_FLAT,   /* Flat floor near bottom {Z}of domain. */
    multifok_scene_object_type_RAMP,   /* Ramp where {Z} increases with {X}. */
    multifok_scene_object_type_DISK,   /* A flat horizontal disk. */
    multifok_scene_object_type_BALL,   /* A ball. */
    multifok_scene_object_type_CONE,   /* A cone. */
    multifok_scene_object_type_PYRA,   /* A square pyramid. */
  } multifok_scene_object_type_t;
  /* Type of scene object.  In the comments below, {ot_{XXXX}} means
    {multifok_scene_object_type_{XXXX}}. */
    
#define multifok_scene_object_type_FRST multifok_scene_object_type_FLAT
#define multifok_scene_object_type_LAST multifok_scene_object_type_PYRA
      
#define WD_RAMP 0.70
  /* With of {ot_RAMP} inclined part relative to scene width. */ 

typedef int32_t multifok_scene_object_ID_t;
  /* Type of an object ID. */
  
#define multifok_scene_object_ID_NONE (-1)
  /* An invalid or undefined ID value. */

typedef struct multifok_scene_object_t
  { multifok_scene_object_ID_t ID;     /* Object ID, used e.g. for patterning. */
    multifok_scene_object_type_t type; /* Type of object. */
    interval_t bbox[3];                /* Bounding box for the object. */
    frgb_t bgGlo;                      /* "Background" value of object's gloss. */
    frgb_t fgGlo;                      /* "Foreground" color of object's gloss. */
    frgb_t bgLam;                      /* "Background" color of object's lambedo. */
    frgb_t fgLam;                      /* "Foreground" color of object's lambedo. */
  } multifok_scene_object_t;
  /* A record {obj} that specifies a an object in the scene.
    The {bbox[0..2]} is a tight 3D bounding box for the object.
    
    If {type} is {ot_BALL} the object is a ball inscribed in the specified
    {bbox}, which should be a cube.
    
    If {type} is {ot_DISK} the object is a horizontal disk. Its {XY}
    projection is inscribed in the rectangle {bbox[0]×bbox[1]}, which
    should be a square, and its {Z} coordinate is the midpoint of
    {bbox[2]}, which should have very small height.
    
    If {type} is {ot_CONE} the object is a cone. Its {XY}
    projection is inscribed in the rectangle {bbox[0]×bbox[1]}, which
    should be a square, and its {Z} coordinates spans {bbox[2]}, with the 
    base at bottom.
    
    If {type} is {ot_PYRA} the object is a square pyramid. Its {XY}
    projection is the rectangle {bbox[0]×bbox[1]}, which
    should be a square, and its {Z} coordinates spans {bbox[2]}, with the 
    base at bottom.
    
    If {type} is {ot_FLAT}, the object is a horizontal rectangle whose {XY}
    projection is {bbox[0]×bbox[1]}, and whose {Z} coordinate is the
    midpoint of {bbox[2]}, which should have very small height.
        
    If {type} is {ot_RAMP}, the object consists of three rectangles that
    span the whole range {bbox[1]} in the {Y} direction and one third of
    the range {bbox[0]} in the {X} direction. The middle rectangle is
    parallel to the {Y} axis and rises from {Zlo} when {X} is {Xlo} to
    {Zhi} when {X} is {Xhi}, where {Xhi-Xlo} is {WD_RAMP} times the
    width of {bbox[0]}, and {Zlo,Zhi} are the limits of {bbox[2]}. The
    first rectangle is horizontal at height {Zlo}, and the third one is
    horizontal at height {Zhi}.
    
    The Lambertian self-color (albedo, lambedo) of the object at each point {q}
    will be a varying combination of colors {obj.bgLam} and {obj.fgLam}.
    The glossy component of its finish will be a varying combination of
    {obj.bgGlo} and {obj.fgGlo}. See {multifok_scene_object_finish} for
    more info. */ 

void multifok_scene_object_print(FILE *wr, uint32_t indent, multifok_scene_object_t *obj);
  /* Writes a readable description of object {obj} to {wr}, in one line,
    preceded by {pref} and folloed by {suff}. */

void multifok_scene_object_finish
  ( multifok_scene_object_t *obj,
    r3_t *q,
    multifok_pattern_double_proc_t *patternGlo,
    frgb_t *sGlo_P,
    multifok_pattern_double_proc_t *patternLam,
    frgb_t *sLam_P
  );
  /* Computes the glossy coefficient {sGlo} and the lambedo 
    (Lambertian albedo, intrinsic color) {sLam} of the surface of object {obj} at the point {q}. 
    returns them in {*sGlo_P} and {*sLam_P}. 
    
    If {obj} is not {NULL}, the procedure evaluates the glossy coefficient triple
    {sGlo} at the point {q} by evaluating {r=patternGlo(dq)}, where
    {dq} is the vector {q - obj.ctr} rotated and/or translated by
    an angle that depends on the object's ID.  Then {sGlo} is
    obtained by interpolating between {obj.bgGlo} and {obj.fgGlo}
    with ratio {rGlo}. 
    
    The lambedo coeff triple {sLam} is similary obtained by
    interpolating between {obj.bgLam} and {obj.fgLam} with the ratio {rLam}
    obtained from {patternLam}.
    
    If {obj} is {NULL}, returns {(0,0,0)} as the normal, {(0,0,0)}
    as the glossy component {sGlo}, and uniform {50%} gray as the 
    lambedo {sLam}. */

bool_t multifok_scene_object_XY_is_inside
  ( multifok_scene_object_t *obj,
    interval_t dom[],
    double margin
  );
  /* Checks whether the {XY} projection of object {obj} 
    is inside the rectangle {dom[0] × dom[1]}.
    
    If {margin} is negative, checks whether the {X} and {Y} coordinates
    of the center of the object )only) lie inside that rectangle.
    
    If {margin} is non-negative, checks whether the entire {XY} bounding
    box {obj->box[0] × obj->box[1]} is inside the rectangle {dom[0] ×
    dom[1]} and at least {margin} from its border. */

bool_t multifok_scene_object_XY_overlap
  ( multifok_scene_object_t *obja,
    multifok_scene_object_t *objb,
    double minSep
  );
  /* The {minSep} parameter must be non-negative. returns {TRUE} iff
    the {XY} projections of the objects {obja} and {objb}
    overlap, or are less than {minSep} of each other.
    
    As a special case, a {ot_RAMP} object is considered to overlap 
    any other object, and a {ot_FLAT} object is considered to overlap
    another object if {Z} range intersects the object's {Z} range.
    Otherwise, the {Z} ranges {obja->box[2]} and {objb->box[2]}
    are ignored. */
    
multifok_scene_object_t multifok_scene_object_background_make
  ( multifok_scene_object_type_t type,
    interval_t dom[],
    bool_t verbose
  );
  /* Generates the background object of the given {type} (either {ot_FLAT} or {ot_RAMP}).
    The {id} field is set to {multifok_scene_object_ID_NONE}.
    The object will have {bgLam} set to white and {fgLam} set to black,
    with zero {fgGlo,bgGlo}. */

multifok_scene_object_t multifok_scene_object_foreground_make
  ( multifok_scene_object_type_t type,
    interval_t bbox[],
    frgb_t *fgGlo,
    frgb_t *bgGlo,
    frgb_t *fgLam,
    frgb_t *bgLam,
    bool_t verbose
  );
  /*  Generates a single object of the specified foreground type
    (not {ot_RAMP} or {ot_FLAT}).
    
    The object will have the given bounding box and finish parameters
    {bgGlo,fgGlo,bgLam,fgLam}. The bounding box may be reduced so that
    it has the proper aspect ratio for the {type}. If any is {NULL},
    {0,0,0} is assumed. The {ID} will be set to
    {multifok_scene_object_ID_NONE}. */

multifok_scene_object_t multifok_scene_object_foreground_throw
  ( interval_t dom[],
    double margin,
    double wMinXY,
    double wMaxXY,
    frgb_t *fgGlo,
    bool_t verbose
  );
  /* Generates a random object of random foreground type (not {ot_RAMP} or {ot_FLAT})
  
    The {Z} range of the object {obj.bbox[2]} will be contained in the
    interval {dom[2]}. If {margin} is non-negative, the {XY} range
    {obj.bbox[0] × obj.bbox[1]} will be within the rectangle
    {dom[0]×dom[1]} and at least {margin} away from its four sides. If {margin}
    is negative, only the {XY} center of the object will be in that
    rectangle, but the object's {XY} projection may extend partially
    outside it.
   
    The width of {obj.bbox[0]} (X) and {obj.bbox[1]} (Y) will be
    randomly chosen in {[wMinXY _ wMaxXY]}, subject to the aspect ratio of
    the object type and to the object fitting in {dom} as specified
    above.
    
    The procedure fails if {dom} is too small for the chosen object
    type, given {wMinXY} and {margin}.
    
    The object will have random contrasting {bgLam} and {fgLam} colors.
    The gloss component {fgGlo} will be as given, and {bgGlo} will be
    {(0,0,0)}. The {ID} will be set to
    {multifok_scene_object_ID_NONE}. */

char *multifok_scene_object_type_to_string(multifok_scene_object_type_t type);
  /* Returns a (non-heap) string version of the type: "FLAT" for 
    {multifok_scene_object_type_FLAT}, "DISK" for 
    {multifok_scene_object_type_DISK}, etc. */

multifok_scene_object_type_t multifok_scene_object_type_from_string(char *str);
  /* Returns the {multifok_scene_object_type_t} value corresponding 
    to the string {str}, which may be "FLAT" for
    {multifok_scene_object_type_FLAT}, "DISK" for 
    {multifok_scene_object_type_DISK}, etc.  Fails
    if {str} is none of those. */

#endif
