/* Objects for a {multifok_scene_t}. */
/* Last edited on 2024-10-24 13:09:17 by stolfi */

#ifndef multifok_scene_object_H
#define multifok_scene_object_H

#define _GNU_SOURCE
#include <stdint.h>

#include <r3.h>
#include <bool.h>
#include <interval.h>
#include <frgb.h>

typedef enum 
  { multifok_scene_object_type_FLAT,   /* Flat floor near bottom {Z}of domain. */
    multifok_scene_object_type_RAMP,   /* Ramp where {Z} increases with {X}. */
    multifok_scene_object_type_DISK,   /* A flat horizontal disk. */
    multifok_scene_object_type_BALL,   /* A ball. */
  } multifok_scene_object_type_t;
  /* Type of scene object. */
    
typedef struct multifok_scene_object_t
  { int32_t ID;                        /* Object ID, used e.g. for patterning. */
    multifok_scene_object_type_t type; /* Type of object. */
    interval_t bbox[3];                /* Bounding box for the object. */
    frgb_t bg;                         /* "Background" color of object texture. */
    frgb_t fg;                         /* "Foreground" color of object texture. */
  } multifok_scene_object_t;
  /* A record {obj} that specifies a an object in the scene.
    The {box} is a tight 3D bounding box for the object.
    
    If {type} is {BALL} the object is a ball inscribed in the specified
    {box}, which should be a cube.
    
    If {type} is {DISK} the object is a horizontal disk. Its {XY}
    projection is inscribed in the rectangle {box[0]×box[1]}, which
    should be a square, and its {Z} coordinate is the midpoint of
    {box[2]}, which should have very small height.
    
    If {type} is {FLAT}, the object is a horizontal rectangle whose {XY}
    projection is {box[0]×box[1]}, and whose {Z} coordinate is the
    midpoint of {box[2]}, which should have very small height.
        
    If {type} is {RAMP}, the object consists of three rectangles that
    span the whole range {box[1]} in the {Y} direction and one third of
    the range {box[0]} in the {X} direction. The middle rectangle is
    parallel to the {Y} axis and rises from {Zlo} when {X} is {Xlo} to
    {Zhi} when {X} is {Xhi}, where {Xlo,Xhi} are {1/3} and {2/3} of the
    way across {box[0]}, and {Zlo,Zhi} are the limits of {box[2]}. The
    first rectangle is horizontal at height {Zlo}, and the third one is
    horizontal at height {Zhi}.

    For all objects, {box[0..2]} is a box that encloses the object.
    
    The color of the object at each point {q} will be a varying
    combination of colors {obj.bg} and {obj.fg} in a way described 
    elsewhere. */ 

void multifok_scene_object_print(FILE *wr, multifok_scene_object_t *obj);
  /* Writes a readable description of object {obj} to {wr}, in one line. 
    Does NOT write the end-of-line. */

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
    
    As a special case, a {RAMP} object is considered to overlap 
    any other object, and a {FLAT} object is considered to overlap
    another object if {Z} range intersects the object's {Z} range.
    Otherwise, the {Z} ranges {obja->box[2]} and {objb->box[2]}
    are ignored. */
    
multifok_scene_object_t multifok_scene_object_background_make(interval_t dom[], bool_t flatFloor);
  /* Generates the background object, of type {FLAT} or {RAMP} depending on {flatFloor}.
    The {id} field is set to {-1}.
    
    The object will have {bg} set to white and {fg} set to black. */

multifok_scene_object_t multifok_scene_object_foreground_throw
  ( interval_t dom[],
    double margin,
    double rMin,
    double rMax,
    bool_t verbose
  );
  /* Generates a random object of random foreground type (not {RAMP} or {FLAT})
  
    The {Z} range of the object {obj.box[2]} will be contained in the
    interval {dom[2]}. If {margin} is non-negative, the {XY} range
    {obj.box[0] × obj.box[1]} too will be within the rectangle{dom[0] ×
    dom[1]} and at least {margin} away from its four sides. If {margin}
    is negative, only the {XY} center of the object will be in that
    rectangle, but the object's {XY} projection may extend partially
    outside it.
   
    The radius of the object (the max half-width of {obj.box[j]} for any
    {j}) will be randomly chosen in {[rMin _ rMax]}, subject to it
    fitting in {dom} as specified above. 
    
    The procedure fails if {dom} is
    too small for the chosen object type, given {rMin} and {margin}.
    
    The object will have random contrasting {bg} and {fg} colors.  The {ID}
    will be set to {-1}.  */

char *multifok_scene_object_type_to_string(multifok_scene_object_type_t type);
  /* Returns a (non-heap) string version of the type: "FLAT" for 
    {multifok_scene_object_type_FLAT}, "DISK" for 
    {multifok_scene_object_type_DISK}, etc. */

#endif
