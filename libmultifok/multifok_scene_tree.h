/* A tree structure to speed up raytracing of a scene. */
/* Last edited on 2024-10-25 08:48:13 by stolfi */

#ifndef multifok_scene_tree_H
#define multifok_scene_tree_H

#define _GNU_SOURCE
#include <stdint.h>

#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <interval.h>
#include <frgb.h>

#include <multifok_image.h>
#include <multifok_frame.h>
#include <multifok_scene.h>
#include <multifok_scene_object.h>

typedef struct multifok_scene_tree_t
  { multifok_scene_object_t *obj; /* The tree's root object. */
    interval_t bbox[3];           /* Bounding box of all objects. */
    int8_t axis;                  /* Axis perp to split plane, 0 or 1. */
    struct multifok_scene_tree_t *sub[2]; /* Subtrees, or NULL. */
  } multifok_scene_tree_t;
  /* Node in an acceleration tree. 
  
    The node stores a pointer to one object {obj}, and to two subtrees
    {sub[0,1]}. The two sets of objects are disjoint and {obj} is in
    neither of them. The {bbox}s of the subtrees may overlap, and
    overlap with {obj}; but the objects in {sub[0]} tend to have lower
    coordinates along the specified {axis} than those in {sub[1]}.
    
    If both subtrees are {NULL} then {axis} is irrelevant. */

multifok_scene_tree_t *multifok_scene_tree_build
  ( int32_t NO,
    multifok_scene_object_t objs[],
    int32_t debug_level
  );
  /* Builds an acceleration tree for the objects {objs[0..NO-1]}.  Will rearrange the
    objects in the process.
    
    If {debug_level} is non-negative, also prints debugging information. */
 
#endif
