/* Test scenes and ray-tracing for {multifok_test}. */
/* Last edited on 2024-10-29 13:18:45 by stolfi */

#ifndef multifok_scene_H
#define multifok_scene_H

#define _GNU_SOURCE
#include <stdint.h>

#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <interval.h>
#include <frgb.h>

#include <multifok_image.h>
#include <multifok_frame.h>
#include <multifok_scene_object.h>

typedef struct multifok_scene_t
  { interval_t dom[3];  /* Clips and contains the scene. */ 
    int32_t NO;         /* Number of objects (including floor or ramp). */
    multifok_scene_object_t *objs;
  } multifok_scene_t;
  /* An object {scene} of this type specifies a set of solid objects
    that can be ray-traced. It has some number {NO} of objects
    {objs[0..NO-1]}.
    
    The {dom} box is the interesting part of the scene.
    
    If flatFloor} is false, the "floor" will be some surface that spans
    the whole range of {Z} coordinates {dom[2]} within the rectangle
    {D}. This shape is extended outside {D} with flat plane(s) as needed
    to ensure that its {Z} coordinates, even outside {D}, are always
    contained in {dom[2]}.
    
    The horizontal projection of some objects may extend outside {D},
    but their {Z} ranges will be contained in {dom[2]}. */ 

multifok_scene_t *multifok_scene_new(interval_t dom[], bool_t verbose);
  /* Creates a {scene} with the given {dom} box, consisting of a floor and no objects.
    The scene will have no objects ({scene.NO=0, scene.objs=NULL}) */

void multifok_scene_throw_objects
  ( multifok_scene_t *scene,
    bool_t floorOnly,
    bool_t flatFloor,
    double rMin, 
    double rMax,
    double minSep, 
    bool_t verbose
  );
  /* Adds to the given {scene} a set of objects picked at random, with centers inside
    the 3D domain {dom = scene.dom}.  
    
    The number of objects is chosen by the procedure. The objects will
    be generated with {multifok_scene_object_throw(dom,rMin,rMax)}.
    
    If {minSep} is negative, the objects may overlap in {X} and {Y} and
    may extend partially outside the rectangle {dom[0]×dom[1]}. If
    {minSep} is zero or positive, the {X} and {Y} projections of the
    objects will be disjoint and separated by at least {minSep} pixels,
    and will be at least {minsep} away from the edges of that rectangle.
    
    In any case, the {Z} coordinates of all objects will be contained in 
    {dom[2]}. 
    
    The procedure also builds a raytrace speedup tree for those objects. */ 
    
r3_t multifok_scene_box_center(interval_t box[]);
  /* Returns the center of the box {box[0] × box[1] × box[2]}. */
  
r3_t multifok_scene_box_radius(interval_t box[]);
  /* Returns a vector {rad} such that {rad.c[j]} is the half-size 
    of {box[j]}. */
  
void multifok_scene_print_box(FILE *wr, char *pref, interval_t box[], char *suff);
  /* Prints {box[0..2]} to {wr}, preceded by {pref} and followed by {suff}. */
    
void multifok_scene_check_object_IDs(int32_t NO, multifok_scene_object_t objs[]);
  /* Checks whether the {ID} fields of {objs[0..NO-1]} are a permutation
    of {0..NO-1}. */

#endif
