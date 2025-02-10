/* Test scenes and ray-tracing for {multifok_test}. */
/* Last edited on 2025-02-10 03:32:55 by stolfi */

#ifndef multifok_scene_H
#define multifok_scene_H

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
  { interval_t dom[3];                 /* Clips and contains the scene. */ 
    uint32_t NO;                       /* Number of objects (incl background). */
    multifok_scene_object_t *objs;     /* The objects. */
  } multifok_scene_t;
  /* An object {scene} of this type specifies a set of solid objects
    that can be ray-traced. It has some number {NO} of objects
    {objs[0..NO-1]}.  In a completed scene, there must be exactly 
    one /floor/ object {objs[0]} with type {ot_FLAT} or {ot_RAMP},
    and zero or more /foreground/ objects of any types except those 
    two.
    
    The {dom} box is the interesting part of the scene.
    
    The {X} and/or {Y} projection of some objects may extend outside
    {dom[0]} and {dom[1]}, but their {Z} ranges will be contained in
    {dom[2]}. */ 

multifok_scene_t *multifok_scene_new(interval_t dom[], bool_t verbose);
  /* Creates a {scene} with the given {dom} box and {pattern} function.
    The scene will have no objects ({scene.NO=0, scene.objs=NULL}) */

void multifok_scene_add_floor
  ( multifok_scene_t *scene,
    multifok_scene_object_type_t type, 
    bool_t verbose
  );
  /* Adds to the given {scene} a floor object of the specified type. The
    scene must currently have no objects.  */ 
    
void multifok_scene_add_foreground_object
  ( multifok_scene_t *scene,
    multifok_scene_object_type_t type,
    interval_t bbox[],
    frgb_t *fgGlo,
    frgb_t *bgGlo,
    frgb_t *fgLam,
    frgb_t *bgLam,
    bool_t verbose
  );
  /* Adds to the given {scene} a single foreground object of type {type}. 
    
    The scene must currently have at least the floor object. The new object
    {type} must not be {ot_FLAT} or {ot_RAMP}.

    The object is generated with {multifok_scene_object_foreground_make}
    (q.v.) with parameters {type,dom,bbox,fgGlo,bgGlo,fgLam,bgLam)}.
    
    The center of {bbox[0]} and {bbox[1]}, as well as the entire range
    {bbox[2]}, must be inside {scene.dom}. */ 

void multifok_scene_throw_foreground_objects
  ( multifok_scene_t *scene,
    double wMinXY, 
    double wMaxXY,
    frgb_t *fgGlo,
    double minSep, 
    bool_t verbose
  );
  /* Adds to the given {scene} a set of objects picked at random, with
    centers inside the 3D domain {dom = scene.dom}. The scene must
    currently have only the floor object.
    
    The number of objects is chosen by the procedure. The objects will
    be generated with 
    {multifok_scene_object_foreground_throw(dom,wMinXY,wMaxXY,fGlo,verbose)}
    (q. v.).
    
    If {minSep} is negative, the objects may overlap in {X} and {Y} and
    may extend partially outside the rectangle {dom[0]×dom[1]}. If
    {minSep} is zero or positive, the {X} and {Y} projections of the
    objects will be disjoint and separated by at least {minSep} pixels,
    and will be at least {minsep} away from the edges of that rectangle.
    
    In any case, the {Z} coordinates of all objects will be contained in 
    {dom[2]}. */ 
    
void multifok_scene_print(FILE *wr, multifok_scene_t *scene);
  /* Prints the scene to {wr}. */
  
void multifok_scene_print_box(FILE *wr, char *pref, interval_t box[], char *suff);
  /* Prints {box[0..2]} to {wr}, preceded by {pref} and followed by {suff}. */
    
void multifok_scene_check_object_IDs(uint32_t NO, multifok_scene_object_t objs[]);
  /* Checks whether the {ID} fields of {objs[0..NO-1]} are a permutation
    of {0..NO-1}. */

#endif
