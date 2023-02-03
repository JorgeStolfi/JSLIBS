/* Test scenes and ray-tracing for {multifok_test}. */
/* Last edited on 2023-01-31 13:24:24 by stolfi */

#ifndef multifok_scene_H
#define multifok_scene_H

#define _GNU_SOURCE
#include <stdint.h>

#include <interval.h>
#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <frgb.h>
#include <float_image.h>

typedef struct multifok_scene_object_t
  { int32_t ID;         /* Object ID, used e.g. for patterning. */
    r3_t ctr;           /* Center of disk/sphere. */
    double rad;         /* Radius of disk/sphere. */
    bool_t flat;        /* True if disk, false if sphere. */
    interval_t bbox[3]; /* Bounding box for the object. */
    frgb_t bg;          /* "Background" color of object texture. */
    frgb_t fg;          /* "Foreground" color of object texture. */
  } multifok_scene_object_t;
  /* An object {obj} of this type pecifies a sphere (if {obj.flat} is false)
    or horizontal disk (if {obj.flat} is true), with center {ctr} and radius {rad}.
    Its color at each point {q} will be a varying combination of colors 
    {obj.bg} and {obj.fg}, depending on the coordinates of {q} and its {ID}. */ 

#define multifok_scene_ZMAX 30.0
  /* Total depth of simulated scene. */

struct multifok_scene_tree_t;
  /* A node in a tree of objets to speed up raytracing. */

typedef struct multifok_scene_t
  { interval_t dom[3];
    /* Disks and balls: */
    int32_t NO;
    multifok_scene_object_t *objs;
    struct multifok_scene_tree_t *tree; /* Raytracing speedup tree. */
    bool_t flatFloor; /* True if the floor is horizontal plane, false if relief: */
  } multifok_scene_t;
  /* An object {scene} of this type specifies a set of solid objects that 
    can be ray-traced.  It has a "floor" surface and some number {NO} of 
    objects {objs[0..NO-1]}.
    
    The {dom} box is the interesting part of the scene. Usually its {XY}
    projection {D = dom[0]×dom[1]} is the part that is projected onto
    the image in vertical cylindrical projection. 

    If {flatFloor} is true, the "floor" will be a horizontal plane at {Z = LO(dom[2])}. 
    
    If flatFloor} is false, the "floor" will be some surface that spans the whole range of {Z} coordinates {dom[2]} within
    the rectangle {D}. This shape is extended outside {D} with flat plane(s)
    as needed to ensure that its {Z} coordinates, even
    outside {D}, are always contained in {dom[2]}.
    
    The horizontal projection of some objects may extend outside {D}, but their {Z} ranges will be
    contained in {dom[2]}; which in turn must be contained in {[0_ZMAX]} where
    {ZMAX=multifok_scene_ZMAX}. */ 

multifok_scene_t *multifok_scene_new(interval_t dom[], bool_t flatFloor, bool_t verbose);
  /* Creates a {scene} with the given {dom} box, consisting of a floor and no objects.  
    The meaning of {flatFloor} is specified under {multifok_scene_t}.
    The scene will have no objects ({NO=0}) */

void multifok_scene_throw_objects
  ( multifok_scene_t *scene,
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

multifok_scene_object_t multifok_scene_object_throw(int32_t ID, interval_t dom[], double margin, double rMin, double rMax);
  /* Generates a random sphere or disk with the given {ID} and center in the given
    3D domain {dom[0..2]}.  The object will have a random center and radius, and random
    but contrasting {bg} and {fg} colors.
    
    The radius {rad} of the object will in principle be chosen randomly
    in {[rMin _ rMax]}.
    
    In any case the object's {Z} projection will be strictly inside
    {dom[2]}. If the given {margin} is non-negative, the {XY} projection
    of the object coordinates will fit entirely within the rectangle
    {dom[0]×dom[1]} and at least {margin} away from its four sides. The
    object's radius may have to be reduced, even below {rMin}, if
    necessary to satisfy these containment constraints. If {margin} is
    negative, the center's {X} and {Y} will be inside that rectabgle,
    but the object's {XY} projection may extend partially ouside it.
    
    In any case, the {dom} must be large enough to make the placement possible. */

void multifok_scene_ray_trace
  ( multifok_scene_t *scene,
    r3_t *p, 
    r3_t *d, 
    bool_t debug,
    multifok_scene_object_t **hob_P, 
    r3_t *hit_pos_P
  );
  /* Traces one ray {R} that goes through the point {p} with the
    direction parallel to {d}, assumed to be not horizontal. Finds the
    object (disk, sphere, or floor) with max {Z} that hits that ray.
    Returns that object in {*hob_P} and the hit point
    {hpt(R)} in {*hpt_P}.  If the 
    
    Uses the tree structure {tr} to speed up the computation. */
    
/* OBJET TEXTURING */

typedef double multifok_scene_pattern_t(double x, double y, double z, int32_t iobj);
  /* Type of a function that computes a grayscale value as a function of 
    an arbitrary /object index/ {iobj} and the point {(x,y,z)}. */
    
frgb_t multifok_scene_compute_hit_color(multifok_scene_object_t *obj, r3_t *q, multifok_scene_pattern_t *pattern);
  /* Computes the color {clr} of the surface of object {obj} at the point {hpt}.
    
    If {obj} is not {NULL}, the procedure evaluates
    {r=pattern(x,y,z,obj.ID)}, where {(x,y,z) = q - obj.ctr}, and
    obtains {clr} by interpolating between {obj.bg} and {obj.fg} with
    the ratio {r}. If {obj} is NULL, computes instead {r =
    pattern(x,y,z,-1)} where {(x,y,z) = q}, and maps {r} lineary to a
    grayscale value from back to white. */
 
#endif
