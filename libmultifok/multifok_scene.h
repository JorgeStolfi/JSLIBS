/* Test scenes and ray-tracing for {multifok_test}. */
/* Last edited on 2023-01-25 20:05:12 by stolfi */

#ifndef multifok_scene_H
#define multifok_scene_H

#define _GNU_SOURCE
#include <stdint.h>

#include <interval.h>
#include <r3.h>
#include <bool.h>
#include <frgb.h>
#include <float_image.h>

typedef struct multifok_scene_object_t
  { r3_t ctr;          /* Center of disk/sphere. */
    double rad;        /* Radius of disk/sphere. */
    bool_t flat;       /* True if disk, false if sphere. */
    frgb_t bg;         /* Background color of object pattern. */
    frgb_t fg;         /* Foreground color of object pattern. */
  } multifok_scene_object_t;
  /* Specifies a sphere or horizontal disk with center {ctr}, radius {rad}, and color {bg[]}. */ 

typedef struct multifok_scene_t
  { interval_t box[3];
    /* Disks and balls: */
    int32_t ND;
    multifok_scene_object_t *obj;
    /* Backplane tilt and colors: */
    bool_t back_tilt;
    frgb_t bg;
    frgb_t fg;
  } multifok_scene_t;
  /* A test scene to generate test images from. It has {ND} objects
    entirely contained in the given {box}. The lowest {Z} of the {box}
    must be non-negative.
    
    It also has a backplane with colors ranging between {bg} and {fg}.
    If {back_tilt} is false, the backplane is horizontal at {Z=0}.
    If {back_tilt} is true, the backplane contains the bottom left
    and top right edges of the {box} .*/ 

typedef void multifok_scene_pattern_t(double x, double y, double z, int32_t iobj, int32_t NC, float fs[]);
  /* Type of a function that computes pixel samples {fs[0..NC-1]} as a function of 
    an arbitrary /object index/ {iobj} and {(x,y,z)}. */

#define multifok_scene_ZMAX 30.0
  /* Total depth of simulated scene. */

multifok_scene_t *multifok_scene_throw
  ( interval_t box[],
    bool_t back_tilt,
    double rMin, 
    double rMax,
    double minSep, 
    bool_t verbose
  );
  /* Generates a test scene with several objects picked at random inside the box {box[0..2]}. 
    
    The backplane will be horizontal at {Z=box[2].end[0]} if {back_tilt} is false,
    or tilted across the {box} if {back_tilt} is true (see the
    {back_tilt} attribute of {multifok_scene_t}). In the second
    case there will be no objects, and the parameters {rMin,rMax,minSep}
    are ignored. The backplane's {bg} and {fg} colors are chosen at
    random.
    
    If {back_tilt} is false, the number of objects is chosen by the
    procedure. The objects will be generated with
    {multifok_scene_object_throw(box,rMin,rMax)}.
    
    If {minSep} is negative, the objects may overlap in {X} and {Y} and
    may extend partially outside the rectangle {box[0]×box[1]}. If
    {minSep} is zero or positive, the {X} and {Y} projections of the
    objects will be disjoint and separated by at least {minSep} pixels,
    and will be at least {minsep} away from the edges of that rectangle.
    
    In any case, the {Z} range of the objects will be stricty contained
    in the interval {box[2]}, which should be contained in {[0_ZMAX]}
    where {ZMAX=multifok_scene_ZMAX}. */ 

multifok_scene_object_t multifok_scene_object_throw(interval_t box[], double margin, double rMin, double rMax);
  /* Generates a random sphere or disk entirely contained in the given
    3D domain {box[0..2]}.  The object will have a random center and radius, and random
    but contrasting {bg} and {fg} colors.
    
    The radius {rad} of the object will in principle be chosen randomly in {[rMin _ rMax]}. 
    
    in any case the object's {Z} projection will be strictly inside {box[2]}.  
    If the given {margin} is non-negative, the {XY} projection of the object 
    coordinates will fit entirely within the rectangle {box[0]×box[1]} and
    at least {margin} away from its four sides.  If {margin} is negative, the center's {X} and {Y}
    will be inside that rectabgle, but the object's {XY} projection may extend
    partially ouside it.  The object's radius may have to be reduced, even below {rMin},
    if necessary to satisfy these containment constraints.
    
    In any case, the {box} must be large enough to make the placement possible. */

void multifok_scene_ray_trace
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    r3_t *p, 
    r3_t *d, 
    bool_t debug,
    frgb_t *ray_clr_P, 
    r3_t *ray_hit_P
  );
  /* Traces one ray {R} that goes through the point {p} with the
    direction parallel to {d}, assumed to be not horizontal. Finds the
    object (disk, spher, or backplane) with max {Z} that hits that ray.
    Returns its color {clr(R)} in {*ray_clr_P} and the hit point
    {hit(R)} in {*ray_hit_P}.
    
    The color {clr(R)} is the pixel {fs[0..2]} returned by
    {pattern(x,y,z,kd,3,fs)}. If {hit(R)} is on a disk or sphere with
    center {ctr}, {kd} is the index if the object that was hit, the
    pattern coordinates {(x,y,z)} are {hit(R)-ctr}. If {hit(R)} is on
    the backplane, {kd} is set to {-1}, {x} and {y} are the coords of
    {hit(R)}, and {z} is set to 0, ignoring the actual {Z} of the
    backplane. */
 
#endif
