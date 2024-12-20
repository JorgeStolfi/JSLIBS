/* Ray-tracing scene objects for multi-focus stereo. */
/* Last edited on 2024-12-05 10:36:43 by stolfi */

#ifndef multifok_scene_object_raytrace_H
#define multifok_scene_object_raytrace_H

#include <stdint.h>

#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <interval.h>
#include <frgb.h>

#include <multifok_image.h>
#include <multifok_scene.h>
#include <multifok_scene_object.h>

double multifok_scene_object_raytrace
  ( multifok_scene_object_t *obj,
    r3_t *p,
    r3_t *d,
    double tMin,
    double tMax,
    bool_t verbose
  );
  /* Checks whether the ray {R} defined by point {p} and unit direction vector {d}
    intersects the upper surface of the object {obj} at 
    a parameter value {tHit} in the range {[tMin _ tMax]}. If it does,
    returns {tHit}, otherwise returns {+INF}. */

/* TRACING SPECIFI OBJECTS */
    
double multifok_scene_object_raytrace_FLAT(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_object_raytrace}, but specialized 
    for a {FLAT} floor object. Returns the ray's {t} parameter 
    if hit, {+INF} otherwise. */
      
double multifok_scene_object_raytrace_RAMP(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_object_raytrace}, but specialized for a {RAMP}
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
double multifok_scene_object_raytrace_DISK(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_object_raytrace}, but specialized for a {DISK}
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
double multifok_scene_object_raytrace_BALL(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_object_raytrace}, but specialized for a {BALL} 
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
double multifok_scene_object_raytrace_CONE(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_object_raytrace}, but specialized for a {CONE} 
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */

    
#endif
