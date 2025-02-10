/* Computing normals of objects. */
/* Last edited on 2025-02-02 19:31:55 by stolfi */

#ifndef multifok_scene_object_normal_H
#define multifok_scene_object_normal_H

#include <stdint.h>

#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <interval.h>

#include <multifok_scene.h>
#include <multifok_scene_object.h>

/* GENERIC PROCEDURES */

r3_t multifok_scene_object_normal(multifok_scene_object_t *obj, r3_t *q, bool_t verbose);
  /* Computes the outwards unit normal to the object {obj} at the point {q},
    assumed to be on its surface (apart from roundoff errors). */

/* TYPE-SPECIFIC PROCEDURES */   
   
r3_t multifok_scene_object_normal_FLAT(interval_t bbox[], r3_t *q, bool_t verbose);   
  /* Like {multifok_scene_object_normal_object}, but specialized 
    for an {ot_FLAT} floor object. Returns the ray's {t} parameter 
    if hit, {+INF} otherwise. */
      
r3_t multifok_scene_object_normal_RAMP(interval_t bbox[], r3_t *q, bool_t verbose);   
  /* Like {multifok_scene_object_normal_object}, but specialized for an {ot_RAMP}
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
r3_t multifok_scene_object_normal_DISK(interval_t bbox[], r3_t *q, bool_t verbose);   
  /* Like {multifok_scene_object_normal_object}, but specialized for an {ot_DISK}
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
r3_t multifok_scene_object_normal_BALL(interval_t bbox[], r3_t *q, bool_t verbose);   
  /* Like {multifok_scene_object_normal_object}, but specialized for an {ot_BALL} 
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
r3_t multifok_scene_object_normal_CONE(interval_t bbox[], r3_t *q, bool_t verbose);   
  /* Like {multifok_scene_object_normal_object}, but specialized for an {ot_CONE} 
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
r3_t multifok_scene_object_normal_PYRA(interval_t bbox[], r3_t *q, bool_t verbose);   
  /* Like {multifok_scene_object_normal_object}, but specialized for an {ot_PYRA} 
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */

#endif
