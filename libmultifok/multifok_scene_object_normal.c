/* See {multifok_scene_object_normal.h}. */
/* Last edited on 2024-10-29 19:09:04 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <r3.h>
#include <frgb.h>
#include <r2.h>
#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsqroots.h>

#include <multifok_scene.h>
#include <multifok_scene_object.h>

#include <multifok_scene_object_normal.h>

#define ot_FLAT multifok_scene_object_type_FLAT
#define ot_RAMP multifok_scene_object_type_RAMP
#define ot_DISK multifok_scene_object_type_DISK
#define ot_BALL multifok_scene_object_type_BALL
#define ot_CONE multifok_scene_object_type_CONE

#define FUDGE (1.0e-6)
  /* Fudge amount to expand bounding boxes to account for roundoff. */

r3_t multifok_scene_object_normal(multifok_scene_object_t *obj, r3_t *q, bool_t verbose)
 {
    switch(obj->type)
      { case ot_FLAT: 
          return multifok_scene_object_normal_FLAT(obj->bbox, q, verbose); break;
        case ot_RAMP:
          return multifok_scene_object_normal_RAMP(obj->bbox, q, verbose); break;
        case ot_DISK:
          return multifok_scene_object_normal_DISK(obj->bbox, q, verbose); break;
        case ot_BALL:
          return multifok_scene_object_normal_BALL(obj->bbox, q, verbose); break;
        case ot_CONE:
          return multifok_scene_object_normal_CONE(obj->bbox, q, verbose); break;
        default: demand(FALSE, "unrecognized object type");
      }
 }
 
r3_t multifok_scene_object_normal_FLAT(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    return (r3_t){{ 0.0, 0.0, 1.0 }};
  }
      
r3_t multifok_scene_object_normal_RAMP(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    /* Grab the {X} range of the plane strip:  */
    /* It should be a bit less than the scene's {X} range: */
    double bXrad = interval_rad(&(bbox[0]));
    double bXlo = bbox[0].end[0] + WD_RAMP*bXrad;
    double bXhi = bbox[0].end[1] - WD_RAMP*bXrad;
    
    if ((q->c[0] <= bXlo) || (q->c[0] >= bXhi))
      { return (r3_t){{ 0.0, 0.0, 1.0 }}; }
    else
      { /* Grab the {Z} range of the plane strip: */
        double bZlo = bbox[2].end[0] + FUDGE; /* In case of roundoff errors. */
        double bZhi = bbox[2].end[1] - FUDGE; /* In case of roundoff errors. */
        
        r3_t nrm = (r3_t){{ -(bZhi-bZlo), 0.0, +(bXhi-bXlo) }};
        (void)r3_dir(&nrm, &nrm);
        return nrm;
      }
  }
      
r3_t multifok_scene_object_normal_DISK(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    return (r3_t){{ 0.0, 0.0, 1.0 }};
  }
      
r3_t multifok_scene_object_normal_BALL(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    r3_t ctr = multifok_scene_box_center(bbox);
    r3_t rad; r3_sub(q, &ctr, &rad); (void)r3_dir(&rad, &rad);
    return rad;
  }
      
r3_t multifok_scene_object_normal_CONE(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    r3_t ctr = multifok_scene_box_center(bbox);
    r3_t tip = (r3_t){{ ctr.c[0], ctr.c[1], bbox[2].end[1] }};
    r3_t rad; r3_sub(q, &tip, &rad);  /* Vector from cone tip to {q}. */
    double h = hypot(rad.c[0], rad.c[1]);
    if (h < 1.0e-12)
      { return (r3_t){{ 0.0, 0.0, 1.0 }}; }
    else
      { r3_t nrm = (r3_t){{ -rad.c[2]*rad.c[0]/h, -rad.c[2]*rad.c[1]/h,  h }};
        assert(fabs(r3_dot(&rad, &nrm)) < 1.0e-8);
        (void)r3_dir(&nrm, &nrm);
        return nrm;
      }
  }


