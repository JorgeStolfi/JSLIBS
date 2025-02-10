/* See {multifok_scene_object_normal.h}. */
/* Last edited on 2025-02-07 06:46:45 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <box.h>
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
#define ot_PYRA multifok_scene_object_type_PYRA

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
        case ot_PYRA:
          return multifok_scene_object_normal_PYRA(obj->bbox, q, verbose); break;
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
        
        r3_t sNrm = (r3_t){{ -(bZhi-bZlo), 0.0, +(bXhi-bXlo) }};
        (void)r3_dir(&sNrm, &sNrm);
        return sNrm;
      }
  }
      
r3_t multifok_scene_object_normal_DISK(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    return (r3_t){{ 0.0, 0.0, 1.0 }};
  }
      
r3_t multifok_scene_object_normal_BALL(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    r3_t ctr; box_center(3, bbox, ctr.c);
    r3_t rad; r3_sub(q, &ctr, &rad); (void)r3_dir(&rad, &rad);
    return rad;
  }
      
r3_t multifok_scene_object_normal_CONE(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    r3_t ctr; box_center(3, bbox, ctr.c);
    r3_t tip = (r3_t){{ ctr.c[0], ctr.c[1], bbox[2].end[1] }};
    r3_t u; r3_sub(q, &tip, &u);  /* Vector from cone tip to {q}. */
    double h = hypot(u.c[0], u.c[1]);
    if (h < 1.0e-12)
      { return (r3_t){{ 0.0, 0.0, 1.0 }}; }
    else
      { r3_t sNrm = (r3_t){{ -u.c[2]*u.c[0]/h, -u.c[2]*u.c[1]/h,  h }};
        assert(fabs(r3_dot(&u, &sNrm)) < 1.0e-8);
        (void)r3_dir(&sNrm, &sNrm);
        return sNrm;
      }
  }
      
r3_t multifok_scene_object_normal_PYRA(interval_t bbox[], r3_t *q, bool_t verbose)
  {
    double cX, rX; interval_mid_rad(&(bbox[0]), &cX, &rX);
    double cY, rY; interval_mid_rad(&(bbox[1]), &cY, &rY);
    double loZ = bbox[2].end[0];
    double hiZ = bbox[2].end[1];

    demand(q->c[2] >= loZ - 1.0e-10, "hit {Z} below base");
    demand(q->c[2] <= hiZ + 1.0e-10, "hit {Z} above apex");

    assert(fabs(rX-rY) < 1.0e-12*(rX+rY));
    double rXY = 0.5*(rX + rY);
    double hZ = hiZ - loZ;

    r3_t a = (r3_t){{ cX, cY, hiZ }}; /* Pyramid tip. */
    r3_t qa; r3_sub(q, &a, &qa);  /* Vector from {a} to {q}. */

    demand(fmax(fabs(qa.c[0]), fabs(qa.c[1])) < rXY + 1.0e-10, "hit outside {XY} projection");

    /* Determine the hit face and its normal: */
    r3_t sNrm;
    if (q->c[2] <= loZ)
      { /* Must have hit base from below: */
        sNrm = (r3_t){{ 0, 0, -1 }};
      }
    else if (fabs(qa.c[0]) > fabs(qa.c[1]))
      { /* Hit a side with normal parallel to the {XZ} plane: */
        sNrm = (r3_t){{ (qa.c[0] < 0 ? -hZ : +hZ), 0.0, rXY }};
      }
    else
      { /* Hit a side with normal parallel to the {YZ} plane: */
        sNrm = (r3_t){{ 0.0, (qa.c[1] < 0 ? -hZ : +hZ), rXY }};
      }
    /* assert(fabs(r3_dot(&qa, &sNrm)) < 1.0e-8); */
    (void)r3_dir(&sNrm, &sNrm);
    return sNrm;
  }


