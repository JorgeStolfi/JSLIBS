/* Ray-tracing test scenes for multi-focus stereo. */
/* Last edited on 2025-01-30 14:59:53 by stolfi */

#ifndef multifok_scene_raytrace_H
#define multifok_scene_raytrace_H

#include <stdint.h>

#include <r3.h>
#include <hr3.h>
#include <bool.h>
#include <interval.h>
#include <frgb.h>

#include <multifok_image.h>
#include <multifok_scene.h>
#include <multifok_scene_object.h>
#include <multifok_scene_tree.h>

void multifok_scene_raytrace
  ( multifok_scene_t *scene,
    multifok_scene_tree_t *tree,
    r3_t *p, 
    r3_t *d, 
    bool_t verbose,
    multifok_scene_object_t **oHit_P, 
    r3_t *pHit_P
  );
  /* Traces one ray {R} that goes through the point {p} with the
    direction parallel to {d}, assumed to be not horizontal. 
    
    Finds the first object (disk, sphere, floor, etc.) hit by that ray.
    Returns that object in {*oHit_P} and the scene coordinates of the
    hit point {pHit(R)} in {*pHit_P}.
    
    If the ray desn't hit any object, returns {NULL} in {*oHit_P}.
    The value of {*pHit_P} will be undefined.
    
    Uses the tree structure {tree} to speed up the computation.*/

void multifok_scene_raytrace_get_ray_bbox
  ( r3_t *p,
    r3_t *d,
    double tMin,
    double tMax,
    interval_t bbox[]
  );
  /* Sets {bbox[0..2]} to the bounding box of the part of the ray
    with parameter value in the range {[tMin _ tMax]}. */

/* SCENE TO IMAGE COORDINATE MAPPING */

r3_t multifok_scene_coords_to_image_coords
  ( r3_t *p3_scene,
    double zFoc,
    interval_t scene_dom[],
    int32_t NX,
    int32_t NY
  );
  /* Given the scene coordinates {p3_scene} of a point,
    returns the corresponding 3D image coordinates {p3_img}.
    
    The mapping will be such that the scene's {XY} domain {scene_dom[0]
    × scene_dom[0]} fits snugly and concentrically in the {XY} image
    domain {[0 _ NX] × [0 _ NY]}, with the same scene-to-image unit
    scale factor in both axes.
    
    The {Z} coordinate of {p3_img} will be the height, in pixels, of {p3_scene} 
    above the the in-focus plane; that is, the {Z} coordinate of {p3_scene}
    minus the nominal position {zFoc} of taht plane, scaled by the same factor. */

r3_t multifok_scene_coords_from_image_coords
  ( r3_t *p3_img,
    int32_t NX,
    int32_t NY,
    interval_t scene_dom[],
    double zFoc
  );
  /* Given the image coordinates {p3_img} of a point in 
    image space, returns the scene coordinates {p_scene)of
    the corresponding point in scene space.  
    It is the inverse of {multifok_scene_coords_to_image_coords},
    if given the same parameters{NX,NY,scene_adom[0..1],zFoc}. */
 
double multifok_scene_pixel_size(int32_t NX, int32_t NY, interval_t scene_dom[]);
  /* Returns the pixel size in scene coordinates. Fails if the scale
    of scene to image coordinates is not the same on both axes. */

#endif
