/* Creates focus-blurred images of a scene. */
/* Last edited on 2024-10-24 03:37:56 by stolfi */

#ifndef multifok_scene_make_frame_H
#define multifok_scene_make_frame_H

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
#include <multifok_raytrace.h>
#include <multifok_scene_raytrace.h>

multifok_frame_t *multifok_scene_make_frame
  ( int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_tree_t *tree,
    multifok_scene_raytrace_pattern_t *pattern,
    double zFoc, 
    double zDep, 
    int32_t HS,
    int32_t KR,
    bool_t verbose,
    multifok_raytrace_debug_pred_t *debug_pix,
    multifok_raytrace_report_ray_proc_t *report_ray
  );
  /* Creates a frame record {fr} with images {fr.sVal}, {fr.shrp}, {fr.hAvg},
    and {fr.hDev} of the given {scene}, with simulated depth-of-focus
    blur. All four images will have {NX} columns of {NY} pixels.
    
    The scene will be implicitly scaled so that the entire image fits in
    its {XY} domain.
    
    The {fr.sVal} image will have 3 channels, and shows the color of the
    scene as seen through a camera or muscope focused at {Z=zFoc} with
    focal depth {zDep}.
    
    The {fr.shrp} image will have 1 channel, and shows the estimated or
    computed sharpness-of-focus indicator at each pixel, a non-negative
    number that is inversely proportional to the square of the spread
    radius of the blurring. Since the latter is at least 1.0 (because of
    pixel antialiasing by sub-sampling) the sharpness is between 1 and
    0, where 1 would indicate a maximally sharp image.
    
    The {fr.hAvg} image will have 1 channel, and shows the average {Z}
    coordinate of the scene at each pixel. Its values will be contained in
    {[zMin _ zMax] = scene.dom[2]} which in turn will be contained in 
    {[multifok_scene_ZMIN _ multifok_scene_ZMIN]}.

    The {fr.hDev} image will have 1 channel, and shows the standard
    deviation of the {Z} coordinate of the scene inside each pixel. It
    should be in the range {[0 _ zdMax]} where {zdMax} is
    {max(|zMin-zFoc|, |zMax-zFoc|)}.

    The value of each pixel of the four images is computed by
    {multifok_raytrace}, using an array of {2*HS+1} by {2*HS+1} sampling
    points for each pixel, and some number {NR} of rays for each sampling
    point {p}.
    
    If {KR} is zero (which must be the case if {zDep=+INF}), the same
    vertical ray direction will be used for all sampling points.
    Otherwise a distinct set of {KR} rays will be used for each sampling
    point. In this case, if there is more than one sampling point
    ({HS>0}) the weighted average RMS deviation of the rays from the
    vertical, over a vertical travel of {zDep} away from the focus
    plane, will be about 1.
    
    The color of the scene at a ray hit point {x,y,z} is obtained by
    computing {r = pattern(x,y,z,kd,3,fs)} where and {kd} is the index
    of the object hit by the ray, or {-1} if the ray hit the background
    surface.
    
    if {verbose} is true, will print to {stderr} some general debugging
    information about the frame building.
    
    The predicate {debug_pix}, if not {NULL}, is called for each pixel
    of the frame. It should return true if detailed debugging
    information is to be printed to {stderr} about that pixel. In that
    case, is {report_ray} is not {NULL} the procedure will also call
    {report_ray} for each ray traced on behalf of that pixel. */

#endif
