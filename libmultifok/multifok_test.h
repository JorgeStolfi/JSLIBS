/* Test tools for {multifok_focus_op} and related funcs. */
/* Last edited on 2023-01-10 20:43:59 by stolfi */

#ifndef multifok_test_H
#define multifok_test_H

#define _GNU_SOURCE
#include <stdint.h>

#include <interval.h>
#include <r3.h>
#include <bool.h>
#include <frgb.h>
#include <float_image.h>

typedef struct multifok_test_disk_t
  { r3_t ctr;          /* Center of disk. */
    double rad;        /* Radius of disk. */
    frgb_t bg;         /* Background color of disk pattern. */
    frgb_t fg;         /* Foreground color of disk pattern. */
  } multifok_test_disk_t;
  /* Specifies a horizontal disk with center {ctr}, radius {rad}, and color {bg[]}. */ 

typedef struct multifok_test_scene_t
  { interval_t box[3];
    int32_t ND;
    multifok_test_disk_t *dsk;
    frgb_t bg;
    frgb_t fg;
  } multifok_test_scene_t;
  /* A test scene to generate test images from. It consists of {ND} disks
    entirely contained in the given {box}, whose lowest {Z} must be non-negative;
    plus a floor plane at {Z=0}. */ 

typedef double multifok_test_pattern_t(double x, double y, double z);
  /* Type of a function that generates a grayscale pattern in {[0_1]} as a function of {(X,Y,Z)}. */

multifok_test_scene_t *multifok_test_scene_throw
  ( interval_t box[],
    int32_t ND,
    double rMin, 
    double rMax,
    double minSep, 
    bool_t verbose
  );
  /* Generates a test scene with up to {ND} disks picked at random
    inside the given {box}. The disks will be generated with
    {multifok_test_disk_throw(box,rMin,rMax)}. The backplane's {bg} and
    {fg} colors are chosen at random.
    
    If {minSep} is negative, the disks may overlap in {X} and {Y}. If
    {minSep} is zero or positive, the {X} and {Y} projections of the
    disks will be disjoint and separated by at least {minSep} pixels. In
    this case the number of disks may be less than {ND}.
    
    In any case the {box} must be strictly wider than {2*rMax} in {X} and {Y}. */ 

multifok_test_disk_t multifok_test_disk_throw(interval_t box[], double rMin, double rMax);
  /* Generates a random disk entirely contained in the given {box},
    with a random center, a  random radius in {[rMin _ rMax]}, and random
    but contrasting {bg} and {fg} colors. The {box} must be strictly wider 
    than {2*rMax} in {X} and {Y}. */

void multifok_test_images_make
  ( int32_t NX, 
    int32_t NY, 
    multifok_test_scene_t *scene,
    multifok_test_pattern_t *pattern,
    double rPix,
    double zFoc, 
    double zDep, 
    int32_t NR,
    float_image_t **cimg_P,
    float_image_t **rimg_P
  );
  /* Returns images {cimg} and {rimg} of the given {scene},
    with simulated depth-of-focus blur.
  
    Both images will have {NX} columns of {NY} pixels. The {cimg} image,
    returned in {*cimg_P}, will have 3 channels, and shows the color of
    the scene as seen through a camera or muscope focused at {Z=zFoc}
    with focal depth {zDep}. The {rimg} image, returned in {*rimg_P},
    will have 1 channel, and shows the estimate blur radius {rBlur} 
    of the image at each pixel.

    Each pixel of the two images is computed by tracing {NR} rays
    through a small disk centered above that pixel with blurry radius
    {rPix} and {Z = zFoc}, and averaging the colors and blur radii of
    the point on the topmost disk hit by each ray.
    
    !!! The fuzzy pixel should have a 2D Hahn disribution extending 1 pix on either side of center. !!!    
    
    The rays will have random directions that deviate from the vertical
    by a distance of about 1 pixel (fuzzy) over a vertical travel of {zDep}.
    
    !!! The tilt of the ray should be limited to {2/zDep} pixels, circular, with apodizing. !!!    
    
    The color of the scene at any point is determined by computing
    {r = pattern(x,y,z)} where {x,y,z} are the coordinates of the hit point,
    and then using {r} as the proportion to mix the {bg} ({r=0}) and
    {fg} ({r=1}) colors of the hit object (disk or backplane).
    */

#endif
