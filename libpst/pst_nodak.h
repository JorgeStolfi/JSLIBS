#ifndef pst_nodak_H
#define pst_nodak_H

/* Extracting sub-images of the N-Spot grayscale chart  and elements thereof. */
/* Last edited on 2024-12-22 22:31:42 by stolfi */

#include <bool.h>
#include <r2.h>
#include <r3.h>
#include <r3x3.h>
#include <float_image.h>
#include <argparser.h>
    

/* GEOMETRIC CORRECTION */

r3x3_t pst_nodak_get_matrix(r2_vec_t* geo_ctr, double geo_radius, uint32_vec_t* img_num, r2_vec_t* img_ctr);
  /* Computes the {3×3} projective map matrix that maps
    coordinates on the N-Spot chart to image coordinates. 
    
    Requires the numbers {img_num.e[..]} and center coordinates
    {img_ctr.e[...]}  of 4 or more spots in the image, not all
    collinear, and a table {geo_ctr.e[k]} of all spot centers in the
    chart, indexed by spot number {k}.

    Chart coordinates are in arbitrary real-world units
    (e.g. millimeters), with the X and Y axes oriented as 
    in mathematics.

    Image coordinates are in pixels, relative to the bottom left
    corner, with the X axis pointing to the right and the Y axis
    pointing up. The pixel with column index {c} and row index {r} in
    the image is assumed to have corners with coordinates {(c,r)} and
    {(c+1,r+1)}. So an image with {NX} columns and {NY} rows covers a
    rectangle with corners {(0,0)} and {(NX,NY)}.
    
    The chart's projection on the image may be in arbitrary
    perspective projection, but should be flat and free from radial
    (barrel/pincushion) distortions --- i.e. straight lines on the
    chart should be straight in {img} too.

    The {geo_radius} parameter is the radius of the chart, used only for debugging

 */

float_image_t *pst_nodak_extract_chart
  ( float_image_t *img,   /* Photo of a scene that includes a N-Spot chart. */
    double rad,           /* Chart radius in chart coordinates. */
    r3x3_t *C2I,          /* Chart-to-image projective map matrix. */
    uint32_t OSZ               /* Width and height of output image (pixels) */
  );
  /* Extracts from image {img} the sub-image of the N-Spot chart,
    correcting for perspective distortion. Requires the {3×3} homogeneous
    projective matrix {C2I} that maps the chart coordinates to image
    coordinates. The extracted image will cover a circle with center
    {(0,0)} and the given {rad} in chart coordinates. The output image will
    have {OSZ×OSZ} pixels. */

float_image_t *pst_nodak_extract_gray_scale
  ( float_image_t *img,    /* Photo of a scene that includes a N-Spot chart. */
    r2_vec_t* geo_ctr,     /* Center of each spot in chart coordinates. */
    double_vec_t* geo_rad, /* Radius of each spot in chart coordinates. */
    r3x3_t *C2I,           /* Chart-to-image projective map matrix. */
    double mrg,            /* Safety margin width in pixels. */
    uint32_t NX,                /* Width of each patch in the output image (pixels) */
    uint32_t NY                 /* Height of each patch in the output image (pixels) */
  );
  /* Extracts from {img} the gray-scale proper contained in the N-spot
    chart. 
    
    Requires the centers and radii of the spots in the chart
    coordinate system, and the homogeneous projective matrix {C2I} that
    maps the chart coordinates to image coordinates.

    The extracted image will be the left-to-right concatenation of
    {NS} rectangular patches, where {NS = geo_ctr->ne} and each patch
    has {NX × NY} pixels.  The sample values of patch number {i} are
    an average of the sample values in spot number {i}. */

    
#endif

    
