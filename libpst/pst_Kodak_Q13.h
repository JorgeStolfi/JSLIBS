#ifndef pst_Kodak_Q13_H
#define pst_Kodak_Q13_H

/* Extracting sub-images of the Kodak Q-13 chart and elements thereof. */
/* Last edited on 2008-12-12 16:21:17 by stolfi */

#include <bool.h>
#include <r2.h>
#include <r3.h>
#include <r3x3.h>
#include <float_image.h>
#include <argparser.h>
    
/* KODAK Q-13 DATA */
  
#define pst_Kodak_Q13_total_width   (203.0)
#define pst_Kodak_Q13_total_height  (60.0)
  /* External X and Y dimensions of the Q-13 (mm) */

#define pst_Kodak_Q13_num_steps   (20)
  /* Number of patches in gray scale, including `white' and `black'. */

#define pst_Kodak_Q13_end_patch_width (11.0)
  /* Width of first and last gray-scale patches (mm). */

#define pst_Kodak_Q13_mid_patch_width  (10.0)
  /* X extent of each gray-scale patch, except first and last (mm). */

#define pst_Kodak_Q13_patch_height  (22.0)
  /* Y extent of each gray-scale patch (mm). */

#define pst_Kodak_Q13_scale_height  (pst_Kodak_Q13_patch_height)
  /* Y extent of gray-scale (mm). */

#define pst_Kodak_Q13_patch_trim_x (1.5)
  /* Safety border to trim from left and right edges of each patch (mm). */

#define pst_Kodak_Q13_patch_trim_y (2.0)
  /* Safety border to trim from top and bottom edges of each patch (mm). */

#define pst_Kodak_Q13_ref_strip_gap (2.0)
  /* Safety distance between patch edges and strips (mm). */

#define pst_Kodak_Q13_ref_strip_height (3.0)
  /* Y extent of reference strips (mm). */

#define pst_Kodak_Q13_body_reflectance (0.195)
  /* Nominal reflectance of the gray background of the Q-13, above the gray-scale. */

/* GEOMETRIC CORRECTION */

void pst_Kodak_Q13_get_matrix(r2_t *tl, r2_t *tr, r2_t *bl, r2_t *br, r3x3_t *P);
  /* Computes the {3 × 3} projective map matrix {P} that maps
    coordinates on the Kodak Q-13 chart to image coordinates. 
    
    Requires the image coordinates of the chart's corners: top left
    {tl}, top right {tr}, bottom left {bl}, and bottom right {br}. The
    {tl}--{tr} edge should be on the lettering side of the chart, and
    {bl}--{br} should be adjacent to the gray patch scale. In both
    cases, the `left' end of the chart is that of the `white' patch.

    Chart coordinates are in millimeters; the X axis runs along the
    gray scale, from from the bottom left corner (adjacent to the
    `white' patch) to the bottom right corner (adjacent to the `black'
    patch). The Y axis rouns across the chart, from the bottom left
    corner to the top left corner (next to the word "Kodak").

    Image coordinates are in pixels. The pixel with column index {c}
    and row index {r} in the image is assumed to have corners with
    coordinates {(c,r)} and {(c+1,r+1)}. So an image with {NX} columns
    and {NY} rows covers a rectangle with corners {(0,0)} and
    {(NX,NY)}.
    
    The chart's image may be in arbitrary perspective projection, but
    should be flat and free from radial (barrel/pincushion)
    distortions --- i.e. straight lines on the chart should be
    straight in {img} too. */

float_image_t *pst_Kodak_Q13_extract_chart
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize      /* Output pixel size (mm). */
  );
  /* Extracts from {img} the image of th Kodak Q-13 chart. Requires
    the homogeneous projective matrix {P} that maps the chart
    coordinates to image coordinates. The extracted image will have
    pixels measuring {pixelSize} mm along each axis. */

float_image_t *pst_Kodak_Q13_extract_gray_scale
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize      /* Output pixel size (mm). */
  );
  /* Extracts from {img} the gray-scale proper contained in the Kodak
    Q-13 chart. Requires the homogeneous projective matrix {P} that
    maps the chart coordinates to image coordinates. The extracted
    image will have pixels measuring {pixelSize} mm along each axis. */

float_image_t *pst_Kodak_Q13_extract_body_strip
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize      /* Output pixel size (mm). */
  );
  /* Extracts from {img} the `body' reference strip, contained in the
    Kodak Q-13 chart just above the gray-scale patches. Requires the
    homogeneous projective matrix {P} that maps the chart coordinates
    to image coordinates. The extracted image will have pixels
    measuring {pixelSize} mm along each axis. */

float_image_t *pst_Kodak_Q13_extract_frame_strip
  ( float_image_t *img,   /* Photo of a scene that includes a Kodak Q-13 chart. */
    r3x3_t *P,            /* Chart-to-image projective map matrix. */
    double pixelSize      /* Output pixel size (mm). */
  );
  /* Extracts from {img} the `frame' reference strip, just below the
    gray-scale patches of the Kodak Q-13 chart. Requires the
    homogeneous projective matrix {P} that maps the chart coordinates
    to image coordinates. The extracted image will have pixels
    measuring {pixelSize} mm along each axis. */

void pst_Kodak_Q13_copy_rectangle
  ( float_image_t *img,
    r3x3_t *P,
    double xlo,
    double xhi,
    double ylo,
    double yhi,
    int x0,
    int y0,
    int SNX,
    int SNY,
    float_image_t *DST
  );
  /* Extracts from {img} the rectangle {R = [xlo _ xhi] × [ylo _ yhi]}, mapped by the 
    homogeneous projective matrix {P}, and stores it into a rectangular 
    subarray of {DST}.
    
    The destination subarray consists of {SNX} columns and {SNY} rows
    of image {DST}, starting with column {x0} and row {y0}. The rectangle {R} is
    implicitly divided into a grid with {SNX} columns and {SNY} rows,
    and each affected pixel of {DST} will get the average of the pixel
    values of {img} within the corresponding grid cell. */

double pst_Kodak_Q13_patch_albedo(int i);
  /* Nominal albedo (reflectivity) of patch number {i}.  The patches
    are numbered from 0 (lightest) to {pst_Kodak_Q13_num_steps - 1}
    (darkest). */

/* PHOTOMETRIC CORRECTION */

typedef enum /* Kind of function basis to use for function approximation: */
  { pst_Kodak_Q13_btype_A,  /* Quadratic sigmoid C1 splines. */
    pst_Kodak_Q13_btype_B   /* Sigmoid elements derived from {erf}. */
  } pst_Kodak_Q13_btype_t;

typedef struct pst_Kodak_Q13_basis_t /* Function basis for function approximation. */
  { pst_Kodak_Q13_btype_t bt;  /* Kind of basis. */
    int N;                 /* Number of element in basis. */
  } pst_Kodak_Q13_basis_t;
  /* A basis for least-squares fitting. The basis has {N} elements.
    Element {N-1} is the unit function. */

void pst_Kodak_Q13_compute_light_map
  ( int c,                     /* Channel of {img} to consider. */
    double noise,              /* Noise level to assume in {img} sample values. */
    float_image_t *grayScale,  /* The extracted and rectified gray-scale patches. */
    float_image_t *bodyStrip,  /* The extracted and rectified body strip, or NULL. */
    float_image_t *frameStrip, /* The extracted and rectified frame strip, or NULL. */
    pst_Kodak_Q13_basis_t B,  /* Function basis to use for fitting. */
    bool_t monotonic,     /* TRUE forces the map to be monotonic. */
    double z[],           /* (OUT) Coefficient vector. */
    double *logVlo,       /* (OUT) Low end of grid. */
    double *logVhi        /* (OUT) Hight end of grid. */ 
  );
  /* Computes the coefficients {z[0..B.N-1]} and the parameters
    {logVlo,logVhi} of a logscale spline function that maps sample values of
    channel {c} of image {img} to actual intensities in a linear
    scale, based on the appearance of the gray-scale image {grayScale}
    extracted from a photo of the Kodak Q-13 chart.
    
    Assumes that the sample values are a reasonably smooth and
    monotonic function of the physical intensities. Assumes also that
    the gray-scale image is entirely valid and free from projected
    shadows or highlights.
    
    Every sample value {V} in the input images must be non-negative,
    and is assumed to represent the interval {[V _ V + noise]}.
    
    The scene's illumination does not have to be uniform, as long as
    it varies smoothly across the images. If {bodyStrip} is not NULL,
    it should be a strip extracted from the chart's image, in the
    medium gray area just above the gray-scale patches. If
    {frameStrip} is not NULL, it should be a strip of the background
    just below the gray-scale patches, assumed to be a matte white
    frame. These strips are used to estimate the variations in the
    intensity of the light incident. The meaning of the parameters
    {z,logVlo,logVhi} is defined internally by {pst_Kodak_Q13_eval_map}
    below.
  */
  
double pst_Kodak_Q13_eval_map
  ( double V, 
    double noise, 
    double logVlo, 
    double logVhi, 
    pst_Kodak_Q13_basis_t B, 
    double z[]
  );
  /* Returns the linear intensity that corresponds to sample value {V}
    according to the light map defined by the parameters {noise,logVlo,logVhi},
    the basis {B}, and the coefficient vector {z[0..B.N-1]}. */
    
void pst_Kodak_Q13_apply_map
  ( float_image_t *img,
    int c, 
    double noise, 
    double logVlo, 
    double logVhi,
    pst_Kodak_Q13_basis_t B,
    double z[]
  );
  /* Maps each sample of channel {c} of image {img} through the value
    map defined by the the parameters {noise,logVlo,logVhi}, the basis {B},
    and the coefficient vector {z[0..B.N-1]}. */
    
double pst_Kodak_Q13_eval_raw_map(double v, pst_Kodak_Q13_basis_t B, double z[]);
  /* Evaluates the logscale correction function {h(v) = log(H(V))}
    determined by {B} and {z}, where {v} is a logscaled and normalized
    version of the sample value {V}. Namely, returns the linear
    combination {pst_Kodak_Q13_eval_basis(v, p, B)} with coefficient
    {z[p]}, for {p} in {0..B.N-1}. */
    
double pst_Kodak_Q13_eval_basis(double v, int p, pst_Kodak_Q13_basis_t B);
  /* Evaluates element {p} of the least-squares fitting basis {B}, at
    the logscaled and normalized argument {v}. Fails if {p} is not in
    {0..B.N-1}. */

#endif

    
