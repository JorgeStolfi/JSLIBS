#ifndef pst_gray_scale_fit_H
#define pst_gray_scale_fit_H

/* Fitting camera transfer function to data from a gray scale chart. */
/* Last edited on 2024-12-24 18:57:27 by stolfi */

#include <bool.h>
#include <r2.h>
#include <r3.h>
#include <r3x3.h>
#include <float_image.h>
#include <argparser.h>

/* PHOTOMETRIC SAMPLE MAPPING */

typedef enum /* Kind of function basis to use for function approximation: */
  { pst_gray_scale_fit_btype_A,  /* Quadratic sigmoid C1 splines. */
    pst_gray_scale_fit_btype_B   /* Sigmoid elements derived from {erf}. */
  } pst_gray_scale_fit_btype_t;

typedef struct pst_gray_scale_fit_basis_t /* Function basis for function approximation. */
  { pst_gray_scale_fit_btype_t bt;  /* Kind of basis. */
    uint32_t N;                 /* Number of element in basis. */
  } pst_gray_scale_fit_basis_t;
  /* A basis for least-squares fitting. The basis has {N} elements.
    Element {N-1} is the unit function. */
  
double pst_gray_scale_fit_eval_map
  ( double V, 
    double noise, 
    double logVlo, 
    double logVhi, 
    pst_gray_scale_fit_basis_t B, 
    double z[]
  );
  /* Returns the linear intensity that corresponds to sample value {V}
    according to the light map defined by the parameters {noise,logVlo,logVhi},
    the basis {B}, and the coefficient vector {z[0..B.N-1]}. */
    
void pst_gray_scale_fit_apply_map
  ( float_image_t *img,
    int32_t c, 
    double noise, 
    double logVlo, 
    double logVhi,
    pst_gray_scale_fit_basis_t B,
    double z[]
  );
  /* Maps each sample of channel {c} of image {img} through the value
    map defined by the the parameters {noise,logVlo,logVhi}, the basis {B},
    and the coefficient vector {z[0..B.N-1]}. */

/* FITTING A MAP TO A GRAY SCALE CHART */

void pst_gray_scale_fit_light_map
  ( int32_t c,                /* Channel of {img} to consider. */
    uint32_t NS,               /* Number of steps in gray scale. */
    double noise,              /* Noise level to assume in {img} sample values. */
    float_image_t *imgScale,   /* The extracted and rectified gray-scale patches. */
    double albScale[],         /* Nominal albedo of each patch. */
    double dX,                 /* X displ between centers of successive patches (pixels). */
    bool_t useSelf,            /* TRUE to use neighboring patches for lighting estimation. */
    float_image_t *imgStrip0,  /* Extracted image of a reference strip, or NULL. */
    double albStrip0,          /* Albedo of {imgStrip0}. */
    double dY0,                /* Y displ from center of {imgScale} to center of {imgStrip0} (pixels). */
    float_image_t *imgStrip1,  /* Extracted image of another reference strip, or NULL. */
    double albStrip1,          /* Albedo of {imgStrip1}. */
    double dY1,                /* Y displ from center of {imgScale} to center of {imgStrip1} (pixels). */
    pst_gray_scale_fit_basis_t B,  /* Function basis to use for fitting. */
    bool_t monotonic,     /* TRUE forces the map to be monotonic. */
    double z[],           /* (OUT) Coefficient vector. */
    double *logVlo,       /* (OUT) Low end of grid. */
    double *logVhi        /* (OUT) Hight end of grid. */ 
  );
  /* Computes the coefficients {z[0..B.N-1]} and the parameters
    {logVlo,logVhi} of a logscale spline function that maps sample
    values of channel {c} to actual intensities in a linear scale,
    based on the appearance of a string of gray scale patches, such as
    those provided by the Kodak Q-13 grayscale reference chart. The
    scene's illumination does not have to be uniform, as long as it
    varies smoothly across the images.
    
    The procedure assumes that the sample values are a reasonably
    smooth function of the physical intensities. If {monotonic} is
    true, also assumes that pixel values vary monotonically with
    intensity.
    
    The main input data is the image {imgScale}. It should consist of
    {NS} sub-images, all with the same number of columns,
    concatenated in the X direction. Sub-image {i} is assumed to be
    the image of a surface with uniform reflectivity {albScale[i]}.
    
    If {imgStrip0} is not NULL, it should contain another sequence of
    {NS} sub-images, also all with the same column count. Each
    sub-image is assumed to be the picture of a surface with
    reflectiveity {albStrip0}, close to patch number {i} of
    {imgScale}. The center of each patch in {imgStrip0} is assumed to
    be {dY0} pixels away in the Y direction from the center of the
    corresponding patch of {imgScale}. Ditto for {imgStrip1},
    {albStrip1}, and {dY1}. These strips are used to estimate the
    variations in the intensity of the light incident along the
    {imgScale}.
    
    All input images should be free from projected shadows or highlights.
    Every sample value {V} in those images must be non-negative,
    and is assumed to represent the interval {[V _ V + noise]}.
    
    The meaning of the parameters {z,logVlo,logVhi} is defined
    internally by {pst_gray_scale_fit_eval_map} below.
  */
    
double pst_gray_scale_fit_eval_raw_map(double v, pst_gray_scale_fit_basis_t B, double z[]);
  /* Evaluates the logscale correction function {h(v) = log(H(V))}
    determined by {B} and {z}, where {v} is a logscaled and normalized
    version of the sample value {V}. Namely, returns the linear
    combination {pst_gray_scale_fit_eval_basis(v, p, B)} with coefficient
    {z[p]}, for {p} in {0..B.N-1}. */
    
double pst_gray_scale_fit_eval_basis(double v, uint32_t p, pst_gray_scale_fit_basis_t B);
  /* Evaluates element {p} of the least-squares fitting basis {B}, at
    the logscaled and normalized argument {v}. Fails if {p} is not in
    {0..B.N-1}. */

#endif

    
