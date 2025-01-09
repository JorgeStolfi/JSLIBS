#ifndef pst_proc_map_H
#define pst_proc_map_H

/* pst_proc_map.h -- procedures for creating procedurally-defined images. */
/* Last edited on 2025-01-08 17:35:13 by stolfi */

#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <float_image.h>

typedef struct pst_proc_map_sampling_t
  { int32_t N;     /* Number of sampling points per axis. */
    double *d; /* Sampling positions relative to center point (indexed {0..NS-1}). */
    double *w; /* Sampling weights (indexed {0..NS-1}). */
  } pst_proc_map_sampling_t;
  /* A one-dimensional sampling kernel. The displacements {d[0..N-1]} are in pixels.
    The weights are in {[0_1]}, not particularly normalized. */

typedef void pst_proc_map_zfunc_t (r2_t p, double *z, r2_t *dz);
  /* A function that returns a height field {Z} and its gradient at a
    given point {p} of the plane.
    
    A function of this type should should return in {*z} the height
    {Z(p)}, which should be a number in the range {[-1 _ +1]}. If {dz}
    is not NULL, it should also return in {*dz} the gradient of {Z}
    with respect to {p}. */

void pst_proc_map_make_images
  ( pst_proc_map_zfunc_t *func,
    int32_t NX,
    int32_t NY,
    pst_proc_map_sampling_t smpZ,
    pst_proc_map_sampling_t smpG,
    bool_t numGrad,
    double maxGDiff,
    double sigmaG,
    double sigmaW,
    float_image_t *IZ, 
    float_image_t *IG,
    float_image_t *IW,
    float_image_t *IN
  );
  /* 
    Create a height map {Z(X,Y)}, the associated gradient map {G(X,Y)
    = dZ/dX(X,Y),dZ/dY(X,Y)}, a reliability weight map {W(X,Y)}, and a
    normal map {N(X,Y)}, from the given function {func}.
    
    All images refer to the same grid {M} with unit square cells,
    spanning a rectangle {[0 _ NX]×[0 _ NY]}.
    
    The height {Z(X,Y)} is nominally computed at the grid vertex
    {(X,Y)} and stored in cell {[X,Y]} of the grid of image {IZ},
    which must be non-NULL and must have {NX+1} columns and {NY+1} rows.
    The sampling is controlled by {smpZ}.
    
    If {IG} is not null, the procedure stores into {*IG} the gradient
    {G(X,Y)} of the height function {Z(X,Y)}. The {numGrad} parameter
    specifies how the gradient is computed. If {numGrad} is false, the
    procedure samples the analytic gradient as computed by {func} in
    the neighborhood of the cell centers, using the sub-sampling
    points and weights defined by {smpG}. If {numGrad} is true, each
    pixel of {IG} is set to a numeric gradient estimate computed by
    averaging the differences between the {IZ} heights at the four
    corners of the corresponding grid cell. In either case, the float image {IG} must
    have two channels (for the X and Y derivatives), {NX} columns, and
    {NY} rows.
    
    If {IW} is not null, the procedure fills {*IW} with the gradient
    reliability weights. If {numGrad} is true, the weights {IW[X,Y]}
    are 1 everywhere. If {numGrad} is false, each weight {IW[X,Y]} is
    computed by comparing the analytic gradient {IG[X,Y]} with the
    numerical gradient computed from the {IZ} map. In either case, the
    image {IW} must have one channel, {NX} columns and {NY} rows.
    
    The parameters {sigmaG} and {sigmaW} control the amount of noise
    to be introduced in the gradient and weight maps, respectively.
    
    If {IN} is not null, the procedure stores into {*IN} the unit
    normal vectors that correspond to the gradients stored in {IG}
    (after adding the noise, if any). The {Z} coordinate of the normal
    vector will be always non-negative. The image {IN} must have three
    channels (for the X, Y, and Z coordinates of the normal vector),
    {NX} columns, and {NY} rows.
    
    The procedure assumes that the interesting part of the function
    occurs when the argument {p=(x,y)} spans the signed unit square
    {U^2=[-1_+1]×[-1_+1]} of the plane. Therefore, when calling the
    {func} procedure, the image domain coordinates {P=(X,Y)} are
    uniformly scaled down and shifted, in such a way that the scaled
    grid {M} snugly surrounds the square {U^2} and is concentric with
    it. The function values returned by {func} are then scaled up by
    the same factors, so that the computed gradient remains consistent
    with the image coordinate system. */
  
pst_proc_map_sampling_t pst_proc_map_make_sampling_tables(int32_t L, int32_t N);
  /* Returns a sampling table {smp} for the given sampling parameters.
    Allocates {smp->d} and {smp->w} and fills them with proper values
    as specified by {L} and {N}. The parameter {L} defines the width
    of the sampling kernel (in pixels) and its degree. Current values
    are {0..3}:
      
      * {L=0} - point sampling, assumes {N=1}.
      * {L=1} - uniform over one pixel.
      * {L=2} - linear tent over 2 pixels.
      * {L=3} - quadratic bell over 3 pixels. */

void pst_proc_map_free_sampling_tables(pst_proc_map_sampling_t *smp);
  /* Reclaims the space used by the internals tables of {*smp} (but NOT
     the {*smp} record itself). */

double pst_proc_map_compute_height_value
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    pst_proc_map_sampling_t *smp,
    double pxPerUnit
  );
  /* Computes the height field defined by {func} at the grid corner point {p}.
    
    Evaluates {f} at a grid of {N}×{N} subsampling points near 
    {p} with displacements {smp->d[0..N-1]/pxPerUnit} and averages them with the
    weights {smp->w[0..N-1]}. The count {N} is {smp->N}, which must be odd. */
  
r2_t pst_proc_map_compute_analytic_pixel_gradient
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    pst_proc_map_sampling_t *smp,
    double pxPerUnit
  );
  /* Computes the mean gradient of the height field defined by {func}
    at the point {p}. 
    
    Evaluates the gradient of {f} at a grid of {N}×{N} subsampling points near 
    {p} with displacements {r*smp->d[0..N-1]/pxPerUnit} and averages them with the
    weights {smp->w[0..N-1]}. The count {N} is {smp->N}, which must be odd. */

r2_t pst_proc_map_compute_numeric_pixel_gradient
  ( double f00,
    double f10,
    double f01,
    double f11
  );
  /* Estimates the gradient of the height in a pixel by numerical
    differences from the four height values at the corners of the
    pixel: {f00} at bottom left, {f10} at bottom right, {f01} at top
    left, and {f11} at top right. */

double pst_proc_map_compute_pixel_weight
  ( pst_proc_map_zfunc_t *func,
    r2_t *dza,
    r2_t *dzn,
    double maxGDiff
  );
  /* Computes the reliability weight of a pixel given the average 
    gradient {*dza} in the pixel and the numerical gradient {*dzn}.
    The weight depends on the Euclidean norm {e} of the difference between the 
    two gradients, and is zero if {e >= maxGDiff}. */

void pst_proc_map_perturb_gradient(r2_t *dz, double sigma);
  /* Adds to each component of {*dz} a random value
    with Gaussian distribution, mean 0 and deviation {sigma}. */

void pst_proc_map_perturb_weight(double *w, double sigma);
  /* Multiplies {*w} by {max(0,1-S)} where {S} is a 
    random value with log-Gaussian distribution with mean of log 0 
    and deviation of log {sigma}. */

void pst_proc_map_function_00(r2_t p, double *z, r2_t *dz); /* 00 "zeflat" Constant function (zero gradient). */
void pst_proc_map_function_01(r2_t p, double *z, r2_t *dz); /* 01 "ramp10" Linear ramp in the X direction. */
void pst_proc_map_function_02(r2_t p, double *z, r2_t *dz); /* 02 "ramp01" Linear ramp in the Y direction. */
void pst_proc_map_function_03(r2_t p, double *z, r2_t *dz); /* 03 "ramp11" Linear ramp in the XY diagonal direction. */
void pst_proc_map_function_04(r2_t p, double *z, r2_t *dz); /* 04 "parabo" Slanted elliptical parabolic hump. */
void pst_proc_map_function_05(r2_t p, double *z, r2_t *dz); /* 05 "spdom1" Buried sphere. */
void pst_proc_map_function_06(r2_t p, double *z, r2_t *dz); /* 06 "pyram5" Pentagonal pyramid. */
void pst_proc_map_function_07(r2_t p, double *z, r2_t *dz); /* 07 "rtcone" Conical mound. */
void pst_proc_map_function_08(r2_t p, double *z, r2_t *dz); /* 08 "ripple" Circular waves. */
void pst_proc_map_function_09(r2_t p, double *z, r2_t *dz); /* 09 "sbabel" Tower of Babel with smooth soulders. */
void pst_proc_map_function_10(r2_t p, double *z, r2_t *dz); /* 10 "hcliff" Diverging parabolic ramps with cliff. */
void pst_proc_map_function_11(r2_t p, double *z, r2_t *dz); /* 11 "sinw01" Sinusoidal wave, low freq. */
void pst_proc_map_function_12(r2_t p, double *z, r2_t *dz); /* 12 "cosw01" Co-sinusoidal wave, low freq. */
void pst_proc_map_function_13(r2_t p, double *z, r2_t *dz); /* 13 "sinw02" Sinusoidal wave, med freq. */
void pst_proc_map_function_14(r2_t p, double *z, r2_t *dz); /* 14 "cosw02" Co-sinusoidal wave, med freq. */
void pst_proc_map_function_15(r2_t p, double *z, r2_t *dz); /* 15 "sinw03" Sinusoidal wave, high freq. */
void pst_proc_map_function_16(r2_t p, double *z, r2_t *dz); /* 16 "cosw03" Co-sinusoidal wave, high freq. */
void pst_proc_map_function_17(r2_t p, double *z, r2_t *dz); /* 17 "cbabel" Tower of Babel with cliff. */
void pst_proc_map_function_18(r2_t p, double *z, r2_t *dz); /* 18 "cbramp" Cubic ramp with cliff on three sides. */
void pst_proc_map_function_19(r2_t p, double *z, r2_t *dz); /* 19 "holes1" Buried sphere with many holes in slope map. */
void pst_proc_map_function_20(r2_t p, double *z, r2_t *dz); /* 20 "cpiece" Three nested stages connected by ramps. */
void pst_proc_map_function_21(r2_t p, double *z, r2_t *dz); /* 21 "cplat3" Three flat round platforms with smooth edges. */
void pst_proc_map_function_22(r2_t p, double *z, r2_t *dz); /* 22 "plat3a" Triangular wall with smooth edges. */
void pst_proc_map_function_23(r2_t p, double *z, r2_t *dz); /* 23 "plat5a" Pentagonal platform with smooth edges. */
void pst_proc_map_function_24(r2_t p, double *z, r2_t *dz); /* 24 "cplatv" Circular platform with sharp edges. */
void pst_proc_map_function_25(r2_t p, double *z, r2_t *dz); /* 25 "cplats" Circular platform with smooth edges. */
void pst_proc_map_function_26(r2_t p, double *z, r2_t *dz); /* 26 "qtbell" Quartic bell. */
void pst_proc_map_function_27(r2_t p, double *z, r2_t *dz); /* 27 "fracto" Fractalish round mountain. */
  /* The built-in height map functions. They assume that the signed unit 
    square {[-1_+1]×[-1_+1]} fits snugly inside the image's domain.
    
    !!! Each function should define a weight too. !!!
    !!! Either that, or let {z=NAN} to get weight 0. !!!
    !!! Each function should be told the resolution of the map? !!! */

#define pst_proc_map_MIN_ZFUNC  0
#define pst_proc_map_MAX_ZFUNC 27
  /* Min and max numbers of {n} in {pst_proc_map_function_{n}}. */

pst_proc_map_zfunc_t *pst_proc_map_function_generic(int32_t n);
  /* Returns a pointer to {pst_proc_map_function_{n}}. */

char* pst_proc_map_function_generic_name(int32_t n);
  /* Returns the name of the function {n}, as per the table above. */
  
/* Parametrized procedure maps: */

void pst_proc_map_function_affine(r2_t *p, r2_t *a, double fa, r2_t *b, double fb, r2_t *c, double fc, double *z, r2_t *dz);
  /* Computes the affine function {z = A*x + B*y + C} that has values {fa,fb,fc}
    at points {a,b,c} (which must noy be collinear).  Also computes its gradient {*dz}. */

void pst_proc_map_function_wave(r2_t *p, r2_t *f, double phase, double *z, r2_t *dz);
  /* Computes the generic wave function {*z = 0.5*sin(2*PI*dot(p,f) + phase)}.
    Also computes its gradient {*dz}. */

void pst_proc_map_function_babel(r2_t *p, double RI, double RO, double N, double EF, double *z, r2_t *dz);
  /* Computes the generic Tower of Babel function {*z} with inner radius {RI},
    outer radius {RO}, a ramp with {N} full turns. Also computes its gradient
    {*dz}.  
    
    If {EF = 0}, the ramp is bounded by a vertical cliff on both
    sides. If {EF > 0} a filling is added adjacent to the cliff, to
    smooth it out. The filling width is {EF} times the bare ramp
    width, so the max value of {EF} is 1. */

void pst_proc_map_function_fractal_mound(r2_t *p, double *z, r2_t *dz);
  /* A round fractalish mound with radius 1, center height 1, background height 0. */
  
void pst_proc_map_function_cubic_ramp(double x, double *z, double *dz);
  /* Computes a smooth cubic ramp {*z} and its derivative {*dz}.
    The ramp is a function of {x} only; it has value 0 and slope 0
    at {x <= -1}, and value {1} and slope 0 at {x = +1}. */

void pst_proc_map_function_round_platform(r2_t *p, double R, double HS, double *z, r2_t *dz);
  /* Computes a flat circular platform function {*z} with height 1 and radius {R}. 
    
    If {HS} is zero the edges are sharp steps; if {HS>0} they are
    soft cubic step extending by distance {HS} on either side of the
    ideal perimeter. Also computes the gradient {*dz}. */

void pst_proc_map_function_polygonal_platform(r2_t *p, int32_t N, double R, double tilt, double S, double *z, r2_t *dz);
  /* Computes a flat polygonal platform function {*z} with height 1. 
    
    The basic outline of the platform will be a regular polygon with
    circumradius {R}, {N} sides, and one vertex rotated {tilt}
    degrees counterclockwise from the {X} axis.
    
    If {S} is zero the edges are sharp steps; if {S>0} they are
    soft cubic step extending by distance {S} outwards of the
    ideal outline. Also computes the gradient {*dz}. */

#endif
