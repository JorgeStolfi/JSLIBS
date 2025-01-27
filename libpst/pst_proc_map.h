#ifndef pst_proc_map_H
#define pst_proc_map_H

/* pst_proc_map.h -- procedures for creating procedurally-defined images. */
/* Last edited on 2025-01-26 18:59:53 by stolfi */

#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <float_image.h>

typedef void pst_proc_map_zfunc_t (r2_t p, double *z, r2_t *dz);
  /* A function that returns a height field {Z} and its gradient at a
    given point {p} of the plane.
    
    A function of this type should should return in {*z} the height
    {Z(p)}, which should be a number in the range {[-1 _ +1]}. If {dz}
    is not NULL, it should also return in {*dz} the gradient of {Z}
    with respect to {p}.  If the height is undefined or uncomputable, 
    it should return {NAN} in both. */
     
typedef struct pst_proc_map_zfunc_props_t 
  { pst_proc_map_zfunc_t *func;  /* The function. */
    int32_t num;                 /* ID number. */
    char *name;                  /* Function name (6 alphanum chars). */
    double maxGrad;              /* Maximum valid gradient value. */
    double maxGDiff;             /* Maximum valid numeric/analytic gradient mismatch. */
  } pst_proc_map_zfunc_props_t;
  /* Main attributes of a procedural map.

    The {maxGrad} and {maxGDiff} attributes are the proper parameters for the
    paint procedures below. */
     
/* AVERAGED HEIGHTS AND GRADIENTS */

/* The procedures in this section compute either the mean height {Z}
  or the mean gradient {(dZ/dX,dZ/dY)} in the neighborhood of a
  given point {p} of {func}'s domain.
  
  The result is obtained by evaluanting {func} at a grid of {NS×NS}
  sampoints (sampling points) around the given point {p} and averaging
  the desired quantity. The count {NS} must be at
  least 2. The sampoints will span about {1/xyScale} on each side 
  of {p}.   The weight of the sampoint with indices {kx,ky}
  will be {ws[kx]*ws[ky]}, for {kx,ky} in {0..NS-1}.  Preferably,
  {NS} should be odd, and the weights {ws[0..NS-1]} should be a 
  partition of unity with stride {NS/2+1}.  See {wt_table_hann_fill}
  for a suitable instance. */ 

void pst_proc_map_compute_height
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    uint32_t NS,
    double ws[],
    double xyScale,
    bool_t debug,
    double *z_P,
    double *wz_P
  );
  /* Computes the height field defined by {func} at the point {p}
    of func's domain.
    
    The result {z}, returned in {*z_P}, is the weighted average of the
    height function {Z(X,Y)} as defined by {func} over the grid of
    {NS×NS} sampoints with weights {ws[0..NS-1]}.  Sampoints
    where {func} returns {NAN} or {±INF} are excluded from the average.
    
    The height is scaled by {xyScale} so that the gradient is 
    independent of {xyScale}.
    
    The procedure also computes a reliability weight {wz}, that is
    returned in {*wz_P}, consisting of the weighted fraction of
    sampoints where {func} returned a finite value. In particular, the
    weight {wz} will be 1.0 if {func} did not return {NAN} or {±INF} at
    any sampoint. */
  
void pst_proc_map_compute_numeric_gradient
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    uint32_t NS,
    double ws[],
    double xyScale,
    double maxGrad,
    bool_t debug,
    r2_t *dnz_P,
    double *wdnz_P
  );
  /* Computes the mean numeric gradient {dnz=(dnzx,dnzy)} of the height field defined by
    {func} at the point {p} of {func's} domain. The gradient is returned
    in {*dnz_P}.
    
    The result is obtained by computing the height {Z(X,Y)} with {func}
    at a grid of {NS×NS} sampoints {(X,Y)} around {p}, with weights
    {ws[0..NS-1]}; computing the numeric derivatives {dnxk=dZ/dX} and
    {dnyk=dZ/dY} by divided differences between pairs of adjacent sampoints; and
    averaging these numeric derivatives.
    
    The procedure computes a reliability weight {wdnz}, that is returned
    in {*wdnz_P}. If any sampoint heights are {NAN} or {±INF}, or any
    computed derivative exceeded {maxGrad} in absolute value, or the
    final gradient is greater than {maxGrad} in modulus, the gradient is
    set to {NAN,NAN} and the weight is zero; otherwise it is 1.0. */
  
void pst_proc_map_compute_analytic_gradient
  ( pst_proc_map_zfunc_t *func,
    r2_t p,
    uint32_t NS,
    double ws[],
    double xyScale,
    double maxGrad,
    bool_t debug,
    r2_t *daz_P,
    double *wdaz_P
  );
  /* Computes the mean analytic gradient {daz=(dazx,dazy)} of the height
    field defined by {func} at the point {p} of {func's} domain. The
    gradient is returned in {*daz_P}.
    
    The gradient is obtained by computing the analyitc derivatives
    {dazxk=dZ/dX} and {dazyk=dZ/dY} at each sampoint, using {func}, and
    taking the average of those valus with the sampoint weights.
    
    The procedure computes a reliability weight {wdnz}, that is returned
    in {*wdnz_P}. If any sampoint heights are {NAN} or {±INF}, or any
    gradient returned by {func} exceeds {maxGrad} in modulus, the gradient is
    set to {NAN,NAN} and the weight is zero; otherwise it is 1.0. */
    
/* CREATING IMAGES FROM PROCEDURAL MAPS */

/* The procedures in this section assumes that the interesting part of the
  function occurs when the argument {p=(x,y)} spans the signed unit square
  {U^2=[-1_+1]×[-1_+1]} of the plane. Therefore, when calling the
  {func} procedure, the image domain coordinates {P=(X,Y)} are
  uniformly scaled and shifted, in such a way that the
  rectangle {[0 _ NX] × [0 _ NY]} is mapped to a rectangle
  that snugly surrounds the square {U^2} and is concentric with
  it. */
    
float_image_t* pst_proc_map_make_height_map
  ( pst_proc_map_zfunc_t *func,
    int32_t NX,
    int32_t NY,
    uint32_t NS,
    double ws[]
  );
  /* Creates a height map {IZ} from the height function {Z(X,Y)} defined
    by the procedure {func}.
    
    The {IZ} image will have 2 channels, {NX+1} columns, and {NY+1} rows
    (for consistency with {pst_proc_map_make_slope_map} below). 
    
    The height value of {IZ[0,x,y]} and the reliability weight
    {IZ[1,x,y]}, for {x} in {0..NX} and {y} in {0..NY}, will be obtained
    by {pst_proc_map_compute_height}, with {p} being the grid VERTEX
    {(x,y)} mapped to {func}'s domain. */

float_image_t* pst_proc_map_make_slope_map
  ( pst_proc_map_zfunc_t *func,
    int32_t NX,
    int32_t NY,
    uint32_t NS,
    double ws[],
    bool_t numGrad, /* Use numeric gradient? */
    double maxGrad,
    double maxGDiff
  );
  /* Computes a slope map {IG} from the gradient of the 
    height function {Z(X,Y)} defined by the procedure {func}.
    
    The slope map image {IG} will have three channels, {NX} columns, and
    {NY} rows.  The procedure stores into {IG[0..2,x,y]} the
    estimated derivatieves {dZ/dX} and {dZ/dY} (numeric or analytic, as
    specified by {numGrad}) and the reliability 
    weight {w}, as computed by {pst_proc_map_compute_gradient} with {p}
    being the CENTER of the pixel {[x,y]} -- that is, {(x+0.5,y+0.5)} --
    mapped to the domain of {func}.
    
    If the */

#define pst_proc_map_MIN_ZFUNC  0
#define pst_proc_map_MAX_ZFUNC 27
  /* Min and max numbers of {n} in {pst_proc_map_function_{n}}. */

pst_proc_map_zfunc_props_t pst_proc_map_function_generic(int32_t n);
  /* Returns a record with the attributes of
    the height function {pst_proc_map_function_{n}} below. */
  
/* PREDEFINED PROCEDURAL HEIGHT MAPS */

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

/* PARAMETRIZED PROCEDURAL HEIGHT MAPS: */

void pst_proc_map_function_affine(r2_t *p, r2_t *a, double fa, r2_t *b, double fb, r2_t *c, double fc, double *z, r2_t *dz);
  /* Computes the affine function {z = A*x + B*y + C} that has values {fa,fb,fc}
    at points {a,b,c} (which must noy be collinear).  Also computes its gradient {*dz}. */

void pst_proc_map_function_wave(r2_t *p, r2_t *f, double phase, double *z, r2_t *dz);
  /* Computes the generic wave function {*z = 0.5*sin(2*PI*dot(p,f) + phase)}.
    Also computes its gradient {*dz}. */

void pst_proc_map_function_babel(r2_t *p, double RI, double RO, double NS, double EF, double *z, r2_t *dz);
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
