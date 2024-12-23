#ifndef pst_normal_map_H
#define pst_normal_map_H

/* pst_normal_map.h -- procedures for working with normal maps. */
/* Last edited on 2024-12-22 12:38:29 by stolfi */

#include <float_image.h>
#include <r2.h>
#include <r3.h> 
#include <r3x3.h> 
#include <argparser.h>

#include <pst_basic.h>

/* NORMAL MAPS
  
  A /normal map/ is a three-channel float-valued image of some
  three-dimensional scene, where the value of each pixel is the
  outwards-pointing unit vector normal to the scene's surface in that
  pixel. */

r3_t pst_normal_map_get_pixel(float_image_t *NRM, int32_t x, int32_t y);
  /* Extracts channels 0..2 of the pixel in column {x}, row {y} of {NRM},
    and returns them as the X, Y, and Z coordinates of an {r3_t}. */

void pst_normal_map_set_pixel(float_image_t *NRM, int32_t x, int32_t y, r3_t *nrm);
  /* Stores the X, Y, and Z coordinates of {nrm} as channels 0..2 pixel in column {x}, 
    row {y} of {NRM}. */

/* CREATING NORMAL MAPS */

typedef r3_t pst_normal_map_proc_t (r2_t *p);
  /* A procedure that computes the normal direction {nrm} at a visible point
    {P} of a some surface, given the projection {p} of that point in some
    plane. 
    
    Both {p} and {nrm} are given in some orthogonal {U,V,W} coordinate
    system such that the projection of point {(u,v,w)} has coordinates
    {(u,v)} (i.e., such that the {W} axis is parallel to the direction
    of projection). The returned normal should have a non-negative {W}
    component.
    
    The procedure should return the null vector {(0,0,0)} if the
    normal direction is not defined at the point {P} (e.g. if {P} is
    at infinity). */ 

r2_t pst_normal_map_scene_pt_from_image_pt(r2_t *xy, r3x3_t *xym_to_uvm); 
  /* Given the coordinates {xy} of a point in the image coordinate
    system (in pixels, with origin at the bottom left corner and Y
    axis pointing up), returns its coordinates {uv} in some scene
    coordinate system.
    
    Assumes that the two coordinate systems are related by a
    projective transformation, which is computed by multiplying the row vector
    {(x,y,1)} by the matrix {xym_to_uvm}, in that order, to yield a
    row vector {(u,v,n)}; and returning the vector {(u/n,v/n)}.  */

r3_t pst_normal_map_image_dir_from_scene_dir(r3_t *uvw, r3x3_t *uvw_to_xyz); 
  /* Given the coordinates {uvw} of a vector in some scene coordinate
    system, returns its coordinates {xyz} in the image coordinate
    system.
    
    Assumes that vectors in the two coordinate systems are related by
    a similarity (rotation plus scaling), which is computed by
    multiplying the row vector {(u,v,w)} by the matrix {uvw_to_xyz},
    in that order. The result is scaled to unit norm by the procedure.
    
    As a special case, if {uvw} is the zero vector {(0,0,0)}, the
    result {xyz} is {(0,0,0)}, too. */

r3_t pst_normal_map_eval
  ( pst_normal_map_proc_t nrmf, /* Normal-computing funtion. */
    double x,                   /* X-coordinate of projected point in image system. */
    double y,                   /* Y-coordinate of projected point in image system. */
    r3x3_t *xym_to_uvm,         /* Affine map of image {X,Y} coords to scene {U,V} coords. */
    r3x3_t *uvw_to_xyz          /* Linear map of normal from {U,V,W} coords to {X,Y,Z} coords. */
  );
  /* Computes the normal direction to a scene's point, in image
    coordinates.  Applies {pst_normal_map_scene_pt_from_image_pt}
    with {xym_to_uvm}, computes the normal direction with {nrmf},
    then applies {pst_normal_map_image_dir_from_scene_dir} with
    {uvw_to_xyz} to the result.
    
    The argument of that transformation is a
    vector {(x,y,m)} where {(x/m,y/m)} is the projected point in image
    coordinates. The result is {(u,v,n)} where {(u/n,v/n)} is the
    projected point in the {U,V,W} system.
    
    The procedure {nrmf} is supposed to take a projected point {(u,v)}
    and return the normal direction as a vector {(u',v',w')}, in that
    same {U,V,W} system. The matrix {uvw_to_xyz} is then used to
    convert the vector {(u',v',w')} back to the {X,Y,Z} coord system.
    The result is normalized to unit Euclidean length and returned.
    
    In both maps, the argument is presented as a row vector, and the
    matrix is multiplied on the right side. */

void pst_normal_map_from_proc
  ( pst_normal_map_proc_t nrmf, /* Normal-computing funtion. */
    int32_t NS,                     /* Order of subsampling grid within each pixel. */
    r3x3_t *xym_to_uvm,         /* Affine map of image {xy} coords to model {uv} coords. */
    r3x3_t *uvw_to_xyz,         /* Linear map of {uvw} coords to normal {xyz} coords. */
    float_image_t *NRM          /* (OUT) Computed normal map. */
  );
  /* Computes a normal map for a scene given a procedure {nrm} that
    computes the normal direction at a given projected point.
    
    The normal vector for each pixel is obtained by averaging the
    surface normals at a grid of {NS*NS} points within that pixel,
    with {pst_normal_map_pixel_avg}, using the matrices {xym_to_uvm} and
    {uvw_to_xyz} to convert the image {X,Y,X} coords to scene's
    {U,V,W} coords and back. */

r3_t pst_normal_map_pixel_avg
  ( pst_normal_map_proc_t nrmf, /* Normal-computing funtion. */
    int32_t x, int32_t y,               /* Pixel indices (coords of lower left corner). */
    int32_t NS,                     /* Order of sub-sampling grid in pixel. */
    r3x3_t *xym_to_uvm,         /* Affine map of image {xy} coords to model {uv} coords. */
    r3x3_t *uvw_to_xyz          /* Linear map of {uvw} coords to normal {xyz} coords. */
  );
  /* Computes the outwards normal vector of a sphere's syrface,
    averaged over the patch that projects onto the pixel in column
    {x}, row {y} --- that is, over the square with diagonal corners
    {(x,y)} and {(x+1,y+1)} in the image coordinate system. 
    
    The average is taken over a grid of {NS×NS} samples inside the
    pixel. If the normal is {(0,0,0)} for any sample point, the pixel
    average is {(0,0,0)}. Otherwise, normals are averaged by
    converting them to slopes {dW/dU,dW/dV} in the {U,V,W} system,
    computing the arithmetic mean of the slopes, and converting the
    result back to a normal vector.
    
    The normal function {nrmf} is evaluated with
    {pst_normal_map_eval}, using the matrices {xym_to_uvm} and
    {uvw_to_xyz} to convert the image {X,Y coords to scene's {U,V},
    and the normal in {U,V,W} coords to scene's {X,Y,Z} coords. */

void pst_normal_map_perturb(float_image_t *NRM, double amt);
  /* Randomly perturbs the direction of the normal vectors in the
    normal map {NRM}, by the amount {amt}, as in {pst_perturb_normal}.
    Null vectors are not disturbed. */

void pst_perturb_normal(r3_t *nrm, double amt);
  /* Randomly perturbs the unit vector {nrm} by the amount {amt},
    preserving its unit length. In particular, {amt=0} means no
    perturbation, while {amt=1} means that the result is a uniformly
    distributed on the unit sphere and independent of the given {nrm}.
    For small {amt}, the perturbation is (to first order) a tangential
    vector perpendicular to {nrm}, with zero mean and root mean square
    length {amt}.  */

/* NORMAL MAPS AND SLOPE MAPS */

r2_t pst_normal_map_slope_from_normal(r3_t *nrm, double maxSlope);
  /* Converts a normal vector {nrm} to a gradient vector (that is, the
    slopes {dZ/dX} and {dZ/dY}). 
    
    The Z coordinate or the normal must be non-negative. The length of
    the normal vector needs not be 1. If necessary, the computed
    gradient is scaled so that its Euclidean norm does not exceed
    {maxSlope}. The result is undefined if {nrm} is the zero
    vector. */

r3_t pst_normal_map_normal_from_slope(r2_t *grd);
  /* Converts a gradient vector {grd} (that is, the slopes {dZ/dX} and
    {dZ/dY}) to an outwards-pointing normal vector, with unit
    Euclidean length. */

float_image_t *pst_normal_map_to_slope_map(float_image_t *NRM, double maxSlope);
  /* Given a three-channel normal map {NRM} returns a two-channel slope
    map {GRD} with the same size, such that the sample in channel {c},
    column {x} and row {y} of {GRD} is the derivative of the surface
    height at along axis {c}, computed from the corresponding normal. */

float_image_t *pst_normal_map_from_slope_map(float_image_t *GRD);
  /* Given a two-channel slope map {GRD}, returns a three-channel 
    normal map {NRM}, with the same size, such that the pixel value in
    column {x} and row {y} of {NRM} is the outward unit vector normal
    to the surface, computed from the corresponding slopes. */

#endif
