#ifndef pst_normal_map_H
#define pst_normal_map_H

/* pst_normal_map.h -- procedures for working with normal maps. */
/* Last edited on 2025-01-19 08:46:42 by stolfi */

#include <float_image.h>
#include <r2.h>
#include <r3.h> 
#include <r3x3.h> 
#include <argparser.h>

#include <pst_basic.h>

/* NORMAL MAPS
  
  A /normal map/ is a channel float-valued image of some
  three-dimensional scene, that, at each pixel, gives the outwards-pointing unit vector normal to the scene's surface in that
  pixel. 
  
  The image must have 3 or 4 channels.  Channels {0..2} of each pixel are the coordinates
  of the normal vector. Channel 3, if present, is a reliability weight for that normal,
  between 0 (normal is undefined) to 1 (normal is as good as it can be). */

r3_t pst_normal_map_get_vector(float_image_t *NRM, int32_t x, int32_t y);
  /* Extracts channels 0..2 of the pixel in column {x}, row {y} of {NRM},
    and returns them as the X, Y, and Z coordinates of an {r3_t}. */

void pst_normal_map_set_vector(float_image_t *NRM, int32_t x, int32_t y, r3_t *nrm);
  /* Stores the X, Y, and Z coordinates of {nrm} as channels 0..2 pixel in column {x}, 
    row {y} of {NRM}. */

float pst_normal_map_get_weight(float_image_t *NRM, int32_t x, int32_t y);
  /* Extracts the weight (channel 3) of the pixel in column {x}, row {y} of {NRM}.
    If {NRM} has only 3 channels, returns {1.0}. */

void pst_normal_map_set_weight(float_image_t *NRM, int32_t x, int32_t y, float w);
  /* Stores {w} channel 3 of the pixel in column {x},  row {y} of {NRM}. 
    The map must have 4 channels.  */

/* CREATING NORMAL MAPS */

r2_t pst_normal_map_scene_pt_from_image_pt(r2_t *imgp, r3x3_t *img2_to_scn2); 
  /* Given the 2D coordinates {imgp} of a point in the image coordinate
    system (in pixels, with origin at the bottom left corner and Y
    axis pointing up), returns its coordinates {scnp} in some scene
    coordinate system.
    
    Assumes that the two coordinate systems are related by a
    projective transformation, which is computed by multiplying the row vector
    {(ix,iy,1)} by the matrix {img2_to_scn2}, in that order, to yield a
    row vector {(sx,sy,sn)}; and returning the vector {(sx/sn,sy/sn)}.  */

r3_t pst_normal_map_image_dir_from_scene_dir(r3_t *scnv, r3x3_t *scn_to_img); 
  /* Given the coordinates {scnv} of a vector in some scene coordinate
    system, returns its coordinates {imgv} in the image coordinate
    system.
    
    Assumes that vectors in the two coordinate systems are related by
    a similarity (rotation plus scaling), which is computed by
    multiplying the row vector {(sx,sy,sz)} by the matrix {scn_to_img},
    in that order. The result is scaled to unit norm by the procedure.
    
    As a special case, if {scnv} is the zero vector {(0,0,0)}, the
    result {imgv} is {(0,0,0)}, too. */

r3_t pst_normal_map_eval
  ( pst_normal_func_t nrmf,   /* Normal-computing funtion. */
    double ix,                /* X-coordinate of projected point in image system. */
    double iy,                /* Y-coordinate of projected point in image system. */
    r3x3_t *img2_to_scn2,     /* Projective map of 2D image coords to 2D scene coords. */
    r3x3_t *scn_to_img        /* Linear map of vectors from 3D scene to 3D image coords. */
  );
  /* Computes the normal direction to a scene's point, in image
    coordinates.  Applies {pst_normal_map_scene_pt_from_image_pt}
    with {img2_to_scn2}, computes the normal direction and weight with {nrmf},
    then applies {pst_normal_map_image_dir_from_scene_dir} with
    {scn_to_img} to the result.
    
    The argument of the first transformation is a homogeneous coordinate
    vector {imgh=[hix,hiy,him]} where {imgp=(hix/him,hiy/him)=(ix,iy)}
    is the projected point in image coordinates. The result is
    {scnh=[hsx,hsy,hsn]} where {scnp=(hsx/hsn,hsy/hsn)=(sx,sy)} is the
    projected point in the scene's system.
    
    The procedure {nrmf} is supposed to take a projected point {(sx,sy)}
    and return the normal direction as a vector {(sx',sy',sz')}, in that
    same {SX,SY,SZ} system. The matrix {scn_to_img} is then used to
    convert the vector {(sx',sy',sz')} back to the {IX,IY,Z} coord system.
    The result is normalized to unit Euclidean length and returned.
    
    In both maps, the argument is presented as a row vector, and the
    matrix is multiplied on the right side. 
    
    !!! Fix this to use a 3D projective map, not a linear map !!! */

void pst_normal_map_from_proc
  ( pst_normal_func_t nrmf,  /* Normal-computing funtion. */
    int32_t NS,              /* Order of subsampling grid within each pixel. */
    r3x3_t *img2_to_scn2,    /* Projective map of 2D image coords to 2D scene coords. */
    r3x3_t *scn_to_img,      /* Linear map of vectors from 3D scene to 3D image coords. */
    float_image_t *NRM       /* (OUT) Computed normal map. */
  );
  /* Computes a normal map for a scene given a procedure {nrm} that
    computes the normal direction at a given projected point.
    
    The normal vector for each pixel is obtained by averaging the
    surface normals at a grid of {NS*NS} points within that pixel,
    with {pst_normal_map_pixel_avg}, using the matrices {img2_to_scn2} and
    {scn_to_img} to convert the normal vector {scnv} from scene coords to 
    image coords and back. 
    
    The map {NRM} must have 3 or 4 channels.  If it has 4 channels,
    the precedure sets channel 3 to the pixel weights returned by
    {pst_normal_map_pixel_avg}. */

void pst_normal_map_perturb(float_image_t *NRM, double amt);
  /* Randomly perturbs the direction of the normal vectors in the
    normal map {NRM}, by the amount {amt}, as in {pst_perturb_normal}.
    Null vectors are not disturbed. 
    
    The image {NRM} must have 3 or 4 channels.  If it has 4 channels,
    the procedure sets channel 3 (the weight) to 0 whenever the normal is 
    not finite or has zero length. */

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
  /* Given a normal map {NRM} returns a slope map {GRD} with 3 channels
    and the same size, such that the samples in channels 0 and 1, column
    {x} and row {y} of {GRD} are the derivatives of the surface height
    along the {X} and {Y} axes, computed from the normal vector in
    channels {0..2} of {NRM}.
    
    The normal map {NRM} must have either 3 or 4 channels.  The resulting
    slope map {GRD} will have 3 channels; channel 2 of {GRD} will be 
    the weight (channel 3) of {NRM}, or 1.0 if there is no such channel.
    However the weight will be 0 if the normal vector is not finite
    or has zero length. */

float_image_t *pst_normal_map_from_slope_map(float_image_t *GRD);
  /* Given a three-channel slope map {GRD}, returns a four-channel 
    normal map {NRM}, with the same size, such that channels {0..2} in
    column {x} and row {y} of {NRM} is the outward unit vector normal
    to the surface, computed from the corresponding slopes; and
    channel 3 of {NRM} is the weight (channel 2) of that pixel in {GRD}. */

#endif
