#ifndef neuromat_eeg_image_basis_H
#define neuromat_eeg_image_basis_H

/* NeuroMat tools to interpolate values over flat scalp image. */
/* Last edited on 2021-08-20 15:51:29 by stolfi */

#define _GNU_SOURCE
#include <r2.h>
#include <r3.h>
#include <float_image.h>

#include <neuromat_eeg.h>
#include <neuromat_image.h>
#include <neuromat_eeg_image.h>
    
typedef enum {
    neuromat_eeg_image_basis_type_SHEPARD0,
    neuromat_eeg_image_basis_type_SHEPARD1,
    neuromat_eeg_image_basis_type_RADIAL,
    neuromat_eeg_image_basis_type_VORONOI,
    neuromat_eeg_image_basis_type_CORREL
  } neuromat_eeg_image_basis_type_t;
#define neuromat_eeg_image_basis_type_LAST neuromat_eeg_image_basis_type_CORREL
  /* The avaliable basis types. */

float_image_t **neuromat_eeg_image_basis_make(int btype, int ne, float_image_t *msk, r2_t pos2D[], r3_t *srad, r2_t *ictr, r2_t *irad);
  /* Computes an interpolating image basis {bas[0..ne-1]} of type
    {btype} for electrodes with schematic positions {pos2D[0..ne-1]}.  
    
    An /interpolating image basis/ for {ne} electrodes is a list of 2D
    images {bas[0..ne-1]} of the scalp such that each image {bas[i]}
    has value 1 at the position of electrode {i} and value 0 at the
    positions of all other electrodes.
    
    Some basis functions require computing distances between points in the scalp.
    Since the actual positions of the electrodes on the scalp are unknown,
    the basis are computed with their /idealized 3D positions/ which are obtained from 
    their schematic 2D
    
    
    3D as if the scalp was
    scaled in all three axes so that it can be approximated by 
    the sphere 
    
    The procedure assumes that the 2D electrode positions {pos2D[0..ne-1]}
    are of the schematic 3D electrode
    positions on the unit sphere to the plane, so that the {+Z} hemisphere projects into the unit
    disk.
    
    
    The basis functions are actually defined on the unit sphere.
    
    The image {msk} must be a single-channel mask for the head outline
    (1 inside, 0 outside, possibly antialiased along the border). Each
    basis element will have the same size as {msk} and will be zero
    where {msk} is zero.
    
    In the images, the basis functions are projected to the plane
    stereographically as above, then stretched and translated so that
    the unit disk (the {+Z} hemisphere) becomes the ellipse with center
    {ctr} and main radii {rad.c[0],rad.c[1]}. The electrode positions,
    thus mapped, should normally be inside the head outline as defined
    by {msk}. */
    
/* SPECIFIC INTERPOLATING BASES */

/* The following procedures paint into the images {bas[0..ne-1]} an interpolating image
  basis for electrodes in the positions {pos2D[0..ne-1]} on the sphere. */

void neuromat_eeg_image_basis_fill_shepard_0(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[]);
  /* Uses the zero-order Shepard method. The basis is smooth and has the partition
    of unity property: all samples are between 0 and 1, and, for any
    pixel, the sum of all basis images at that pixel is 1.  The gradient
    of every element is zero at every electrode position, and the 
    same will be true for the interpolated image. */
    
void neuromat_eeg_image_basis_fill_shepard_1(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[]);
  /* Uses a first-order version of the Shepard method. The basis is smooth and its values sum to 1 at any pixel, 
    but may overshoot or undershoot the interval {[0_1]}.  It is an interpolating basis.  The gradient
    of every element at every electrode position depends on the positions of the nearby electrodes. */
    
void neuromat_eeg_image_basis_fill_radial(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[], double rho);
  /* The basis is a linear combination of gauss-sinc radial elements with
    ratio {sigma/tau} equal to {rho}. It usually has sample 
    values less than 0 and greater than 1. */
    
void neuromat_eeg_image_basis_fill_voronoi(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[]);
  /* The basis is essentially the Voronoi diagram of the electrodes in image form, after 
    expanding the axes so that the ellipse becomes a circle. The values are
    between 0 and 1 and form a partition of unity. */
     
void neuromat_eeg_image_basis_fill_correl(float_image_t *msk, int ne, float_image_t *bas[], r2_t *ctr, r2_t *rad, r3_t pos3D[], double rho);
  /* The basis is computed so that pairs of interpolated values have correlations that
    decay with their distance {d} as {1 - rho*d}. The value of {rho} must be between 0
    and {-2/dmax} where {dmax} is the diameter of the domain (max distance between two 
    scalp pixels). */
   
typedef double neuromat_image_basis_elem_t(int j, r3_t *p);
  /* A procedure that evaluates a generic basis element number {j} at an
    arbitrary point {p} of the unit sphere. */

void neuromat_eeg_image_basis_lagrangian_matrix(int ne, r3_t pos3D[], neuromat_image_basis_elem_t *elem, double *L);
  /* Given an arbitrary function basis {elem} with {ne} elements and 
    a list of electrode positions {pos3D[0..ne-1]}, returns in {L}
    the mixing matrix that turns {elem} into an interpolating basis at
    those points.
    
    The procedure assumes that {elem(j,p)} returns the value of element {j} of 
    of the given basis at the point {p} of the unit sphere,
    for all {j} in {0..ne-1}. It also assumes that {*L} is a preallocated vector with {ne^2} elements,
    representing a matrix with {ne} by {ne} elements,linearized by rows.
    
    On return, the desired interpolation basis {interp} is obtained by
    
      { interp(i,p)  = SUM{ L[i,j]*elem(j,p) : j \in 0..ne-1 } } */

#endif
