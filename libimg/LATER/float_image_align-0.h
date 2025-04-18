#ifndef float_image_align_H
#define float_image_align_H

/* Tools for translational alignment between float-valued images. */
/* Last edited on 2024-12-04 23:19:53 by stolfi */ 

#include <bool.h>
#include <float_image.h>

void fial_optimal_shifts
  ( int32_t nims, 
    float_image_t *iml, 
    double dxmax, 
    double dymax, 
    double epsilon,
    double rx, 
    double ry, 
    double *dx, 
    double *dy
  );
  /* Takes a list {iml[0..nims-1]} of float-valued images, and
    computes a joint alignment for them with maximal quality.
    
    A joint alignment is a displacement {(dx[i],dy[i])} for each image
    {iml[i]}.  The quality of an alignment is measured by comparing
    each shifted image to the arithmetic average of all shifted images.
    
    The procedure will consider displacments in the range 
    {[-dxmax _ +dxmax]} for {dx} and {[-dymax _ +dymax]} for {dy}.
    
    The displacements will have nominal accuracy {epsilon}. The
    shifted images are compared with a Gaussian weight window 
    with deviations {rx} and {ry} along the X and Y axis, respectively.
    
    This procedure is best called with images, whose width is about
    {6*rx + 2*dxmax} and whose height is about {6*ry + 2*dymax}. */

void fial_multiscale_optimal_shifts
  ( int32_t nims, 
    float_image_t *iml, 
    double dxmax, 
    double dymax, 
    double epsilon, 
    double rx, 
    double ry, 
    int32_t mult, 
    alignment_t *opt, 
    int32_t *noptP
  );
  /* Takes a list {iml[0..nims-1]} of float-valued images, and
    computes up to {mult} distinct good joint alignments {opt[0..*noptP-1]}.
    
    The remaining parameters are as specified for {fial_optimal_shifts}. */

#endif
