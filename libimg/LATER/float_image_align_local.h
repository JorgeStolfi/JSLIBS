#ifndef float_image_align_local_H
#define float_image_align_local_H

/* Tools for finding the optimum alignment of two or more images at given points. */
/* Last edited on 2024-12-05 10:29:21 by stolfi */ 

#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <r2x2.h>
#include <i2.h>
#include <float_image.h>

typedef bool_t float_image_align_local_report_t(int32_t ni, r2_t p[], double )

void float_image_align_local
  ( int32_t ni,               /* Number of images to align. */
    float_image_t *img[],     /* The images. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each image. */
    double crad,              /* Radius of comparison neighborhood (pixels). */
    r2x2_t L[],               /* Local linear mapping for each image, or {NULL}. */
    r2_t arad[],              /* Max alignment adjustment for each image. */
    bool_t quadopt,           /* True uses quadratic optimization, false uses exhausive enum. */ 
    double tol,               /* Desired precision. */
    float_image_align_local_report_t *rep  /* Reporting function, or {NULL}.
  );
  /* 
    The vectors {img,p,L,arad} must all have {ni} elements. The images {img[0..ni-1]}
    must all have the same number of channels.
    
    On input, each point {p[i]} should be within the domain of the image
    {img[i]}. The procedure tries to adjust those points so that the
    images are as similar as possible in the neighborhood of those
    points.  The radius of that neighborhood will be about 2 pixels.
    
    If {L} is not {NULL}, the neighborhood of each point {p[i]} is
    assumed to be deformed by the corresponding linear map {L[i]} before
    comparison. That is, for any {i1,i2}, and any small enough 2-vector
    {z}, point {p[i1]+z*L[i1]} of image {img[i1]} is assumed to
    coorespond to point {p[i2]+z*L[i2]} of image {img[i2]}. If {L} is
    {NULL} the {L[i]} are assumed to be the identity.
    
    For the details, see {r2_align.h}. In the nomenclature of that
    module, the list {p} is an alignment vector, and the changes made by
    the procedure are an adjustment vector in the ellipsoid {\RF}.

    That is, the output {p} will be a point of the basic ellipsoid {\RE}
    defined by the input vector {p} as the center and by {arad} as the
    radius vector.  In particular, if any coordinate {arad[i].c[j]} is
    zero, the coordinate {p[i].c[j]} will not be changed.
    
    Moreover, the adjustment made by the procedure will be /balanced/,
    meaning that, for each {j} in {0..1} the changes in {p[i].c[j]},
    summed over all {i}, are zero.  
    
    It follows that the procedure will be a no-op unless there is at least
    one {j} in {0..1} such that {arad[i].c[j]} is nonzero for at least two
    indices {i}. 
    
    The procedure will stop as soon as it believes that it has found an
    alignment vector that is within distance {tol} of the optimum.
    
    The procedure will try to minimize some measure of the mismatch of 
    the images in the neighborhoods defined by {p} and {L}
    
    If {rep} is not {NULL}, the procedure will call {rep(ni,q,F2)} every
    time it evaluates a candidae alignment {q[0..ni-1]}. The parameter
    {F2} will be a quadratic measure of the mismatch of the images in
    the neighborhood of those points. */

#endif
