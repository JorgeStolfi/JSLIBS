#ifndef float_image_average_H
#define float_image_average_H

/* Averaging float images inside extended regions. */
/* Last edited on 2024-11-23 05:55:28 by stolfi */ 

#include <bool.h>
#include <ix_reduce.h>
#include <r2.h>
#include <r2x2.h>
#include <r3x3.h>
#include <interval.h>
#include <float_image.h>
   
void float_image_average_parallelogram
  ( float_image_t *img, 
    ix_reduce_mode_t red,
    r2_t *p, 
    r2x2_t *J,
    bool_t avg,
    int order,
    float f[],
    bool_t debug
  );
  /* Evaluates the image {img} at the point {p = (x,y)}. The value is an
    average or integral of the values of samples of {img} around {p},
    using a pixel sampling kernel whose domain has been transformed by
    the Jacobian {J}. The result is returned in {v[0..chns-1]} where
    {chns=img->sz[0]}.
    
    The kernel weights are such that if {p} is stepped over a regular
    unit grid linearly mapped by the matrix {J}, the resulting kernels 
    form a partition of unity over the sample points. 
    
    If {avg} is true, the result is a weighted average of the samples in
    the kernel domain. This is appropriate, for example, it {img} is a
    digital photo.
    
    If {avg} is false, the result is simply the weighted sum
    of the sample values. This is appropriate, for example, it {img} is a 
    density distribution.
    
    The subsamples used to compute the result are interpolated in {img}
    using C0 bilinear interpolation (if {order} is 0) or C1 bicubic
    interpolation (if {order} is 1). Note that the latter may overshoot
    the interval {{0_1]}. 
    
    Any pixel index that falls outside its valid range {0..N-1} in {img} is
    reduced with {ix_reduce(i,N,red)}; if the result is {-1} the
    pixels are omitted from the average. */
 
void float_image_average_persp_rectangle
  ( float_image_t *img, /* Image to sample. */
    ix_reduce_mode_t red, /* Index reduction method. */ 
    interval_t tbox[],  /* Rectangle in true coordinates. */
    r3x3_t *T2I,        /* Projective map from true coords to image coords. */
    double mrg,         /* Safety border width (pixels). */
    double avg[]        /* (OUT) Pixel average. */
  );
  /* Computes the average {avg[0..NC-1]} of the pixels in image {img}
    within a quadrilateral region {Q}, where {NC=img->sz[0]}. The
    quadrilateral is defined by a recangle {tbox[0] × tbox[1]} in true
    coordinates, and the true-to-image matrix {T2I}. 
    
    The average considers only pixels that are inside the
    quadrilateral {Q} and at least {mrg} pixels away from the boundary
    of the quadrilatera; if there are no such pixels, sets {val} to
    {NAN}s.
    
    Any pixel index that falls outside its valid range {0..N-1} in {img} is
    reduced with {ix_reduce(i,N,red)}; if the result is {-1} the
    pixels are omitted from the average. */

void float_image_average_persp_disk
  ( float_image_t *img, /* Input image. */
    ix_reduce_mode_t red, /* Index reduction method. */ 
    r2_t *ctr,          /* Disk center in true coordinates. */
    double rad,         /* Disk radius in true coordinates. */
    r3x3_t *T2I,        /* True-to-image projective map matrix. */
    double mrg,         /* Safety border width (pixels). */
    float avg[]         /* (OUT) average disk color. */
  );
  /* Computes the average {avg[0..NC-1]} of the pixels in image {img}
    within an ellipse {E}, where {NC=img->sz[0]}. The ellipse is defined
    by the center {ctr} and radius {rad} of a disk in true
    coordinates, and the true-to-image matrix {T2I}.

    The average considers only pixels that are inside the ellipse {E}
    and at least {mrg} pixels away from its border in the image. The
    procedure returns all {NAN}s in {avg} if there are no such pixels
    (that is, if the ellipse is too small of too narrow).
    
    Any pixel index that falls outside its valid range {0..N-1} is
    reduced with {ix_reduce(i,N,red)}; if the result is {-1} the
    pixels are omitted from the average. */

#endif
