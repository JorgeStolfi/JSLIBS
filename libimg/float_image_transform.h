#ifndef float_image_transform_H
#define float_image_transform_H

/* Tools for projective and barrel/pincushion image correction. */
/* Last edited on 2024-12-04 23:32:16 by stolfi */ 

#include <r2.h>
#include <r2_extra.h>
#include <r2x2.h>
#include <r3x3.h>
#include <bool.h>
#include <ix_reduce.h>
#include <float_image.h>

void float_image_transform_all
  ( float_image_t *iimg,     /* Input image. */
    ix_reduce_mode_t red,      /* Index reduction method. */ 
    r2_map_jacobian_t *map,  /* Output-to-input coordinate transformation. */
    float undef,             /* Sample value for undefined output pixels. */
    bool_t avg,              /* TRUE to average pixels, FALSE to add them. */
    int32_t order,               /* Interpolation order. */
    r2_pred_t *debugp,       /* Predicate, tells whether pixel should be debugged (NULL means "never"). */
    float_image_t *oimg      /* Output image. */
  );
  /* Applies the transformation {map} to the image {iimg} and
    stores the result in {oimg}. 
    
    The {map} procedure should compute the *inverse* of the desired
    image transformation. More precisely, the value of the transformed
    image {timg} at point {op} is assumed to be the value of {iimg} at
    point {ip}, where {ip} is computed by {ip=op; map(&ip,&J)}.
    
    The {order} parameter determines how the value of {iimg} at a
    given point is obtained from nearby samples. See
    {float_image_transform_interpolation_HELP_INFO} below.
    
    Point coordinates are relative to the standard {float_image}
    coordinate system for the relevant image. Namely, the first (X)
    coordinate increases with the column index, the second (Y)
    coordinate increases with the row index, pixels are squares with
    side 1, the image domain is the rectangle {[0 _ NX]×[0 _ NY]}, and
    the pixel in column 0, row 0 adjacent to the origin. Note that
    pixel centers have half-integer coordinates.
    
    If {avg} is true, each pixel of {oimg} is computed by a suitable
    average of the values of {timg} at several sampling points 
    over the pixel.  That is, the pixel of {oimg} is the 
    average of several values of {iimg} around the point {ip} corresponding
    to the pixel's center.  If {avg} is FALSE, the procedure computes a
    weighted integral, so as to approximately preserve the total
    integral of {iimg} over its domain. This option is appropriate,
    for instance, when {iimg} is a probability distribution, or a
    weight mask.
    
    For the averaging, the procedure assumes that the given {map} can
    be approximated by an affine (1st degree) map, which is determined
    by the output pixel center {op}, the corresponding point {ip} in
    the input domain, and the Jacobian matrix {J} returned by the
    {map} procedure. This approximation is assumed to be valid within
    a couple of pixels of the point {op}.
    
    If {map} returns {(NAN,NAN)} when applied to a pixel center
    {op}, the procedure assums that {op} has no corresponding point in
    the input domain, and sets all pixel's samples to the given
    {undef} value.

    The procedure implicitly extends the input image {ip} to infinity by
    replication of the pixels around its edges. This assumption is
    relevant when {map} returns a point {ip} that is not {(NAN,NAN)}
    but lies outside the domain of {iimg}.
    
    If the predicate {debug} is not null, the procedure calls {debug(p)}
    for every output pixel center {p}. If the result is TRUE, the
    procedure prints debugging information about the computation of that
    output pixel. */
    
void float_image_transform_sub
  ( float_image_t *iimg,     /* Input image. */
    ix_reduce_mode_t red,      /* Index reduction method. */ 
    r2_map_jacobian_t *map,  /* Output-to-input coordinate transformation. */
    float undef,             /* Sample value for undefined output pixels. */
    bool_t avg,              /* TRUE to average pixels, FALSE to add them. */
    int32_t order,               /* Interpolation order. */
    int32_t x0,                  /* First output image column. */
    int32_t y0,                  /* First output image row. */
    int32_t NX,                  /* Number of output image columns. */
    int32_t NY,                  /* Number of output image rows. */
    r2_pred_t *debugp,       /* Predicate, tells whether pixel should be debugged (NULL means "no"). */
    float_image_t *oimg      /* Output image. */
  );
  /* Same as {float_image_transform_all}, but copies only the 
    rectangular sub-image of {oimg} with pixel indices in the range
    {[x0..x0+NX-1]×[y0..y0+NY-1]}. */
  
void float_image_transform_get_pixel
  ( float_image_t *img,     /* Input image. */
    ix_reduce_mode_t red,     /* Index reduction method. */ 
    int32_t col,                /* Column index of output pixel. */
    int32_t row,                /* Row index of output pixel. */
    r2_map_jacobian_t *map, /* Output-to-input coordinate transformation. */
    float undef,            /* Sample value for undefined output pixels. */
    bool_t avg,             /* TRUE to average pixels, FALSE to add them. */
    int32_t order,              /* Interpolation order. */
    float f[],              /* Output pixel. */
    bool_t debug            /* If TRUE, prints debugging info. */
  );
  /* Computes the value of the pixel on column{col} and row {row} of a
    transformed image. The result is an average (if {avg} is true) or partition
    integral (if {avg} is false) of image values in the neighborhood
    of the the pixel's center. Uses {map} to find the corresponding
    points in the input image. Returns the result in {f[0..chns-1]}
    where {chns=img->sz[0]}.
    
    Pixels are interpolated using C0 bilinear interpolation (if
    {order} is 0) or C1 bicubic interpolation (if {order} is 1). Note
    that the latter may overshoot the interval {{0_1]}.
    
    If {map} returns {(NAN,NAN)}, the procedure sets all {f} samples
    to {undef}. */

#define float_image_transform_interpolation_HELP_INFO \
  "The input image value at a given point is obtained from nearby" \
  " samples by C0 bilinear interpolation (if the order is 0) or C1" \
  " bicubic interpolation (if the order is 1). Note that the latter may" \
  " overshoot the interval of the original samples."

/* PERSPECTIVE TRANSFORMATIONS

  The procedures below require a {3×3} projective map matrix {T2I}
  from some system of /true coordinates/ to /image coordinates/, that
  is, points in the domain of a given input image {iimg}. */
  
void float_image_transform_copy_persp_rectangle
  ( float_image_t *iimg,
    ix_reduce_mode_t red, /* Index reduction method. */ 
    double xlo,         /* Min X in true coords. */
    double xhi,         /* Max X in true coords. */
    double ylo,         /* Min Y in true coords. */
    double yhi,         /* Max Y in true coords. */
    r3x3_t *T2I,        /* Projective map from true coords to image coords. */
    float undef,        /* Defaut for undefined source pixels. */
    bool_t avg,         /* TRUE to compute average. */
    int32_t order,          /* Interpolation order to use. */
    int32_t x0,             /* First output image column. */
    int32_t y0,             /* First output image row. */
    int32_t NX,             /* Number of output image columns. */
    int32_t NY,             /* Number of output image rows. */
    r2_pred_t *debugp,  /* Predicate, tells whether pixel should be debugged (NULL means "no"). */
    float_image_t *oimg /* Output image. */
  );
  /* Copies a quadrilateral region {Q} from image {iimg} to a rectangular
    sub-image of image {oimg}.

    The receiving sub-image of {oimg} is defined by the lower indices
    {x0,y0}, and its dimensions in pixels {NX,NY}. The quadrilateral
    {Q} is the rectangle {R = [xlo_xhi] × [ylo_yhi]} in the in some
    /true coordinate system/, mapped by the {3×3} projective matrix
    {T2I}.

    The true rectangle {[xlo_xhi] × {ylo_yhi]} is mapped linearly to
    the rectangle {[x0_x0+NX} × {y0_y0+NY]} in the domain of {dst}.
    
    This procedure works no matter how small the {img} region is. */

#endif
