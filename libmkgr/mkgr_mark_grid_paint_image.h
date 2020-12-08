#ifndef mkgr_mark_grid_paint_image_H
#define mkgr_mark_grid_paint_image_H

/* mkgr_mark_grid_paint_image.h - functions to draw grid of marks as PS files. */
/* Last edited on 2020-11-29 18:28:28 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <r2.h>
#include <float_image.h>

#include <mkgr_mark.h>
#include <mkgr_mark_grid.h>

void mkgr_mark_grid_paint_image
  ( float_image_t *img,
    int32_t chns,           /* Number of channels to paint into. */
    int32_t ch[],           /* The channels, or {NULL} for the trivial ones. */
    mkgr_mark_grid_t *gr,   /* The list of marks. */
    double scale,           /* Scale factor from grid coordinates and sizes to pixels. */
    r2_t *org,              /* Pixel coordinates of grid origin, or {NULL}. */
    int32_t m               /* Subsampling factor. */
  );
  /* Draws the marks described in {gr} to the float image {img}. The
    marks are overlaid on top of the existing contents, with opaque
    paint but antialiasing along borders. The mark coordinates in {gr}
    are multiplied by {scale} and (if {org} is not {NULL}) shifted by
    {org}, and interpreted as displacements from the low corner of the
    image's domain, in pixel units.
    
    The parameter {m} specifies the subsampling order for antialiasing,
    as explained in {float_image_paint.h}.
  
    The image is assumed to cover the rectangle {[0 _ NX]x[0 _ NY]} of the 
    coordinate plane, with the lower corner at {(0,0)}; where {NX} and {NY} are
    the column and row counts. 
    
    The R, G, and B components of the marks' colors are painted in
    channels {ch[0..chns-1]}, repeating cyclically if {chns} is more
    than 3. If {chns} is zero, the procedure is a no-op. If {ch} is
    {NULL}, assumes {ch[i]=i} for all {i} */

#endif
