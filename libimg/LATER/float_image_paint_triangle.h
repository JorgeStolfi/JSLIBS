#ifndef float_image_paint_triangle_H
#define float_image_paint_triangle_H

/* Tools for drawing into float images. */
/* Last edited on 2021-08-28 23:27:45 by stolfi */ 

#include <bool.h>
#include <float_image.h>

double float_image_paint_triangle
  ( float_image_t *A, 
    int c,           /* Channel. */
    double xa, double ya,    /* Triangle corner. */
    double xb, double yb,    /* Triangle corner. */
    double xc, double ycc,   /* Triangle corner. */
    double hwd,              /* Radius of pen tip. */
    float vfill,             /* Ink value for filling. */                              
    float vdraw,             /* Ink value for stroking. */  
    int m                    /* Subsampling parameter. */
  );
  /* Draws into channel {c} of image {A} a triangle with corners {(xa,ya),(xb,yb),(xc,yc)}.
    The sign of {hwd} is ignored. */

#endif
