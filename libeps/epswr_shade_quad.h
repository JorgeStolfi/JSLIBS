/* Smooth shading of quadrilaterals. */
/* Last edited on 2022-10-20 06:53:49 by stolfi */

#ifndef epswr_shade_quad_H
#define epswr_shade_quad_H

#include <stdio.h>
#include <stdint.h>
#include <epswr.h>

void epswr_shade_quadrilateral
  ( epswr_figure_t *epsf,
    double x00, double y00, double R00, double G00, double B00,
    double x01, double y01, double R01, double G01, double B01,
    double x10, double y10, double R10, double G10, double B10,
    double x11, double y11, double R11, double G11, double B11,
    int32_t ns         /* Number of subdivisions. */
  );
  /* Fills a quadrilateral, given the corner coordinates
    {(x00,y00),(x01,y01),(x01,y10),(x11,y11)}, and their respective
    colors {(R00,G00,B00),.. (R11,G11,B11)}.
    
    The corners must be given in row-by-row order (NOT ccw order).
    
    If {ns == 0} the quadrilateral is painted solid with the average
    of the four colors.
    
    If {ns > 0} the quadrilateral is subdivided into {(ns+1)^2}
    smaller quadrilaterals by bilinear interpolation. Each of these is
    then painted with the solid color interpolated bilinearly 
    from the corner colors at its barycenter.
    
    If {ns < 0} the quadrilateral is painted with smooth bilinearly
    interpolated colors, using the using the Postscript 3.0 {shfill}
    operator with {ShadingType = 7} (2D bicubic Bézier patch) with
    the control points adjusted for bilinear interpolation. */

    
/* !!! Add routine for bicubic quadrangular shading given corners and 1st tensor derivatives. !!! */

#endif
