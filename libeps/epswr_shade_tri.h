/* Smooth shading of triangles. */
/* Last edited on 2009-08-25 23:27:36 by stolfi */

#ifndef epswr_shade_tri_H
#define epswr_shade_tri_H

#include <stdio.h>
#include <epswr.h>

void epswr_shade_triangle
  ( epswr_figure_t *epsf,
    double xa, double ya, double Ra, double Ga, double Ba,
    double xb, double yb, double Rb, double Gb, double Bb,
    double xc, double yc, double Rc, double Gc, double Bc,
    int ns         /* Number of subdivisions. */
  );
  /* Fills a triangle, given the corner coordinates
    {(xa,ya),(xb,yb),(xc,yc)} and the respective colors
    {(Ra,Ga,Ba),(Rg,Gg,Bg),(Rc,Gc,Bc)}.
    
    If {ns == 0} the triangle is painted solid with the average of the
    three colors.
    
    If {ns > 0} the triangle is subdivided into {(ns+1)^2} smaller
    triangles and each is painted solid with the color interpolated at
    its barycenter.
    
    If {ns < 0} the triangle is painted with Gouraud shading (linearly
    interpolated colors) using the Postscript 3.0 {shfill} operator
    with {ShadingType = 4} (triangular Gouraud shading). */
    
/* !!! Add routine for cubic triangular shading given corners and 1st derivatives. !!! */

#endif
