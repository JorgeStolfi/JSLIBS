/* Tools for plotting Delaunay and Voronoi diagrams */
/* Last edited on 2017-06-21 00:30:59 by stolfilocal */

#ifndef delaunay_plot_POV_H
#define delaunay_plot_POV_H

#include <quad.h>
#include <delaunay.h>
#include <bool.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

void delaunay_plot_POV_triangles(FILE *wr, quad_arc_t e, double height[], char* texture);
  /* Writes to {wr} a triangle mesh whose corners are the sites of the delaunay triangulation {e},
    with Z coordinates stored in {height} referenced by the index attribute of the sites.
    The mesh triangles are coloured with {texture}. */
 
void delaunay_plot_POV_skirt(FILE *wr, quad_arc_t e, double height[], char* texture);
  /* Same above, but writes only a skirt surrounding the triangle mesh */

#endif
