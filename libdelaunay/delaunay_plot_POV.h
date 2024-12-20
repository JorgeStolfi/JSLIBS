/* Tools for plotting Delaunay and Voronoi diagrams */
/* Last edited on 2024-12-05 10:25:14 by stolfi */

#ifndef delaunay_plot_POV_H
#define delaunay_plot_POV_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <bool.h>
#include <delaunay.h>
#include <quad.h>

void delaunay_plot_POV_triangles(FILE *wr, quad_arc_t e, double height[], char* texture);
  /* Writes to {wr} a triangle mesh whose corners are the sites of the delaunay triangulation {e},
    with Z coordinates stored in {height} referenced by the index attribute of the sites.
    The mesh triangles are coloured with {texture}. */
 
void delaunay_plot_POV_skirt(FILE *wr, quad_arc_t e, double height[], char* texture);
  /* Same above, but writes only a skirt surrounding the triangle mesh */

#endif
