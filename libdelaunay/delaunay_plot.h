/* Tools for plotting Delaunay and Voronoi diagrams */
/* Last edited on 2017-06-21 01:08:31 by stolfilocal */

#ifndef delaunay_plot_H
#define delaunay_plot_H

#include <quad.h>
#include <delaunay.h>
#include <bool.h>
#include <i3.h>

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

void plot_delaunay (quad_arc_t e, delaunay_site_t *st, int nsites, char *prefix, bool_t eps);
  /* Writes the plot file named "{prefix}-doc.ps" or "{prefix}-{page}.eps". */
  
void draw_delaunay_edge (quad_arc_t e, void *closure);
void draw_voronoi_edge (quad_arc_t e, void *closure);

#endif
