/* Tools for plotting Delaunay and Voronoi diagrams */
/* Last edited on 2024-12-05 10:25:10 by stolfi */

#ifndef delaunay_plot_H
#define delaunay_plot_H

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include <i3.h>
#include <bool.h>
#include <epswr.h>
#include <quad.h>

#include <delaunay.h>

/* In the following procedures, the sites of the Delaunay/Voronoi diagram are {st[0..nsites-1]}. */

void delaunay_plot_diagram(epswr_figure_t *eps, quad_arc_t e, delaunay_site_t st[], int nsites);
  /* Draws into {eps} a Delaunay diagram, given one of its quad-arcs {e} and the 
    table {st} with the site coordinates. */
  
void delaunay_plot_edge(epswr_figure_t *eps, quad_arc_t e, delaunay_site_t st[], int nsites);
  /* Deaws the egde {e} of the Delaunay diagram. */
  
void delaunay_plot_voronoi_edge(epswr_figure_t *eps, quad_arc_t e, delaunay_site_t st[], int nsites);
  /* Draws the edge {e} of the Voronoi diagram. */
  
void delaunay_plot_site(epswr_figure_t *eps, delaunay_site_t *si);
  /* Draws the site {si} of the Voronoi/Delaunay diagram. */
  
epswr_figure_t *delaunay_plot_new_figure
  ( char *prefix,
    char *tag,
    int32_t capLines, 
    delaunay_site_t st[], 
    int nsites
  );
  /* Creates a new EPS file called "{prefix}_{NNNNN}_{tag}.eps" suitable
     for a Delaunay/Voronoi diagram with the sites {st[0..nsites-1]};
     where {NNNNN} is {nsites} formatted as "%05d".
     
     The figure will have space for {capLines} lines of caption text
     below the plot area proper. The text area is defined to be that
     space, and the font is set to "Courier" 10 pt, so caption lines can
     be added by calling {epswr_text(eps,text,FALSE,hAlign,TRUE,FALSE)}. */

#endif
