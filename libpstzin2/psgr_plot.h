/* Plotting {psgr_t} structures. */ 
/* Last edited on 2025-04-27 11:21:01 by stolfi */

/* Created by Rafael F. V. Saracchini */

#ifndef psgr_plot_H
#define psgr_plot_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <epswr.h>
#include <r2.h>
#include <i2.h>

#include <psgr_types.h>
#include <psgr.h>

void psgr_plot(epswr_figure_t *eps, psgr_t* gr, double pixelSize, double fontSize, double vertexRadius);
  /* Plots the graph {gr} into the Encapsulated Postcript figure {eps}.
    
    If the {fontSize} is positive, also plots the vertex index besides each vertex,
    and the edge label besides each edge, with a font with nominal height {fontSize}.
    
    Vertices are plotted as circles of radius {vertexRadius} filled with
    a color that depends on the {.vmark} field: white for 0, red for 1,
    yellow for 2, green for 3, blue for 4.
    
    The {pixelSize} argument should be the desired size (in mm) of the
    pixels on the plot. Vertex pixel indices are converted to plot
    coordinates by multiplying {pixelSize} into them. The {vertexRadius}
    is assumed to be in mm. The font size is assumed to be in pt. */

void psgr_plot_named(char *fname, psgr_t *gr, double pixelSize, double fontSize, double vertexRadius);
  /* Creates a file called "fname" with the Encapsulated Postcript plot of the graph {gr}.
    The figure size will be large enough to include all vertices and labels.
    The plot istself is created with {psgr_plot}. */

void psgr_plot_path(epswr_figure_t *eps, double pixelSize, i2_t *ipo, psgr_path_t P, i2_t *ipd);
  /* Plots into {eps} the path from point {ipo} through nodes {P.v[0..P.n-1]} to point {ipd},
    using the current pen settings.
    
    The {X} and {Y} components of the arguments {ipo,ipd} and of the
    nodes {P.v[0..P.n-1]} are assumed to be in units of image pixels. In
    the plot, one pixel spacing will correpond to {pixelSize} mm. 
    
    If {P.reverse} is true, the nodes are taken in the reverse order
    {P.v[P.n-1]} to {P[0]}.  In any case the polygonal line will start
    at {ipo} and end at {ipd}. */

r2_t psgr_plot_coords_from_pixel(i2_t ip, double pixelSize);
  /* Converts the integer pixel indices {ip} to float coordinates and
    scales them by {pixelSize}.  The pixel indices must be non-negative. */

r2_t psgr_plot_coords_from_path(int32_t k, i2_t *ipo, psgr_path_t P, i2_t *ipd, double pixelSize);
  /* Returns the plot coordinates (in mm) of a vertex of the polygonal
    line from point {ipd} through path {P} to point {ip2}. The index {k}
    should be in the range {-1..P.n}, where {-1} signifies {ipo}, {P.n}
    signifies {ipd}, and any value in {0..P.n-1} signifies either {P.v[k]}
    (if {P.reverse} is false) or {P.v[P.n-1-k]} (if {P.reverse} is true). In any
    case, the resulting coordinates are converted from pixels to mm by
    multiplying with {pixelSize}. */

#endif


