/* Plotting {pst_gr_t} structures. */ 
/* Last edited on 2025-03-14 06:38:01 by stolfi */

/* Created by Rafael F. V. Saracchini */

#ifndef pst_gr_plot_H
#define pst_gr_plot_H

#include <stdio.h>
#include <stdint.h>

#include <pst_gr.h>
#include <bool.h>
#include <epswr.h>

void pst_gr_plot(epswr_figure_t *eps, pst_gr_t* gr, double fontSize, double vertexRadius);
  /* Plots the graph {gr} into the Encapsulated Postcript figure {eps}.
    
    If the {fontSize} is positive, also plots the vertex index besides each vertex,
    and the edge label besides each edge, with a font with nominal height {fontSize}.
    
    Vertices are plotted as circles of radius {vertexRadius} filled with
    a color that depends on the {.vmark} field: white for 0, red for 1,
    yellow for 2, green for 3, blue for 4.
    
    The vertex coordinates and {vertexRadius} are assumed to be in mm.
    The font size is assumed to be in pt. */

void pst_gr_plot_named(char *fname, pst_gr_t *gr, double fontSize, double vertexRadius);
  /* Creates a file called "fname" with the Encapsulated Postcript plot of the graph {gr}.
    The figure size will be large enough to include all vertices and labels.
    The plot istself is created with {pst_gr_plot}. */

void pst_gr_path_plot(epswr_figure_t *eps, r2_t *po, pst_gr_path_t P, r2_t *pd);
  /* Plots into {eps} the path from {po} through vertices {P.v[0..P.n-1]} to {pd},
    using the current pen settings. */

#endif


