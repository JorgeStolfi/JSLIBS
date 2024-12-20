/* ift_plot.h - Postscript plotting of IFT path forests. */
/* Last edited on 2024-12-05 10:29:07 by stolfi */

#ifndef ift_plot_H
#define ift_plot_H

#include <stdio.h>

#include <bool.h>
#include <frgb.h>
#include <epswr.h>

#include <ift.h>

/* 
  The following procedures assume that the image is 
  a grid with square cells of unit side whose
  lower left corner is at plot coordinates {(0,0)}. */

void ift_plot_pixel
  ( epswr_figure_t *eps, 
    ift_graph_t *G, 
    ift_pixel_index_t col, 
    ift_pixel_index_t row,
    frgb_t *rgb,
    int outline
  );
  /* 
    Fills the specified pixel with the specified color.
    If {outline} is true, also draws its outline with the current line parameters. */

void ift_plot_node
  ( epswr_figure_t *eps, 
    ift_graph_t *G, 
    ift_pixel_index_t col, 
    ift_pixel_index_t row,
    double radius,
    frgb_t *rgb
  );
  /* 
    Fills and draws a round node centered on the specified pixel,
    with the current line parameters. The radius is in millimeters. */

void ift_plot_arc
  ( epswr_figure_t *eps, 
    ift_graph_t *G, 
    ift_pixel_index_t col1, 
    ift_pixel_index_t row1,
    ift_pixel_index_t col2, 
    ift_pixel_index_t row2,
    int arrow
  );
  /* 
    Draws an arc from pixel (col1, row1) to pixel (col2, row2).
    If {arrow} is TRUE, also draws the arrowhead. */

void ift_plot_pixel_values
  ( epswr_figure_t *eps, 
    ift_graph_t *G,
    frgb_t rgb[],
    double whiten
  );
  /* 
    Paints the pixels of {G} with the colors {rgb[0..G.nodes-1]}, mixed with white color 
    in the proportion {(1-whiten):whiten}. */

void ift_plot_forest_edges
  ( epswr_figure_t *eps, 
    ift_graph_t *G
  );
  /* 
    Plots the edges of the optimum-path forest defined
    by the {P} map in {G}, using the current pen settings. */

void ift_plot_forest_nodes
  ( epswr_figure_t *eps, 
    ift_graph_t *G,
    bool_t roots,
    double radius,
    frgb_t *rgb
  );
  /* 
    Plots the root nodes (if {roots=1}) or non-root nodes (if {roots = 0})
    of the optimum-path forest defined by the {P} map in {G}, as
    circles of the specified radius, using the current fill color and
    pen settings. */

#endif

