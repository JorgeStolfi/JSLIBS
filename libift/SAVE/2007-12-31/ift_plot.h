/* ift_plot.h - Postscript plotting of IFT path forests. */
/* Last edited on 2007-12-26 16:20:54 by stolfi */

#ifndef ift_plot_H
#define ift_plot_H

#include <ift.h>
#include <jspnm.h>
#include <jspnm_image.h>
#include <pswr.h>

/* 
  The following procedures assume that the image is 
  a grid with square cells of unit side whose
  lower left corner is at plot coordinates `(0,0)'. */

void ift_plot_pixel(
    PSStream *ps, 
    ImageGraph *G, 
    PixelIndex col, 
    PixelIndex row,
    double r, double g, double b,
    int outline
  );
  /* 
    Fills the specified pixel with the specified color.
    If `outline' is true, also draws its outline with the current line parameters. */

void ift_plot_node(
    PSStream *ps, 
    ImageGraph *G, 
    PixelIndex col, 
    PixelIndex row,
    double radius,
    double r, double g, double b
  );
  /* 
    Fills and draws a round node centered on the specified pixel,
    with the current line parameters. The radius is in millimeters. */

void ift_plot_arc(
    PSStream *ps, 
    ImageGraph *G, 
    PixelIndex col1, 
    PixelIndex row1,
    PixelIndex col2, 
    PixelIndex row2,
    int arrow
  );
  /* 
    Draws an arc from pixel (col1, row1) to pixel (col2, row2).
    If `arrow' is TRUE, also draws the arrowhead. */

void ift_plot_pixel_values(
    PSStream *ps, 
    ImageGraph *G,
    double whiten
  );
  /* 
    Paints the pixels of the image, mixed with white color 
    in the proportion `(1-whiten):whiten'. */

void ift_plot_forest_edges(
    PSStream *ps, 
    ImageGraph *G
  );
  /* 
    Plots the edges of the optimum-path forest defined
    by the `P' map in `G', using the current pen settings. */

void ift_plot_forest_nodes(
    PSStream *ps, 
    ImageGraph *G,
    int roots,
    double radius,
    double r, double g, double b
  );
  /* 
    Plots the root nodes (if `roots=1') or non-root nodes (if `roots = 0')
    of the optimum-path forest defined by the `P' map in `G', as
    circles of the specified radius, using the current fill color and
    pen settings. */

#endif

