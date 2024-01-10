/* Common routines for all kinds of FOI test plots (EPS version) */

#ifndef PLTE_H
#define PLTE_H

#include <stdio.h>

void plte_begin_figure (
    FILE *psfile,
    double xmin, double xmax,
    double ymin, double ymax,
    int xn, int yn
  );
  /* Initializes an EPS figure file: writes page header line, */
  /* sets coordinate system, clip path, */
  /* defines new Postscript operators and constants, etc. */
  /* The plotting area is divided implicitly into a grid */
  /* of /xn/ by /yn/ rectangular "cells". */

void plte_begin_section (FILE *psfile, char *title);
  /* Starts a new section of a plot. The title is a comment */

void plte_set_pen (
    FILE *psfile,
    double gray,
    double width,
    double dashlength,
    double dashspace
  );
  /* Sets pen parameters and ink color for line drawing. */
  /* Dimensions are in mm */

void plte_draw_segment(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb
  );
  /* Draws segment from (xa,ya) to (xb,yb) with current pen and gray */

void plte_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi
  );
  /* Draws the outline of thegiven rectangle using the current pen. */

void plte_fill_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray
  );
  /* Fills given rectangle with given gray color */

void plte_fill_and_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double gray
  );
  /* Fills rectangle with given gray, then */
  /* draws its outline with current pen. */

void plte_fill_triangle(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double gray
  );
  /* Fills triangle /abc/ with given gray level. */

void plte_fill_grid_cell(FILE *psfile, int xi, int yi, double gray);
  /* Fills the given cell of the current cell grid with */
  /* the given gray level. */

void plte_draw_coord_line (FILE *psfile, char axis, double coord);
  /* Draws a reference line perpendicular to the given axis */
  /* at the given coordinate value. */

void plte_draw_grid_lines(FILE *psfile);
  /* Draws the grid lines with the current pen and gray level. */

void plte_end_section (FILE *psfile);
  /* Ends a section of a plot. */

void plte_draw_frame (FILE *psfile);
  /* Draws a frame around the plotting area */

void plte_end_figure (FILE *psfile);
  /* Finalizes figure file. */

#endif
