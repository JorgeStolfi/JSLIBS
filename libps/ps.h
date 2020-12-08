/* OBSOLETE - Routines for creating PostScript files with graphics and text. */
/* Last edited on 2007-12-26 14:53:52 by stolfi */

#ifndef ps_H
#define ps_H

#include <stdio.h>

/* THIS INTERFACE IS OBSOLETE AND IS PROVIDED ONLY FOR THE USE OF
  LEGACY PROGRAMS.  NEW PROGRAMS SHOULD USE {pswr.h} INSTEAD.

  There are two kinds of Postscript (PS) files, `encapsulated'
  (EPS) and `non-encapsulated' (non-EPS).
  
  An EPS file contains a single figure, of arbitrary size. It has a
  "%%BoundingBox" line at the beginning, no "%%Page" lines, and no
  explicit "show" command. It is usually meant to be included in other
  documents, and may confuse some printers.
  
  A non-EPS file contains a complete document, that uses the entire
  area of the page, and possibly multiple pages It has "%%Page"
  comments, explicit "show" commands at the end of each page, and no
  "%%BoundingBox" comment. It is usually meant to be printed by
  itself, and cannot be easily included in other documents. */

/*** ENCAPSULATED FIGURES ***/

void ps_begin_figure 
  ( FILE *psFile,
    double hmin, double hmax,
    double vmin, double vmax
  );
  /* Initializes an Encapsulated PostScript file, by writing to it a
    suitable preamble.
    
    The procedure writes the file's preamble, which defines some
    auxiliary Postscript operators and constants, sets up the caption
    font, etc. 
    
    The procedure also writes the parameters {hmin,hmax,vmin,vmax} (in
    pt) as the file's "%%BoundingBox" comment. */

void ps_end_figure (FILE *psFile);
  /* Finalizes an Encapsulated PostScript figure.
    The file is flushed but is left open. */
   
/*** NON-ENCAPSULATED MULTIPAGE DOCUMENTS ***/

void ps_begin_document (FILE *psFile, const char *paperSize);
  /* Initializes {psFile} to contain a non-encapsulated, possibly
    multi-page Postscript document.
    
    The procedure writes the file preamble, which defines some
    auxiliary Postscript operators and constants, defines the 
    caption font, etc.

    The parameter {paperSize} indicates the paper size of the
    document. Valid values are: "a1" "a2" "a3" "a4" "letter".
    Other values may not be recognized by the printer.
    [R. Chencarek dec/2002]
    
    The client must call {ps_begin_page} and {ps_end_page} around each
    page of the document, and {ps_set_window} before each separate
    drawing on the same page. */

void ps_begin_page(FILE *psFile, const char *page);
  /* Starts a new page of a PS document.  
  
    The {page} string is the logical page number, shown at the bottom
    of the page and as the first argument of the "%%Page" directive.
    If NULL or "", the page number is not printed, and defaults to the
    sequential page number, counting from 1. */

void ps_end_page(FILE *psFile);
  /* Finalizes a page: Writes page trailer line, etc. */

void ps_end_document (FILE *psFile);
  /* Finalizes a multipage document. The file is flushed but is left open. */
 
void ps_get_paper_dimensions(const char *paperSize, double *xpt, double *ypt);
  /* Sets {*xpt} and {*ypt} to the dimensions of the specified paper type,
    in points. Knows about US sizes "letter", "ledger", "tabloid",
    "legal", "executive", and the ISO "A" sizes (from "4A0" to
    "A10"). */

/* COORDINATE SETUP */
  
void ps_set_window
  ( FILE *psFile,
    double xmin, double xmax,
    double ymin, double ymax,

    double hmin, double hmax,
    double vmin, double vmax,

    int xn, int yn
  );  
  /* Sets up the coordinate system and clip path for a new drawing
    inside the figure (EPS file) or the current page (non-EPS file).

    Client coordinates will range over {[xmin _ xmax] × [ymin _ ymax]}.
    The nominal plotting area is {[hmin _ hmax] × [vmin _ vmax]}, in pt,
    relative to the lower left corner of the figure (EPS) or the
    current page (non-EPS). The plotting scales {dh/dx} and {dv/dy} must
    be equal.

    The plotting area is divided implicitly into a grid of {xn} by {yn}
    rectangular "cells", which are used by some of the plotting commands
    below. */

/*** DRAWING COMMANDS ***/

void ps_set_pen 
  ( FILE *psFile,
    double R, double G, double B,
    double width,
    double dashlength,
    double dashspace
  );
  /* Sets pen parameters and ink color for line drawing.
    Dimensions are in *millimeters* */

void ps_draw_segment
  ( FILE *psFile,
    double xa, double ya,
    double xb, double yb
  );
  /* Draws segment from {(xa,ya)} to {(xb,yb)} with current pen and color */

void ps_draw_tic
  ( FILE *psFile, 
    char axis, 
    double xc, double yc, 
    double ticsz,
    double align 
  );
  /* Draws a tic mark (short segment) at coordinates {(xc,yc)}. The
    segment will be perpendicular to the given {axis} ('X' or 'Y') and
    its length will be {ticsz} millimeters, irrespective of the
    current scale) The segment will extend {align*ticsz} mm in the
    negative direction, and {(1-align)*ticsz} mm in the positive
    direction. */
    
void ps_draw_curve
  ( FILE *psFile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double xd, double yd
  );
  /* Draws a Bezier arc with given control points, using the current
    pen and color. */
  
void ps_draw_coord_line (FILE *psFile, char axis, double coord);
  /* Draws a reference line perpendicular to the given axis 
    at the given coordinate value. */

void ps_draw_grid_lines(FILE *psFile);
  /* Draws the cell boundaries with the current pen and color. */

void ps_draw_frame (FILE *psFile);
  /* Draws a frame around the plotting area. (The frame 
    will extend half a line width outside the nominal bounding box.) */

void ps_draw_rectangle
  ( FILE *psFile,
    double xlo, double xhi,
    double ylo, double yhi
  );
  /* Draws the outline of the given rectangle using the current pen. */

void ps_fill_rectangle
  ( FILE *psFile,
    double xlo, double xhi,
    double ylo, double yhi,
    double R, double G, double B
  );
  /* Fills given rectangle with given color */

void ps_fill_and_draw_rectangle
  ( FILE *psFile,
    double xlo, double xhi,
    double ylo, double yhi,
    double R, double G, double B
  );
  /* Fills rectangle with given color, then draws its outline with
    current pen. */

void ps_draw_circle
  ( FILE *psFile,
    double xc, double yc, double radius
  );
  /* Draws the circle with given center and radius, using the current
    pen and color. */
  
void ps_fill_circle
  ( FILE *psFile,
    double xc, double yc, double radius,
    double R, double G, double B
  );
  /* Fills the circle with given center and radius, using the 
    given color. */
 
void ps_fill_and_draw_circle
  ( FILE *psFile,
    double xc, double yc, double radius,
    double R, double G, double B
  );
  /* Fills the circle with given center and radius, using the given
    color, then draws its outline, using the current pen and color. */
  
void ps_draw_dot
  ( FILE *psFile,
    double xc, double yc, double radius
  );
  /* Same as ps_draw_circle, but the radius is in millimeters,
    irrespective of the current scale. */
  
void ps_fill_dot
  ( FILE *psFile,
    double xc, double yc, double radius,
    double R, double G, double B
  );
  /* Same as ps_fill_circle, but the radius is in millimeters,
    irrespective of the current scale. */
 
void ps_fill_and_draw_dot
  ( FILE *psFile,
    double xc, double yc, double radius,
    double R, double G, double B
  );
  /* Same as ps_fill_and_draw_circle, but the radius is in
    millimeters, irrespective of the current scale. */

void ps_fill_and_draw_lune
  ( FILE *psFile,
    double xc, double yc, double radius, double tilt,
    double R, double G, double B
  );
  /* Fills the lune with given center, radius, and tilt, using the
    given color, then draws its outline, using the current pen and
    color. */
  
void ps_fill_and_draw_slice
  ( FILE *psFile,
    double xc, double yc, double radius, double start, double stop,
    double R, double G, double B
  );
  /* Fills the pie slice with given center, radius, and angle range
    (in degrees), using the given color, then draws its outline, 
    using the current pen and color. */
    
void ps_draw_polygon
  ( FILE *psFile,
    double x[], double y[],
    int n
  );
  /* Draws the contour of the polygon {(x[1],y[1]),.. (x[n],y[n]} 
    using the current pen and color. */
  
void ps_fill_polygon
  ( FILE *psFile,
    double x[], double y[],
    int n,
    double R, double G, double B
  );
  /* Fills the polygon {(x[1],y[1]),.. (x[n],y[n])} with the given
    color. */
  
void ps_fill_and_draw_polygon
  ( FILE *psFile,
    double x[], double y[],
    int n,
    double R, double G, double B
  );
  /* Fills the polygon {(x[1],y[1]),.. (x[n],y[n])} with the given color,
    then draws its contour using the current pen and color. */
  
void ps_fill_triangle
  ( FILE *psFile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double R, double G, double B
  );
  /* Fills triangle {a,b,c} with given color. */

void ps_fill_grid_cell(FILE *psFile, int xi, int yi, double R, double G, double B);
  /* Fills the given cell of the current cell grid with the given
    color. */

/* TEXT */

void ps_set_label_font(FILE *psFile, const char *font, double size);
  /* Sets the name and point size of the font to be used by
    {ps_put_label}. */

void ps_put_label
  ( FILE *psFile, 
    const char *text, 
    double x, double y, 
    double xalign, double yalign
  );
  /* Prints {label} at point {(x,y)}, using the current label font
    size. The parameter {xalign} (resp. {yalign}) specifies which
    point of the string's bounding box will end up at {(x,y)}: 0.0
    means the left (resp. bottom) side, 1.0 means the right (resp.
    top) side. Use (0.5, 0.5) to get the box centered at {(x,y}). */

void ps_add_caption(FILE *psFile, const char *txt);
  /* Adds a caption text below the drawing, *outside* the nominal
    bounding box. */

void ps_comment(FILE *psFile, const char *title);
  /* Writes a comment line to the given file. */

/*** LOW-LEVEL HACKS ***/

void ps_put_text(FILE *psFile, const char *text, const char *newline);
  /* Writes a text string to {psFile}, in Postscript form, properly
    quoting any special chars and parentheses. Replaces any embedded
    '\n' by the given {newline} string --- which could be, for
    example, ") show mynl (". */

#endif
