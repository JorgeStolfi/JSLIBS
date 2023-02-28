/* epswr_dev.h - draw EPS files in Device coordinates. */
/* Last edited on 2023-02-21 12:14:28 by stolfi */

#ifndef epswr_dev_H
#define epswr_dev_H
/*
  This module provides tools to simplify the generation of 
  Encapsulated Postscript (EPS) files containing graphics and text.

  An EPS file contains a single figure, of arbitrary size. It has a
  "%%BoundingBox" line at the beginning, no "%%Page" lines, and no
  explicit "showpage" command. It is usually meant to be included in
  other documents. Some printers and viewers may not be prepared to
  handle an isolated EPS file.  In what follows, we assume that an
  EPS file contains a single ``page,'' whose dimensions are 
  given by the "%%BoundingBox". 
 
  TYPICAL USAGE
  
  The typical usage for this interface is 
  
    { epswr_figure_t *eps = epswr_dev_new_figure(wr, ...);   }
    {   epswr_dev_segment(eps, ...)                          }
    {     ... other drawing commands ...                  }
    { epswr_dev_end_figure(eps);                             }
 */
  
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <epswr.h>
#include <bool.h>

epswr_figure_t *epswr_dev_new_figure
  ( FILE *wr,
    double hSize,   /* Figure width (in pt). */
    double vSize,   /* Figure height (in pt). */
    bool_t verbose  /* TRUE to print diagnostics. */
  );
  /* Writes to file {wr} the suitable preamble for an EPS figure, with
    the given dimensions. The figure's bounding box will be {[0 _ hSize]
    × [0 _ vSize]}. Both dimensions are in points (pt, 1/72 of an inch).
    
    Coordinate and dimension parameters for the functions in this 
    interface are expressed in a /Device coordinate system/, 
    whose origin at the bottom left corner of the bounding box.
    The first coordinate axis ({h}) is horizontal pointing from left to 
    right.  The second coordinate axis ({v}) is vertical, pointing up.
    Both coordinates are in points (pt).
  
    An {epswr_figure_t} object {eps} also specifies /Client coordinate system/
    that may be used by other interfaces.  This interface does not
    use that system. It is undefined on a newly created {epswr_figure_t}.
    and also after calling the plot window changing functions below. */

void epswr_dev_end_figure(epswr_figure_t *eps);
  /* Properly terminates the Encapsulated Postscript file, by writing
    the postamble, and closes the file. */

void epswr_dev_get_figure_size
  ( epswr_figure_t *eps,  /* Picture stream. */
    double *hSizeP,        /* OUT: Total width of figure (in pt). */
    double *vSizeP         /* OUT: Total height of figure (in pt). */
  );
  /* Sets {*hSizeP} and {*vSizeP} to the total width and
    height (in pt) of the given EPS figure. */
    
/* PLOT WINDOW */

/* An {epswr_figure_t} object {eps} keeps track of a rectangle on the
  Device coordinate plane, called the /plot window/. This rectangle is
  used in some operations like {epswr_dev_frame} and
  {epswr_dev_grid_cell}, and is a cliiping path for most other
  operations.
  
  The plot window is represented internally by its coordinate bounds
  {eps.hMin,eps.hMax,eps.vMin,eps.vMax}, in absolute Device
  coordinates, which are mirrored in Postscript variables of the same
  names. */

void epswr_dev_set_window
  ( epswr_figure_t *eps,
    double hMin, double hMax,
    double vMin, double vMax,
    bool_t relative
  );
  /* Sets the plot window (and default clip path) to the rectangle {R
    = [hMin _ hMax] × [vMin _ vMax]}, in Device coordinates (in pt).
    
    If {relative} is false, the coordinates are assumed to be
    relative to the lower left corner of the figure's bounding box.
    If {relative} is true, they are interpreted relative to the
    lower left corner of the current plot window.
    
    In any case, the procedure requires {hMin < hMax}, {vMin < vMax},
    and that the rectangle {R} be contained in the 
    figure's bounding box.  
    
    The procedure also causes the Client coordinate system to become
    undefined. */
  
 void epswr_dev_get_window
  ( epswr_figure_t *eps,
    double *hMinP, double *hMaxP,
    double *vMinP, double *vMaxP
  );
  /* Sets the variables {*hMinP,*hMaxP,*vMinP,*vMaxP} to the bounds
    of the current plot window (in pt, relative to the lower left
    corner of the EPS figure). */

void epswr_dev_shrink_window
  ( epswr_figure_t *eps, 
    double dhMin, double dhMax, 
    double dvMin, double dvMax
  );
  /* Shrinks the device plottng window by displacing the left edge by {dhMin} ,
    the right edge by {dhMax}, the bottom edge by {dvMin}, and the top edge 
    by {dvMax}, all inwards.  The displacements should be in Device 
    coordinates (pt).  Fails if the size of the plot window becomes zero
    or negative.  */

void epswr_dev_set_window_to_grid_cell
  ( epswr_figure_t *eps, 
    double hMin, double hMax, int32_t col, int32_t cols, 
    double vMin, double vMax, int32_t row, int32_t rows
  );
  /* Conceptually divides the the rectangle {[hMin _ hMax] × [vMin _ vmax]}
    (in absolute Device coordinates) into a grid of rectangular
    cells with {cols} columns, numbered {0..cols-1} from left to right,
    and {rows} rows, numbered {0..rows-1} from bottom to top;
    and sets the plot window to the cell in column {col} and row {row},
    using {epswr_dev_set_window}.
    
    The geometry and position of the text area
    (see {epswr_set_text_device_geometry} is not changed. */

void epswr_dev_set_pen
  ( epswr_figure_t *eps,
    double R, double G, double B,
    double pswidth,
    double psdashLength,
    double psdashSpace
  );
  /* */

void epswr_dev_segment
  ( epswr_figure_t *eps,
    double psxa, double psya,
    double psxb, double psyb
  );
  /* Draws segment from {(psxa,psya)} to {(psxb,psyb)}. */
  
void epswr_dev_curve
  ( epswr_figure_t *eps,
    double psxa, double psya,
    double psxb, double psyb,
    double psxc, double psyc,
    double psxd, double psyd
  );
  /* Draws a Bezier arc with given control points. */

void epswr_dev_coord_line
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double pspos
  );
  /* Draws a reference line PERPENDICULAR to the given axis 
    at the given coordinate value, extending across the whole
    plot window. */

void epswr_dev_coord_lines
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double psstart,
    double psstep
  );
  /* Draws equally spaced reference lines PERPENDICULAR to the given axis 
    at the given coordinate value, extending across the whole
    plot window.  One of the lines will have coordinate {psstart},
    and the lines will be {psstep} apart on that axis.  Only lines
    that interesect the plot window will be drawn. */

void epswr_dev_axis
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double pspos, 
    double pslo, 
    double pshi
  );
  /* Draws an axis line PARALLEL to the given axis, spanning the
    cordinate interval {pslo _ pshi}, on that axis, and positioned at
    coordinate {pspos} on the opposite axis. */

void epswr_dev_frame (epswr_figure_t *eps);
  /* Draws the outline of the current plot window. The outline 
    will be drawn half-inside, half-outside the area. */

/* FIGURE FILLING & DRAWING COMMANDS */

/* For all commands of this section, if the {fill} argument is true,
  the specified figure is painted with the current fill color; then,
  if {draw} is true, the outline of the figure is drawn with the
  current pen parameters.
  
  However, if the {R} component of the current fill color is
  negative, {±INF}, or {NAN}, the color is considered "invisible",
  and the shapes will not be filled, irrespective of the {fill}
  argument.
  
  The fill color is initially 50% gray. */

void epswr_dev_set_fill_color(epswr_figure_t *eps, double R, double G, double B);
  /* Defines the fill color for subsequent filling operations on
    {eps}. If {R} is negative,{±INF}, or {NAN}, the color is
    interpreted as "invisible". */

void epswr_dev_rectangle
  ( epswr_figure_t *eps,
    double psxlo, double psxhi,
    double psylo, double psyhi,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the given rectangle. */
  
void epswr_dev_triangle
  ( epswr_figure_t *eps,
    double psxa, double psya,
    double psxb, double psyb,
    double psxc, double psyc,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the triangle with corners {a,b,c}. */
  
void epswr_dev_quadrilateral
  ( epswr_figure_t *eps,
    double psx00, double psy00,
    double psx01, double psy01,
    double psx10, double psy10,
    double psx11, double psy11,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the quadrilateral with corners
   {(x00,y00),(x01,y01),(x10,y10),(x11,y11)}. IMPORTANT: The corners
   must be given in row-by-row order, not CCW order. */
  
/* !!! Polygons should do the right thing with dashes. !!! */
   
void epswr_dev_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double psx[], double psy[],
    int32_t n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  );
  /* Fills and/or outlines the polygon with corners 
   {(psx[0],psy[0]),.. (psx[n-1],psy[n-1])}.
   
   If {closed} is true the polygon is implicitly closed
   with the segment from {(psx[n-1],psy[n-1])} to {(psx[0],psy[0])}.
   Otherwise the procedure draws an open polgonal line.
   
   Filling is applied only if {closed} and {fill} are true.
   A point is considered to be inside the polygon iff its winding
   number is odd (when {evenOdd} is TRUE) or nonzero (when {evenOdd}
   is FALSE). */
  
void epswr_dev_rounded_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double psx[], double psy[],
    int32_t n,
    double psrad,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  );
  /* Same as {epswr_polygon}, but corners are rounded off with arcs
    of circles with radius {psrad}.
    
    The rounding of a corner will never remove an entire side,
    although the rounding of two adjacent corners may do so (and
    thus generate a straight joining line with negative length). If
    the rounding of a corner would require removing more than one
    entire side, the rounding of that corner is suppressed. In
    particular, this happens whenever one of the adjacent sides has
    zero length (i.e. two consecutive corners coincide) of the
    internal angle at the corner is 0 or 360 degrees.
    
    If {closed} is true the procedure assumes that the midpoint of
    the segment from point {n-1} to point 0 is not removed by
    rounding, and will start drawing from there. If {closed} is
    false the drawing starts at point 0 and ends at point {n-1}
    ignoring rounding. */
    
void epswr_dev_bezier_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double psx[], double psy[],
    int32_t n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  );
  /* Fills and/or outlines a curvilinear polygon defined by {n}
    Bézier arcs. Arc number {i} has control points
    {(psx[4*i+j],psy[4*i+j])}, for {i} in {0..n-1} and {j} in {0..3}. If
    the end of each arc is not the beginning of the next arc, adds a
    connecting straight line. A point is considered to be inside the
    polygon iff its winding number is odd (when {evenOdd} is TRUE)
    or nonzero (when {evenOdd} is FALSE). */
  
void epswr_dev_circle
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the circle centered at {(psxc,psyc)} with
    radius {psrad}. */

void epswr_dev_lune
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad, double pstilt,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the lune with given center, radius, and
    tilt (in radians CCW). (A "lune" is the intersection of two circles with the
    given radius, and with their centers spaced {psrad} apart.) */
  
void epswr_dev_slice
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad,
    double start, double stop,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the pie slice with given center, 
    radius, and angle range (in degrees). */

/* MARKER FILLING & DRAWING COMMANDS */

/* All commands in this section, the size of the figure is specified
  in millimeters, irrespecive of the current Client-to-Device scale.
  The specified size does not include the line width that is used
  when {draw} is true. */
    
void epswr_dev_dot
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad,
    bool_t fill, bool_t draw
  );
  /* Same as {epswr_circle}, except that the radius is in
    millimeters, not Client units. */
  
void epswr_dev_tic
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double psxc, double psyc, 
    double ticSize,
    double align 
  );
  /* Draws a tic mark (short segment) at coordinates {(psxc,psyc)}. The
    segment will be perpendicular to the given {axis} and its length
    will be {ticSize} millimeters, irrespective of the current scale)
    The segment will extend {align*ticSize} mm in the negative
    direction, and {(1-align)*ticSize} mm in the positive
    direction. */
  
void epswr_dev_cross
  ( epswr_figure_t *eps, 
    double psxc, double psyc, double psrad, bool_t diag,
    bool_t draw
  );
  /* Draws a cross centered at {(psxc,psyc)} with radius {psrad}. If
    {diag} is false the cross has vertical and horizontal branches;
    if {diag} is true the cross is rotated 45 degres. The radius is
    in millimeters, irrespective of the current scale. The mark has
    no interior. */
  
void epswr_dev_asterisk
  ( epswr_figure_t *eps, 
    double psxc, double psyc, double psrad,
    bool_t draw
  );
  /* Draws asn asterisk consisting of four crossed strokes at
    {(psxc,psyc)} with radius {psrad}. The {radius} is in millimeters,
    irrespective of the current scale. The mark has no interior. */

void epswr_dev_square
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or draws a square with center {psxc,psyc} and
    circum-radius {psrad}. The radius is in millimeters, irrespective
    of the current scale. */

void epswr_dev_diamond
  ( epswr_figure_t *eps, 
    double psxc, double psyc,
    double psxRad, double psyRad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or draws a diamond with center {psxc,psyc}, width {2*psxRad}
    and height {2*psyRad}. The last two parameters are in millimeters,
    irrespective of the current scale. */
    
void epswr_dev_arrowhead 
  ( epswr_figure_t *eps,
    double psxa, double psya, double psxb, double psyb,
    double pswidth, double pslength, 
    double fraction,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines a triangular head for an arrow with base
    at {a = (psxa,psya)} and tip at {b = (psxb,psyb)}. The head will have
    the specified {width} and {length} and its tip will be positioned at the
    given {fraction} of the way from {a} to {b}. */
    
/* GRID LINES AND GRID CELLS */

void epswr_dev_grid_lines(epswr_figure_t *eps, int32_t cols, int32_t rows);
  /* The plot window is implicitly divided into a grid of rectangular
    cells, with {cols} columns and {rows} rows. Draws all grid cell boundaries
    with the current pen and color. */

void epswr_dev_grid_cell
  ( epswr_figure_t *eps, 
    int32_t col, int32_t cols,
    int32_t row, int32_t rows,
    bool_t fill, bool_t draw
  );
  /* The plot window is implicitly divided into a grid of rectangular
    cells, with {cols} columns and {rows} rows. Fills and/or outlines cell 
    in column {col} and row {row} of that grid.  Cell {[0,0]}
    lies at the bottom left corner. */

/* LABELS */

void epswr_dev_set_label_font(epswr_figure_t *eps, const char *font, double size);
  /* Sets the name and point size of the font to be used by {epswr_label}. */

void epswr_dev_label
  ( epswr_figure_t *eps, 
    const char *text, 
    const char *strut, 
    double psx, double psy, 
    double rot,
    bool_t clipped,
    double hAlign, double vAlign,
    bool_t fill, bool_t draw
  );
  /* Prints {label} at point {(psx,psy)}, using the current label font
    and font size.
    
    Treats the outline of the characters as a path that is to be
    filled with the current fill color and/or stroked with the
    current pen.
    
    The parameter {hAlign} specifies which point of the string's
    bounding box will end up at {(psx,psy)}: 0.0 means the left side,
    1.0 means the right side. Other values of are
    interpolated/extrapolated, so that {hAlign = 0.5} the label will
    be horizontally centered at {(psx,psy)}.  
    
    The parameter {vAlign} controls the vertical alignment of the
    bounding box, with 0.0 meaning bottom and 1.0 meaning top. However
    the vertical extent of the box will be that of the {strut} string
    (which must not be empty) rather than that of {text}. Thus, for
    example, if {strut} is "x", then {vAlign=0} means the font's
    baseline and {vAlign=1} means the top of lowercase letters (the
    /x-height/). If {strut} is "Rg", then {vAlign=0} means the bottom of
    the "g" descender and {vAlign=1} means the top of uppercase letters.
    To get {vAlign} be relative to the bounding box of {text}, use
    {strut=text}.
    
    In any case, after the alignment is applied, the label will be
    rotated by {rot} degrees counterclockwise around the point
    {(psx,psy)}.
    
    If {clipped} is TRUE, the label will be clipped to the 
    current plot area.  Otherwise it may extend over the 
    whole figure.
    
    Beware that the label's bounding box is computed from the
    character outlines, not from the nominal character boxes. In
    particular, leading and trailing blanks are ignored. */

/* RUNNING TEXT */

void epswr_dev_set_text_geometry
  ( epswr_figure_t *eps, 
    double hMin, double hMax, 
    double vMin, double vMax,
    double rot
  );
  /* Sets the nominal text margins for {epswr_text} to the rectangle {[hMin_hMax] x
    [vMin x vMax]} in Device coordinates. The actual text will be
    rotated by {rot} degrees counterclockwise relative to the center
    of that rectangle. The next call to {epswr_dev_text} will write a line
    flush at the top of the rectangle.
    
    !!! Should have a {relative} parameter. !!! */

void epswr_dev_set_text_font(epswr_figure_t *eps, const char *font, double size);
  /* Sets the name and point size of the font to be used by {epswr_text}. */

void epswr_dev_text
  ( epswr_figure_t *eps, 
    const char *text, 
    bool_t clipped,
    double hAlign, 
    bool_t fill, bool_t draw
  );
  /* Writes one or more lines of string {text} after the lines already
    written. 
    
    Line breaks ("\n" or "\r\n" or just "\r") in the text are honored; a string
    that contains {k} line breaks generates {k+1} successive lines of
    text.  Note that an empty {text} generates one blank line,
    and a {text} that ends in with a line break generates a
    blank line after the preceding text.  In any case, the {text}
    MUST be followed by a zero byte ('\000').
    
    The position of the text is determined by the rectangle
    {[hMin_hMax] x [vMin x vMax]} defined in the last call to
    {epswr_set_text_geometry} and by the text lines
    already written since that call. Successive lines are written
    top to bottom in that rectangle, with vertical spacing
    determined by the current font's nominal size.
    
    More precisely, the top of the first line of {text} will be {dy}
    points below the top of the last line previously written, where
    {dy} is the current nominal text font size (as set by
    {epswr_set_text_font}); or at top of the text rectangle, if no lines
    have been written since the last call to
    {epswr_dev_set_text_geometry}.
    
    The parameter {hAlign} defines the alignment of each line of
    {txt} relative to the nominal margins {hMin} and {hMax}. If
    {hAlign} is 0.0, the left edge of the box will be at {hMin}. If
    {hAlign} is 1, the right edge of the box will be at {hMaz}.
    Other values are interpolated/extrapolated; so {hAlign=0.5}
    means that the box will be centered beween {hMin} and {hMax}.
    
    The position of each line of {txt} is determined without taking
    the rotation {rot} into account. The line is then written
    rotated by {rot} degrees counterclockwise around the center of
    the text rectangle.
    
    If {clipped} is TRUE, the text will be clipped to the 
    current plot area.  Otherwise it may extend over the 
    whole figure.
    
    !!! The vertical spacing {dy} should be an independent parameter? !!!
    */

/* MISCELLANEOUS */ 

void epswr_dev_set_verbose(epswr_figure_t *eps, const bool_t verbose);
  /* Sets an internal diagnostic flag of {eps},initially FALSE, to
    {verbose}. When that flag is set, some functions print diagnostic
    output to {stderr}. */

void epswr_dev_comment(epswr_figure_t *eps, const char *title);
  /* Writes a Postscript comment line to the underlying file.
    It is intended for documentation and debugging
    purposes only; it has no visible effect on the printed 
    page/figure. */

void epswr_dev_show_stack(epswr_figure_t *eps, int32_t code);
  /* Issues a command that, when the EPS file is interpreted
    by {gs} (Ghostscript) or similar software, will
    print the contents of the Postscript stack to {stderr}.
    
    The {code} is an arbitrary integer that can be used 
    to identify the instance of the operation; it will appear
    at the top of the printed stack, but will  be deleted
    after the printout. */

void epswr_dev_flush (epswr_figure_t *eps);
  /* Flushes the underlying file. */

/* HEADERS AND OTHER INTERNAL OPS */

void epswr_dev_add_font(const char *font, int32_t *nfontsP, char ***fontsP);
  /* Appends a copy of the given {font} name string to the font list
    {(**fontsP)[0..(*nfontsP)-1]}. */

void epswr_dev_write_window_set_cmds
  ( FILE *wr,
    double hMin, double hMax,
    double vMin, double vMax
  );
  /* Writes to {wr} some commands that set the Postscript
    variables {hMin,hMax,vMin,vMax} to the given values,
    and sets the clip path to that rectangle. */

void epswr_dev_write_grid_set_cmds(FILE *wr, int32_t hGridN, int32_t vGridN);
  /* Writes to {wr} some commands that set the Postscript
    variables {hGridN,vGridN} that define the cell grid */

void epswr_dev_write_label_font_set_cmds(FILE *wr, const char *font, double size);
  /* Writes to {wr} commands that create an instance of the font
    named {font}, scaled by {size} pt, and saves it into the Postscript
    variable {textFont}. */

void epswr_dev_write_text_font_set_cmds(FILE *wr, const char *font, double size);
  /* Writes to {wr} commands that create an instance of the font
    named {font}, scaled by {size} pt, and saves it into the  Postscript
    variable {labelFont}. */

void epswr_dev_write_ps_string(FILE *wr, const char *text);
  /* Writes the first line of the given text out to {wr},
    in Postscript form, properly quoting any special chars
    and parentheses.
    
    Printing ASCII characters are written verbatim. Tabs '\t' are
    replaced by two spaces. Backslashes and parentheses are escaped with
    a '\'. Intreprets a zero or line break character {'\000'=NUL,
    '\012'='\n'=LF, '\015'='\r'=CR}) as the end of the string. Other
    characters are written as '\' followed by the 3-digit octal
    value. */

void epswr_dev_write_fill_color_set_cmds(FILE *wr, double fc[]);
  /* Writes the Postoscript command that saves {fc[0..2]}
    as the RGB fill color for subsequent figure commands. */

#endif
