/* Tools for generating Encapsulated Postscript graphics files. */
/* Last edited on 2020-10-27 18:49:36 by jstolfi */

#ifndef epswr_H
#define epswr_H

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
  
    { epswr_figure_t *eps = epswr_new_figure(wr, ...);   }
    {   epswr_segment(eps, ...)                          }
    {     ... other drawing commands ...                 }
    { epswr_end_figure(eps);                             }
 */
  
#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>

/* EPS FIGURE OBJECT */

typedef struct epswr_def_figure_t epswr_figure_t;
  /* A {epswr_figure_t} is an opaque object that can be used to write
    Postscript commands into an EPS file. It keeps track of some
    elements of the Postscript graphics state (native, or defined my
    this module) to reduce the file size and provide a hopefully more
    intuitive semantics.

    The graphics operations below, such as {epswr_segment}, take an
    {epswr_figure_t} as argument, and write the equivelent
    Postscript language commands to those file(s).
    
    COORDINATE SYSTEMS
    
    There are three main coordinate systems used by this interface. 
    The /Device/ coordinates are always in pt, measured from the lower
    left corner of the whole EPS figure.  These coordinates are usually
    named {h} (horizontal) and {v} (vertical).
    
    The /Client/ coordinates are in units chosen by the user and
    have arbitrary origin. They are usually named {x} (horizontal)
    and {y} (vertical). For technical reasons, the scaling factors
    {dx/dh} and {dy/dv} must always be equal.
    
    All Postscript commands written to the file use the Device
    coordinate system. The conversion from Client to Device coordinates
    is made by this module, before writing the data to the file.
    
    An {epswr_figure_t} also keeps track of two axis-aligned rectangles,
    the /Device window/ (in Device coordinates) and the /Client window/
    (in Client coordinates). The mapping from Client to Device
    coordinates is the unique biaffine map that takes the Client window
    to the Device window.  Most drawing commands are clipped 
    to this rectangle.
    
    The Device window is mirrored in the Postscript file by four Postscript
    variables {hMin,hMax,vMin,vMax}. */

typedef enum { epswr_axis_HOR = 0, epswr_axis_VER = 1 } epswr_axis_t;

/* STARTING A NEW EPS FIGURE */

epswr_figure_t *epswr_new_figure
  ( FILE *wr,
    double hPlotSize,  /* Initial plot window width (in pt). */
    double vPlotSize,  /* Initial plot window height (in pt). */
    double leftMargin,   /* Extra margin at left (pt). */
    double rightMargin,  /* Extra margin at right (pt). */
    double botMargin,    /* Extra margin at bottom (in pt). */
    double topMargin,    /* Extra margin at top (in pt). */
    bool_t verbose     /* TRUE to print diagnostics. */
  );
  /* Writes to file {wr} the suitable preamble for an EPS
    figure.
    
    The plot window initially will be a rectangle {hPlotSize} by
    {vPlotSize} points. The figure's bounding box will be that rectangle
    plus extra margins of {leftMargin, rightMargin, botMargin,
    topMargin} points on the four sides.
    
    The origin of the Device coordinate system will be at 
    the bottom left corner of the bounding box.  Therefore,
    the lower left corner of the plot window will be
    initially at Device coordinates {(leftMargin,botMargin)}.
    
    The Client coordinates are initially set to coincide with
    the Device coordinates (in pt), except that the origin will
    be the lower left corner of the plot window. */

void epswr_end_figure(epswr_figure_t *eps);
  /* Properly terminates the Encapsulated Postscript file, by writing
    the postamble, and closes the file. */

void epswr_get_figure_size
  ( epswr_figure_t *eps,  /* Picture stream. */
    double *hSizeP,        /* OUT: Total width of figure (in pt). */
    double *vSizeP         /* OUT: Total height of figure (in pt). */
  );
  /* Sets {*hSizeP} and {*vSizeP} to the total width and
    height (in pt) of the given EPS figure, including 
    its margins. */

/* PLOT WINDOW */

/* An {epswr_figure_t} object {eps} keeps track of a rectangle on the
  Device coordinate plane, called the /plot window/. This rectangle is
  used in some operations like {epswr_frame} and {epswr_grid_cell}, and
  is a cliiping path for most other operations.
  
  The plot window is represented internally by its coordinate bounds
  {eps.hMin,eps.hMax,eps.vMin,eps.vMax}, in absolute Device
  coordinates, which are mirrored in Postscript variables of the same
  names.  */

void epswr_set_device_window
  ( epswr_figure_t *eps,
    double hMin, double hMax,
    double vMin, double vMax,
    bool_t relative
  );
  /* Sets the plot window (and default clip path) to the rectangle {R =
    [hMin _ hMax] × [vMin _ vMax]}, in Device coordinates (pt).
    
    If {relative} is false, the coordinates are assumed to be relative
    to the lower left corner of the figure's bounding box. If {relative}
    is true, they are interpreted relative to the lower left corner of
    the current plot window.
    
    In any case, the procedure requires {hMin < hMax}, {vMin < vMax},
    and that the rectangle {R} be contained in the figure's bounding
    box.
    
    The Client to Device mapping is reset so that the Client units are
    points and the Client origin is the lower left corner of the plot
    window. The Client should use {epswr_set_client_window} to change
    the Client to Device coordinate mapping. The geometry and position
    of the text area (see {epswr_set_text_device_geometry} is not
    changed. */

void epswr_get_device_window
  ( epswr_figure_t *eps,
    double *hMinP, double *hMaxP,
    double *vMinP, double *vMaxP
  );
  /* Sets the variables {*hMinP,*hMaxP,*vMinP,*vMaxP} to the bounds of
    the current plot window in Device coordinates (pt, relative to the
    lower left corner of the EPS figure). */

void epswr_shrink_device_window
  ( epswr_figure_t *eps, 
    double dhMin, double dhMax, 
    double dvMin, double dvMax
  );
  /* Shrinks the Device plottng window by displacing the left edge by
    {dhMin}, the right edge by {dhMax}, the bottom edge by {dvMin}, and
    the top edge by {dvMax}, all inwards. Fails if the size of the plot
    window becomes zero.
    
    The Client to Device mapping is reset as in
    {epswr_set_device_window} The geometry and position of the text area
    (see {epswr_set_text_device_geometry} is not changed. */
    
void epswr_set_device_window_to_grid_cell
  ( epswr_figure_t *eps, 
    double hMin, double hMax, int ih, int nh, 
    double vMin, double vMax, int iv, int nv
  );
  /* Conceptually divides the rectangle {[hMin _ hMax] × [vMin _ vmax]}
    (in absolute Device coordinates) into a grid of rectangular
    cells with {nh} columns, numbered {0..nh-1} from left to right,
    and {nv} rows, numbered {0..nv-1} from bottom to top;
    and sets the plot window to the cell in column {ih} and row {iv},
    using {epswr_set_device_window}.
    
    The Client to Device mapping is reset as in
    {epswr_set_device_window}. The geometry and position of the text area
    (see {epswr_set_text_device_geometry} is not changed. */ 

/* CLIENT COORDINATE SYSTEM */

/* Most coordinate parameters and some size parameters accepted by
  functions in this interface are assumed to be in a /Client coordinate
  system/. Each Client coordinate (generally denoted {x} or {y}) is
  related to the corresponding Device coordinate ({h} or {v},
  respectively) by an affine trasnformation.
  
  The current mapping from Client to Device coordinates is represented
  internally by four other variables
  {eps.xMin,eps.xMax,eps.xMin,eps.xMax}, which are the Client
  coordinates that correspond to the above Device coordinate bounds.
  These four variables will be called here the /Client window/. */

void epswr_set_client_window
  ( epswr_figure_t *eps,
    double xMin, double xMax,
    double yMin, double yMax
  );
  /* Sets the Device-to-Client coordinate mapping based on the 
    current plot window.
    
    The mapping is such that Client abscissas {xMin} and {xMax}
    are mapped to {hMin} and {hMax}, and the Client ordinates
    {yMin} and {yMax} are mapped to {vMin} and {vmax}; where
    {[hMin _ hMax] x [vMin _vMax]} is the current plotting window
    in Device coordinates. 
    
    Currently, the absolute scale factors of the two mappings must be
    positive and equal, apart from floating point rounding errors. To
    ensure that constraint, the procedure may first reduces either the
    width or the height of the plot window, if needed, while preserving
    its center, so that the Device-to-Client scale factors {dh/dx} and
    {dv/dy} are equal.
    
    The Device coordinates and dimensions of the text area (see
    {epswr_set_text_geometry} are not changed.  Note that its
    Client coordinates and dimensions generally will. */

void epswr_get_client_window
  ( epswr_figure_t *eps,
    double *xMinP, double *xMaxP,
    double *yMinP, double *yMaxP
  );
  /* Sets the variables {*xMinP,*xMaxP,*yMinP,*yMaxP} to the bounds
    of the current plot window (in Client coordinates). */

void epswr_set_window
  ( epswr_figure_t *eps,
    double hMin, double hMax,
    double vMin, double vMax,
    bool_t relative,
    double xMin, double xMax,
    double yMin, double yMax
  );
  /* Equivalent to
    {epswr_set_device_window(eps,hMin,hMax,vMin,vMax,relative)}
    followed by
    {epswr_set_client_window(eps,xMin,xMax,yMin,yMax)}. */

/* CLIENT TO DEVICE COORDINATE CONVERSION */ 

void epswr_x_to_h_coord(epswr_figure_t *eps, double x, double *hP);
void epswr_y_to_v_coord(epswr_figure_t *eps, double y, double *vP);
  /* These procedures convert a Client coordinate ({x} or {y}) to the 
    corresponding absolute Device coordinate ({*hP} or {*vP}). */

void epswr_h_to_x_coord(epswr_figure_t *eps, double h, double *xP);
void epswr_v_to_y_coord(epswr_figure_t *eps, double v, double *yP);
  /* These procedures convert an absolute Device coordinate ({h} or {v}) to the 
    corresponding Client coordinate ({*xP} or {*yP}). */

void epswr_x_to_h_dist(epswr_figure_t *eps, double dx, double *dhP);
void epswr_y_to_v_dist(epswr_figure_t *eps, double dy, double *dvP);
  /* These procedures convert a Client distance along a coordinate axis
    ({dx} or {dy}) to the corresponding Device distance ({*dhP} or
    {*dvP}). */

void epswr_h_to_x_dist(epswr_figure_t *eps, double dh, double *dxP);
void epswr_v_to_y_dist(epswr_figure_t *eps, double dv, double *dyP);
  /* These procedures convert a Device distance along a coordinate axis
    ({dh} or {dv}) to the corresponding Client distance ({*dxP} or
    {*dyP}). */

void epswr_xy_to_hv_dist(epswr_figure_t *eps, double dxy, double *dhvP);
void epswr_hv_to_xy_dist(epswr_figure_t *eps, double dhv, double *dxyP);
  /* These procedures convert between Client and Device distances along
    arbitrary directions. The result will be correct for displacements
    in any direction only if the absolute scale factors on both axes are
    the same. Otherwise, the result will be intermediate between the
    correct values for displacement among the two axes. */

/* DRAWING COMMANDS */

void epswr_set_pen 
  ( epswr_figure_t *eps,
    double R, double G, double B,
    double width,
    double dashLength,
    double dashSpace
  );
  /* Sets pen parameters and ink color for line/outline drawing.
    Dimensions are in *millimeters*. Fails if any of {R,G,B} is
    {NAN}; otherwise clips them to the range [0_1]. */

void epswr_segment
  ( epswr_figure_t *eps,
    double xa, double ya,
    double xb, double yb
  );
  /* Draws segment from {(xa,ya)} to {(xb,yb)}. */
  
void epswr_curve
  ( epswr_figure_t *eps,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double xd, double yd
  );
  /* Draws a Bezier arc with given control points. */

void epswr_coord_line
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double pos
  );
  /* Draws a reference line PERPENDICULAR to the given axis 
    at the given coordinate value, extending across the whole
    plot window. */

void epswr_axis
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double pos, 
    double lo, 
    double hi
  );
  /* Draws an axis line PARALLEL to the given axis, spanning the
    cordinate interval {lo _ hi}, on that axis, and positioned at
    coordinate {pos} on the opposite axis. */

void epswr_frame (epswr_figure_t *eps);
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

void epswr_set_fill_color(epswr_figure_t *eps, double R, double G, double B);
  /* Defines the fill color for subsequent filling operations on
    {eps}. If {R} is negative,{±INF}, or {NAN}, the color is
    interpreted as "invisible". */

void epswr_rectangle
  ( epswr_figure_t *eps,
    double xlo, double xhi,
    double ylo, double yhi,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the given rectangle. */
  
void epswr_triangle
  ( epswr_figure_t *eps,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the triangle with corners {a,b,c}. */
  
void epswr_quadrilateral
  ( epswr_figure_t *eps,
    double x00, double y00,
    double x01, double y01,
    double x10, double y10,
    double x11, double y11,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the quadrilateral with corners
   {(x00,y00),(x01,y01),(x10,y10),(x11,y11)}. IMPORTANT: The corners
   must be given in row-by-row order, not CCW order. */
  
/* !!! Polygons should do the right thing with dashes. !!! */
   
void epswr_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double x[], double y[],
    int n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  );
  /* Fills and/or outlines the polygon with corners 
   {(x[0],y[0]),.. (x[n-1],y[n-1])}.
   
   If {closed} is true the polygon is implicitly closed
   with the segment from {(x[n-1],y[n-1])} to {(x[0],y[0])}.
   Otherwise the procedure draws an open polgonal line.
   
   Filling is applied only if {closed} and {fill} are true.
   A point is considered to be inside the polygon iff its winding
   number is odd (when {evenOdd} is TRUE) or nonzero (when {evenOdd}
   is FALSE). */
  
void epswr_rounded_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double x[], double y[],
    int n,
    double rad,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  );
  /* Same as {epswr_polygon}, but corners are rounded off with arcs
    of circles with radius {rad}.
    
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
    
void epswr_bezier_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double x[], double y[],
    int n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  );
  /* Fills and/or outlines a curvilinear polygon defined by {n}
    Bézier arcs. Arc number {i} has control points
    {(x[4*i+j],y[4*i+j])}, for {i} in {0..n-1} and {j} in {0..3}. If
    the end of each arc is not the beginning of the next arc, adds a
    connecting straight line. A point is considered to be inside the
    polygon iff its winding number is odd (when {evenOdd} is TRUE)
    or nonzero (when {evenOdd} is FALSE). */
  
void epswr_circle
  ( epswr_figure_t *eps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the circle centered at {(xc,yc)} with
    radius {rad}. */

void epswr_lune
  ( epswr_figure_t *eps,
    double xc, double yc, double rad, double tilt,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the lune with given center, radius, and
    tilt (in degrees CCW). (A "lune" is the intersection of two circles with the
    given radius, and with their centers spaced {rad} apart.) */
  
void epswr_slice
  ( epswr_figure_t *eps,
    double xc, double yc, double rad,
    double start, double stop,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the pie slice with given center, 
    radius, and angle range (in degrees). */

/* MARKER FILLING & DRAWING COMMANDS */

/* All commands in this section, the size of the figure is specified
  in millimeters, irrespecive of the current Client to-Device scale.
  The specified size does not include the line width that is used
  when {draw} is true. */
    
void epswr_dot
  ( epswr_figure_t *eps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  );
  /* Same as {epswr_circle}, except that the radius is in
    millimeters, not Client units. */
  
void epswr_tic
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double xc, double yc, 
    double ticSize,
    double align 
  );
  /* Draws a tic mark (short segment) at coordinates {(xc,yc)}. The
    segment will be perpendicular to the given {axis} and its length
    will be {ticSize} millimeters, irrespective of the current scale)
    The segment will extend {align*ticSize} mm in the negative
    direction, and {(1-align)*ticSize} mm in the positive
    direction. */
  
void epswr_cross
  ( epswr_figure_t *eps, 
    double xc, double yc, double rad, bool_t diag,
    bool_t draw
  );
  /* Draws a cross centered at {(xc,yc)} with radius {rad}. If
    {diag} is false the cross has vertical and horizontal branches;
    if {diag} is true the cross is rotated 45 degres. The radius is
    in millimeters, irrespective of the current scale. The mark has
    no interior. */
  
void epswr_asterisk
  ( epswr_figure_t *eps, 
    double xc, double yc, double rad,
    bool_t draw
  );
  /* Draws asn asterisk consisting of four crossed strokes at
    {(xc,yc)} with radius {rad}. The {radius} is in millimeters,
    irrespective of the current scale. The mark has no interior. */

void epswr_square
  ( epswr_figure_t *eps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or draws a square with center {xc,yc} and
    circum-radius {rad}. The radius is in millimeters, irrespective
    of the current scale. */

void epswr_diamond
  ( epswr_figure_t *eps, 
    double xc, double yc,
    double xRad, double yRad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or draws a diamond with center {xc,yc}, width {2*xRad}
    and height {2*yRad}. The last two parameters are in millimeters,
    irrespective of the current scale. */
    
void epswr_arrowhead 
  ( epswr_figure_t *eps,
    double xa, double ya, double xb, double yb,
    double width, double length, 
    double fraction,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines a triangular head for an arrow with base
    at {a = (xa,ya)} and tip at {b = (xb,yb)}. The head will have
    the specified {width} and {length} (in millimeters, irrespective
    of the current plot scale) and its tip will be positioned at the
    given {fraction} of the way from {a} to {b}. */
    
/* GRID LINES AND GRID CELLS */

void epswr_grid_lines(epswr_figure_t *eps, int nh, int nv);
  /* The current plot window is implicitly divided into a grid of rectangular
    cells, with {xn} columns and {yn} rows. Draws all grid cell boundaries
    with the current pen and color.  
    
    The first and last lines will extend beyond the current plot window
    by half of the current pen width. */

void epswr_grid_cell
  ( epswr_figure_t *eps, 
    int ih, int nh,
    int iv, int nv,
    bool_t fill, bool_t draw
  );
  /* The curren plot window is implicitly divided into a grid of
    rectangular cells, with {nh} columns and {nv} rows. Fills and/or
    outlines the cell in column {ih} and row {iv} of that grid. Cell
    {[0,0]} lies at the bottom left corner.
    
    The borders of the cells in the first and last column or row will
    extend beyond the current plot window by half of the current pen
    width. */

/* LABELS */

void epswr_set_label_font(epswr_figure_t *eps, const char *font, double size);
  /* Sets the name and point size of the font to be used by {epswr_label}. */

void epswr_label
  ( epswr_figure_t *eps, 
    const char *text, 
    double x, double y, 
    double rot,
    bool_t clipped,
    double hAlign, double vAlign,
    bool_t fill, bool_t draw
  );
  /* Prints {label} at point {(x,y)}, using the current label font
    and font size.
    
    Treats the outline of the characters as a path that is to be
    filled with the current fill color and/or stroked with the
    current pen.
    
    The parameter {hAlign} specifies which point of the string's
    bounding box will end up at {(x,y)}: 0.0 means the left side,
    1.0 means the right side. Other values of are
    interpolated/extrapolated, so that {hAlign = 0.5} the label will
    be horizontally centered at {(x,y)}.  The parameter 
    {vAlign} controls the vertical alignment of the bounding 
    box, with 0.0 meaning bottom and 1.0 meaning top.
    
    In any case, after the alignment is applied, the label will be
    rotated by {rot} degrees counterclockwise around the point
    {(x,y)}.
    
    If {clipped} is TRUE, the label will be clipped to the 
    current plot area.  Otherwise it may extend over the 
    whole figure.
    
    Beware that the label's bounding box is computed from the actual
    outlines of the character in the string, not from the nominal
    character boxes. In particular, leading and trailing blanks are
    ignored.
    
    !!! Add {hoff,voff} in millimeters !!! */

/* RUNNING TEXT */

void epswr_set_text_client_geometry
  ( epswr_figure_t *eps, 
    double xMin, double xMax, 
    double yMin, double yMax,
    double rot
  );
  /* Sets the nominal text margins for {epswr_text} to the rectangle {[xMin _ xMax] x
    [yMin x yMax]} in Client coordinates. 
    
    The {epswr_figure_t} object {eps} will keep track of the current
    /topline/ coordinate {yT}, initially at {yMax}. Each line of text
    written by {epswr_text} will be drawn with the top of its bounding
    box at {yT}, and then {yT} will be reduced by an amount equal to the
    current point size.  The coordinates {xMin} and {xMax} are used to
    define the alignment of each line; see {epswr_text}.
    
    The actual text will be rotated by {rot} degrees counterclockwise
    relative to the center of that rectangle. The {yMin} parameter is
    used only to define that center.
    
    !!! Should use the baseline of the text instead of the char's bounding box. !!! */

void epswr_set_text_font(epswr_figure_t *eps, const char *font, double size);
  /* Sets the name and point size of the font to be used by {epswr_text}. */

void epswr_text
  ( epswr_figure_t *eps, 
    const char *text, 
    bool_t clipped,
    double hAlign, 
    bool_t fill, bool_t draw
  );
  /* Writes one or more lines of string {text} after the lines already
    written. 
    
    Line breaks ("\n" or "\r\n" or just "\r") in the text are honored; a
    string that contains {k} line breaks generates {k+1} successive
    lines of text. Note that an empty {text} generates one blank line,
    and a {text} that ends with a line break generates a blank line
    after the preceding text. In any case, the {text} MUST be delimited
    by a zero byte ('\000').
    
    The position of the text is determined by the rectangle {[hMin_hMax]
    x [vMin x vMax]} defined in the last call to
    {epswr_set_text_geometry} and by the text lines already written
    since that call. Successive lines are written top to bottom in that
    rectangle, with vertical spacing determined by the current font's
    nominal size.
    
    More precisely, the bottom of the first line of {text} will be {dy}
    points below the bottom of the last line previously written, where
    {dy} is the current nominal text font size (as set by
    {epswr_set_text_font}); or at top of the text rectangle, if no lines
    have been written since the last call to
    {epswr_set_text_device_geometry} Successive baselines are spaced
    {dy} points down.
    
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
    
    If {clipped} is TRUE, the text will be clipped to the current plot
    area. Otherwise it may extend over the whole figure.
    
    !!! The vertical spacing {dy} should be an independent parameter? !!!
    */

/* MISCELLANEOUS */ 

void epswr_set_verbose(epswr_figure_t *eps, const bool_t verbose);
  /* Sets an internal diagnostic flag of {eps}, initially FALSE, to
    {verbose}. When that flag is set, some functions print diagnostic
    output to {stderr}. */

void epswr_comment(epswr_figure_t *eps, const char *title);
  /* Writes a Postscript comment line to the underlying file.
    It is intended for documentation and debugging
    purposes only; it has no visible effect on the printed 
    page/figure. */

void epswr_show_stack(epswr_figure_t *eps, int code);
  /* Issues a command that, when the EPS file is interpreted
    by {gs} (Ghostscript) or similar software, will
    print the contents of the Postscript stack to {stderr}.
    
    The {code} is an arbitrary integer that can be used 
    to identify the instance of the operation; it will appear
    at the top of the printed stack, but will  be deleted
    after the printout. */

void epswr_flush (epswr_figure_t *eps);
  /* Flushes the underlying file. */

/* UTILITIES */

#define epswr_MAX_SIZE (72000000.0)
  /* Max figure size in pt (1 million inches), just for paranoia. */

#define epswr_mm (72.0/25.4)
  /* One millimeter in Postscript points. */

void epswr_check_param
  ( const char *name, 
    double z, 
    double zMin, 
    double zMax
  );
  /* Checks whether {z} is in the range {[zMin __ zMax]}; 
    aborts with error if not. */

double epswr_round_to_nice(double x);
  /* Rounds the absolute value of {x} up to a nice value (namely
    0.20, 0.25, 0.50, or 1.00 times a power of 10). The sign is
    preserved. If {x} is already one of those nice values (or too
    close to one), it should be rounded up to it. As special cases,
    if {x} is 0, returns 0; if {|x|} is too large, returns
    {±INF}. */

double epswr_round_dim_mm(double x_mm, double dpi);
  /* Interprets {x_mm} as a dimension in millimeters and rounds it
    to the nearest non-zero even integer multiple of the printer's
    dot size, as specified by {dpi} (dots per inch). Howeverm, if
    {x_mm} is exactly zero, the result is zero. If {dpi} is zero,
    returns {x_mm} itself. */

#endif
