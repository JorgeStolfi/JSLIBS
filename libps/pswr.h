/* Tools for generating Postscript graphics files. */
/* Last edited on 2019-05-06 06:25:22 by jstolfi */

#ifndef pswr_H
#define pswr_H

/*
  This module provides tools to simplify the generation of 
  Postscript (PS) files containing graphics and text. 
  
  The output can be a single file containing a standalone Postscript
  document, or a sequence of files, each containing one Encapsulated
  Postscript (EPS) figure.

  STANDALONE POSTSCRIPT DOCUMENTS

  A standalone (non-EPS) file contains a single complete document,
  with any number of pages, meant to be printed on paper of some
  standard size, such as "letter" or "a4". A standalone PS file has a
  "%%Page" structuring comment before each page, explicit "showpage"
  commands at the end of each page, and (usually) no "%%BoundingBox"
  specification. It is usually meant to be printed by itself, and
  cannot be easily included as a figure of some other document.

  ENCAPSULATED POSTSCRIPT

  An EPS file contains a single figure, of arbitrary size. It has a
  "%%BoundingBox" line at the beginning, no "%%Page" lines, and no
  explicit "showpage" command. It is usually meant to be included in
  other documents. Some printers and viewers may not be prepared to
  handle an isolated EPS file.  In what follows, we assume that an
  EPS file contains a single ``page,'' whose dimensions are 
  given by the "%%BoundingBox".

  CANVAS

  To simplify the description, we will use the term /canvas/ to mean
  either a single page of a non-EPS file, or a single EPS figure.
 
  TYPICAL USAGE
  
  The typical usage for this interface is 
  
    { PSStream *ps = pswr_new_stream(...);               }
    { for each picture, do                               }
    {     pswr_new_picture(ps, xMin, xMax, yMin, yMax);  }
    {     pswr_segment(ps, ...)                          }
    {     ... other drawing commands ...                 }
    { pswr_close_stream(ps);                             }

  With this setup, each call to {pswr_new_picture} will start a new
  canvas: that is, a new page if {ps} is a non-EPS stream, and a new
  figure if {ps} is an EPS stream. In the second case, the call to
  {pswr_new_picture} will close the current EPS file, if any, and open
  a new one, with a sequentially numbered name.
  
  MULTIPLE PICTURES PER PAGE/FIGURE
  
  This module also provides facilities to draw multiple pictures on
  the same canvas. At any time, the client may call
  {pswr_set_canvas_layout} to specify that the canvas is to be divided
  into an array of pictures, with specified number of rows and
  columns; with extra blank space between pictures, including space
  below each picture for captions.

  In that case, each call to {pswr_new_picture} will advance to the
  the next available slot in the current canvas, in occidental reading
  order. Once all slots in the current canvas have been used, the next
  call to {pswr_new_picture} will start a new canvas, as above.
  
 */
  
#include <stdio.h>

#include <pswr_def.h>
#include <bool.h>

typedef enum { HOR = 0, VER = 1 } pswr_axis_t;

/* POSTCRIPT STREAM OBJECT */

typedef struct PSStreamPrivate PSStream;
  /* A {PSStream} is an opaque object that can be used to generate a
    Postscript (PS) document file, or a sequence of Encapsulated
    Postscript (EPS) figure files, containing graphics and/or text.
    The choice between the two is made at the stream's creation.

    The graphics operations below, such as {pswr_segment}, take a
    {PSStream} as argument, and write the necessary Postscript
    language commands to those file(s).

    After creating a {PSStream} of the desired type, the client need
    only specify the start of each new picture, and issue drawing
    commands into it. The {PSStream} routines will issue the
    approrpiate commands to either start a new page of a standalone PS
    file, or start a new separate EPS file, as needed. */

PSStream *pswr_new_stream
  ( char *prefix,       /* Output filename prefix.                      [Changed on 2009-01-05] */
    FILE *file,         /* Optional open file handle.                   [Changed on 2009-01-05] */
    bool_t eps,         /* TRUE for EPS figures, FALSE for PS document. [Changed on 2009-01-05] */
    char *docName,      /* Document name for PS output.                 [Changed on 2009-01-05] */
    char *paperSize,    /* Paper size ("letter", "a3", etc.). */
    bool_t landscape,   /* Paper orientation. */
    double hCanvasSize, /* Total canvas width (in pt). */
    double vCanvasSize  /* Total canvas height (in pt). */
  );
  /* Allocates a new {PSStream} record and initializes it appropriately.  
    
    If {eps} is false, the pictures will be written to the given {file},
    or, if {file} is NULL, to new file called "{prefix}{docName}.ps".
    Either way, the file will be a single standalone printable document,
    possibly with multiple pages numbered from 1. If the {papersize}
    string is neither {NULL} nor empty, the page dimensions are set
    from it and the {landscape} argument (see {pswr_get_paper_dimensions} below),
    and {hCanvasSize,vCanvasSize} are ignored.  Otherwise the page dimensions are set 
    to {hCanvasSize,vCanvasSize}. In any case, the {landscape} flag is used to set the
    document's printing orientation.
    
    If {eps} is TRUE, the EPS figures will be written out as a
    sequence of files called "{prefix}{page}.eps", where {page} is a
    six-digit sequential figure number or a client-given page name
    (see {pswr_new_canvas}). However, if the {file} argument is not
    NULL, then the first EPS figure will be written to that file. The
    {docName}, {paperSize}, and {landscape} arguments are ignored,
    and the EPS bounding box size is set to {hCanvasSize} by {vCanvasSize}.
    
    By default, a new stream will have one picture per canvas, using
    the whole canvas area except for margins and separators. In
    standalone (non-EPS) files, the margins are 1 inch wide all around
    the page, plus extra space for 5 lines of caption (50pt) under the
    picture. In EPS files, the margin is 4pt all around the figure,
    and there is no extra space for caption.
    
    NOTE: The meaning of {prefix} was slightly changed on 2009-01-05.
    For ".ps" output, the name now is "{prefix}{docName}.ps" instead
    of just "{prefix}.ps". For ".eps" output, the figure names are now
    "{prefix}{page}.eps" instead of {prefix}-{page}.eps". The order of
    the first parameters was changed to make the old uses stand
    out. */

void pswr_close_stream(PSStream *ps);
  /* Terminates any pictures that have been written to {ps}, 
    flushes and closes any open files. Then frees all internal
    storage used by {ps}. */
  
bool_t pswr_is_eps(PSStream *ps);
  /* Returns TRUE iff {ps} is an Encapsulated Postscript stream. */

/* STARTING A NEW PICTURE */

void pswr_new_picture
  ( PSStream *ps,              /* Postscript picture stream. */
    double xMin, double xMax,  /* Client X plotting range. */
    double yMin, double yMax   /* Client Y plotting range. */
  );
  /* Signals the start of a new picture. The plot window is set (with
    {pswr_set_window}) to the next available picture slot in the
    current canvas. A new canvas is started if necessary.
     
    Regardless of the {eps} flag, the plotting scale parameters are
    set up so that the rectangle {[xMin _ xMax] × [yMin _ yMax]} in
    client coordinates will fit centered in the selected picture slot
    (excluding its margins), with equal scale factors on both axes. */

/* CANVAS SIZE */

void pswr_set_canvas_size
  ( PSStream *ps,       /* Picture stream. */
    double hCanvasSize,   /* Width of canvas (in pt). */
    double vCanvasSize    /* Height of canvas (in pt). */
  );
  /* If {ps} is an EPS stream, changes the canvas size (i.e. the
    figure size) to the given dimensions in pt. The current EPS file,
    if any, is closed. If {ps} is a standalone document, the operation
    has no effect (the canvas size is determined by the document paper 
    size given at stream creation time). */

void pswr_get_canvas_size
  ( PSStream *ps,       /* Picture stream. */
    double *hCanvasSize,  /* OUT: Width of canvas (in pt). */
    double *vCanvasSize   /* OUT: Height of canvas (in pt). */
  )  ;
  /* Sets {*hCanvasSize} and {*vCanvasSize} to the current width and
    height (in pt) of the document's pages (if {ps} a standalone
    document) or current EPS figure (if {ps} an EPS stream). */

/* PLOTTING WINDOW */

void pswr_set_canvas_window
  ( PSStream *ps,
    double hMin, double hMax,
    double vMin, double vMax
  );
  /* Sets the nominal plotting area in the current canvas to an
    arbitrary rectangle in the current canvas, ignoring the picture
    layout.

    After this call, the nominal plotting area will be {[hMin _ hMax]
    × [vMin _ vMax]}, in pt, relative to the lower left corner of the
    figure (EPS) or the current page (non-EPS).
    
    The client's plotting rectangle is reset to the same rectangle.
    That is, after this call the client's coordinates will be
    interpreted as distances in {pt} from the lower left corner
    of the canvas.  The client should use {pswr_set_client_window}
    to change the scale.
    
    The caption cursor will be positioned to the first text line under
    the window. Except for {pswr_add_caption} and {pswr_frame}, all
    graphics commands will be clipped to the plotting window.
    
    The cell grid (see {pswr_set_grid} below) is reset to a single
    cell covering the whole window.
    
    The procedure also marks all slots in the current canvas full, so
    that the next call to {pswr_new_picture} will start a new
    canvas. */

void pswr_set_client_window
  ( PSStream *ps,
    double xMin, double xMax,
    double yMin, double yMax
  );
  /* Sets the client's plotting window to the given rectangle in
    client coordinates. After this call, client coordinates will range
    over {[xMin _ xMax] × [yMin _ yMax]}. The actual ploting are in
    the canvas is not afected. IMPORTANT: The aspect ratio of the
    client window must be equal to that of the canvas, so that the
    plotting scales {dh/dx} and {dv/dy} for both axes (in pt per
    client unit) must be the same. */

void pswr_set_window
  ( PSStream *ps,
    double xMin, double xMax,
    double yMin, double yMax,

    double hMin, double hMax,
    double vMin, double vMax
  );
  /* Equivalent to {pswr_set_canvas_window(ps,hMin,hMax,vMin,vMax)} 
    followed by {pswr_set_client_window(ps,xMin,xMax,yMin,yMax)}. */

void pswr_get_paper_dimensions(const char *papersize, bool_t landscape, double *xpt, double *ypt);
  /* Sets *xpt and *ypt to the dimensions of the specified paper type,
    in points. Knows about US sizes "letter", "ledger", "tabloid",
    "legal", "executive", and the ISO "A" sizes (from "4A0" to
    "A10"). 
    
    If {landscape} is false, assumes that the paper is in "portrait"
    orientation, witht {*xpt <= *ypt}.  If {landscape} is true, assumes
    the "landscape" orientation, thus swapping {*xpt} and {*ypt} . */

/* MULTIPLE PICTURES PER PAGE/FIGURE  */

void pswr_set_canvas_layout
  ( PSStream *ps,         /* Picture stream. */
    double hPicSize,      /* Width of each picture (in pt). */
    double vPicSize,      /* Height of each picture (in pt). */
    bool_t adjustPicSize, /* TRUE to fit {hPicSize,vPicSize} to canvas size. */
    double hPicMargin,    /* Left/right margin for each picture (in pt). */
    double vPicMargin,    /* Top/bottom margin for each picture (in pt). */
    int captionLines,     /* Number of caption lines below each picture. */
    int hPicCount,        /* Number of pictures in each row. */   
    int vPicCount         /* Number of pictures in each column. */
  );
  /* Changes the number and layout of pictures per canvas (non-EPS
    page or EPS figure). After this call, each canvas will contain
    {vPicCount} rows of pictures, with {hPicCount} pictures each.
    
    Each picture will occupy a rectangle with proportions {hPicSize ×
    vPicSize}. It will be surrounded by a blank margin with
    dimensions {hPicMargin,vPicMargin}, and below the picture there
    will be additional space for {captionLines} lines of caption
    (at 10pt for line).
    
    If {adjustPicSize} is FALSE, then {hPicSize} and {vPicSize} are
    assumed to be actual dimensions in pt. If {adjustPicSize} is true,
    these parameters will be automatically scaled, by the same factor,
    so that the whole picture array fits as snugly as possible in the
    canvas, with a margin of 1 inch for documents,
    or 4pt for encapsulated figures. Similarly, any picture count
    ({hPicCount} and/or {vPicCount}) which is zero is computed automatically
    according to this same criterion.
    
    If any pictures have already been started in the current canvas,
    the latter is flushed, and a new canvas is started. */

void pswr_new_canvas(PSStream *ps, const char *pageName);
  /* Starts a new canvas, and resets the plot window to the whole
    canvas.
    
    For non-EPS streams, the {pageName} string is used as the first
    field of the "%%Page" structuring comment.
    
    For EPS streams, the function closes the current EPS file (if any)
    and opens a new one. the new file is called
    "{prefix}{pageName}.eps", where {prefix} is the argument given to
    {pswr_new_stream}. If {pageName} is NULL or empty, it defaults to
    a six-digit sequential canvas number (starting from "000000"). */

void pswr_sync_canvas(PSStream *ps, const char *pageName);
  /* Equivalent to {pswr_new_canvas}, if any {pswr_new_picture}
    or {pswr_set_window} has been issued to the current canvas.
    Otherwise it is a no-op. */

void pswr_fill_row(PSStream *ps);
  /* Marks all picture slots in the current row as full, so that the
    next call to {pswr_new_picture} will start in the next row
    (possibly in the next canvas). */

void pswr_fill_canvas(PSStream *ps);
  /* Marks all picture slots in the current canvas as full, so that the
    next call to {pswr_new_picture} will start in a new canvas. */

/* DRAWING COMMANDS */

void pswr_set_pen 
  ( PSStream *ps,
    double R, double G, double B,
    double width,
    double dashLength,
    double dashSpace
  );
  /* Sets pen parameters and ink color for line/outline drawing.
    Dimensions are in *millimeters*. Fails if any of {R,G,B} is {NAN};
    otherwise clips them to the range [0_1]. */

void pswr_segment
  ( PSStream *ps,
    double xa, double ya,
    double xb, double yb
  );
  /* Draws segment from {(xa,ya)} to {(xb,yb)}. */
  
void pswr_curve
  ( PSStream *ps,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double xd, double yd
  );
  /* Draws a Bezier arc with given control points. */

void pswr_coord_line (PSStream *ps, pswr_axis_t axis, double pos);
  /* Draws a reference line PERPENDICULAR to the given axis 
    at the given coordinate value, extending across the whole
    plot window. */

void pswr_axis(PSStream *ps, pswr_axis_t axis, double pos, double lo, double hi);
  /* Draws an axis line PARALLEL to the given axis, spanning the
    cordinate interval {lo _ hi}, on that axis, and positioned at
    coordinate {pos} on the opposite axis. */

void pswr_frame (PSStream *ps);
  /* Draws the outline of the current plotting area. The outline 
    will be drawn half-inside, half-outside the area. */

/* FIGURE FILLING & DRAWING COMMANDS */

/* For all commands of this section, if the {fill} argument
  is true, the specified figure is painted with the 
  current fill color; then, if {draw} is true, the outline
  of the figure is drawn with the current pen parameters. 
  
  However, if the {R} component of the current fill color is negative,
  {±INF}, or {NAN}, the color is considered "invisible", and the
  shapes will not be filled, irrespective of the {fill} argument.
  
  The fill color is initially 50% gray and is reset to that 
  value at every new canvas. */

void pswr_set_fill_color(PSStream *ps, double R, double G, double B);
  /* Defines the fill color for subsequent filling operations on {ps}.
    If {R} is negative,{±INF}, or {NAN}, the color is interpreted as
    "invisible". */

void pswr_rectangle
  ( PSStream *ps,
    double xlo, double xhi,
    double ylo, double yhi,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the given rectangle. */
  
void pswr_triangle
  ( PSStream *ps,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the triangle with corners {a,b,c}. */
  
void pswr_quadrilateral
  ( PSStream *ps,
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
   
void pswr_polygon
  ( PSStream *ps,
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
  
void pswr_rounded_polygon
  ( PSStream *ps,
    bool_t closed,
    double x[], double y[],
    int n,
    double rad,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  );
  /* Same as {pswr_polygon}, but corners are rounded off with arcs of circles 
    with radius {rad}.  
    
    The rounding of a corner will never remove an entire side, although
    the rounding of two adjacent corners may do so (and thus generate a
    straight joining line with negative length). If the rounding of a
    corner would require removing more than one entire side, the
    rounding of that corner is suppressed. In particular, this happens
    whenever one of the adjacent sides has zero length (i.e. two
    consecutive corners coincide) of the internal angle at the corner is
    0 or 360 degrees.
    
    If {closed} is true the procedure assumes that the midpoint of the
    segment from point {n-1} to point 0 is not removed by
    rounding, and will start drawing from there. If {closed} is false
    the drawing starts at point 0 and ends at point {n-1} ignoring
    rounding. */
    
void pswr_bezier_polygon
  ( PSStream *ps,
    bool_t closed,
    double x[], double y[],
    int n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  );
  /* Fills and/or outlines a curvilinear polygon defined by {n} Bézier
   arcs. Arc number {i} has control points {(x[4*i+j],y[4*i+j])}, for
   {i} in {0..n-1} and {j} in {0..3}. If the end of each arc is not
   the beginning of the next arc, adds a connecting straight line. A
   point is considered to be inside the polygon iff its winding number
   is odd (when {evenOdd} is TRUE) or nonzero (when {evenOdd} is
   FALSE). */
  
void pswr_circle
  ( PSStream *ps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the circle centered at {(xc,yc)} with radius {rad}. */

void pswr_lune
  ( PSStream *ps,
    double xc, double yc, double rad, double tilt,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the lune with given center, radius, and tilt.
    (A "lune" is the intersection of two circles with the 
    given radius, and with their centers spaced {rad} apart.) */
  
void pswr_slice
  ( PSStream *ps,
    double xc, double yc, double rad,
    double start, double stop,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines the pie slice with given center, 
    radius, and angle range (in degrees). */

/* MARKER FILLING & DRAWING COMMANDS */

/* All commands in this section, the size of the figure is specified in millimeters,
  irrespecive of the current plotting scale.  The specified size does not
  include the line width that is used when {draw} is true. */
    
void pswr_dot
  ( PSStream *ps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  );
  /* Same as {pswr_circle}, except that the radius is in millimeters,
    irrespective of the current scale. */
  
void pswr_tic
  ( PSStream *ps, 
    pswr_axis_t axis, 
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
  
void pswr_cross
  ( PSStream *ps, 
    double xc, double yc, double rad, bool_t diag,
    bool_t draw
  );
  /* Draws a cross centered at {(xc,yc)} with radius {rad}. If {diag}
    is false the cross has vertical and horizontal branches; if {diag}
    is true the cross is rotated 45 degres. The radius is in
    millimeters, irrespective of the current scale. The mark has
    no interior. */
  
void pswr_asterisk
  ( PSStream *ps, 
    double xc, double yc, double rad,
    bool_t draw
  );
  /* Draws asn asterisk consisting of four crossed strokes at
    {(xc,yc)} with radius {rad}. The {radius} is in millimeters,
    irrespective of the current scale. The mark has no interior. */

void pswr_square
  ( PSStream *ps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or draws a square with center {xc,yc} and circum-radius
    {rad}. The radius is in millimeters, irrespective of the current
    scale. */

void pswr_diamond
  ( PSStream *ps, 
    double xc, double yc,
    double xRad, double yRad,
    bool_t fill, bool_t draw
  );
  /* Fills and/or draws a diamond with center {xc,yc}, width {2*xRad}
    and height {2*yRad}. The last two parameters are in millimeters,
    irrespective of the current scale. */
    
void pswr_arrowhead 
  ( PSStream *ps,
    double xa, double ya, double xb, double yb,
    double width, double length, 
    double fraction,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines a triangular head for an arrow with 
    base at {a = (xa,ya)} and tip at {b = (xb,yb)}. The head will have 
    the specified {width} and {length} (in millimeters, 
    irrespective of the current plot scale) and its tip will be positioned
    at the given {fraction} of the way from {a} to {b}. */
    
/* GRID LINES AND GRID CELLS */

void pswr_set_grid(PSStream *ps, int xn, int yn);
  /* The plot window is implicitly divided into a grid of rectangular
    cells, with {xn} columns and {yn} rows. These cells are used by
    {pswr_grid_lines} and {pswr_grid_cell} below. */

void pswr_grid_lines(PSStream *ps);
  /* Draws all grid cell boundaries with the current pen and color. */

void pswr_grid_cell
  ( PSStream *ps, 
    int xi, int yi,
    bool_t fill, bool_t draw
  );
  /* Fills and/or outlines cell {[xi,yy]} of the curent grid.
     Cell {[0,0]} lies at the bottom left corner. */
    
/* TEXT PRINTING */

void pswr_set_label_font(PSStream *ps, const char *font, double size);
  /* Sets the name and point size of the font to be used by pswr_label. */

void pswr_label
  ( PSStream *ps, 
    const char *text, 
    double x, double y, 
    double rot, double xAlign, double yAlign
  );
  /* Prints {label} at point {(x,y)}, using the current label font size.
    The label will be rotated by {rot} degrees counterclockwise around
    the point {(x,y)}.
    
    The parameter {xAlign} (resp. {yAlign}) specifies which point of the
    string's bounding box will end up at {(x,y)}: 0.0 means the left
    (resp. bottom) side, 1.0 means the right (resp. top) side.
    
    Beware that the label's bounding box is computed from the
    character outlines, not from the nominal character boxes. In
    particular, leading and trailing blanks are ignored. */

void pswr_fill_draw_label
  ( PSStream *ps, 
    const char *text,
    double x, double y, 
    double rot, double xAlign, double yAlign,
    bool_t fill, bool_t draw
  );
  /* Same as {pswr_label}, but treats the outline of the characters
    as a path that is to be filled with the current fill color
    and/or stroked with the current pen. */

void pswr_add_caption(PSStream *ps, const char *txt, double xAlign);
  /* If {ps.captionLines} is zero, does nothing. If {ps.captionLines}
    is positive, adds a line of caption text below the current plot
    window, *outside* the nominal bounding box. Successive calls
    append succesive lines to the caption, even if the stated limit
    {ps.captionLines} is exceeded. Newlines in {txt} are honored. Each
    line is aligned as specified by {xAlign}: 0.0 = left aligned, 0.5
    = centered, 1.0 = right aligned. */

/* MISCELLANEOUS */ 

void pswr_set_verbose(PSStream *ps, const bool_t verbose);
  /* Sets an internal diagnostic flag of {ps},initially FALSE, to
    {verbose}. When that flag is set, some functions print diagnostic
    output to {stderr}. */

void pswr_comment (PSStream *ps, const char *title);
  /* Writes a Postscript comment line to the underlying file.
    It is intended for documentation and debugging
    purposes only; it has no visible effect on the printed 
    page/figure. */

void pswr_flush (PSStream *ps);
  /* Flushes the underlying file. */

double pswr_round_to_nice(double x);
  /* Rounds the absolute value of {x} up to a nice value (namely 0.20, 0.25,
    0.50, or 1.00 times a power of 10). The sign is preserved. If {x}
    is already one of those nice values (or too close to one), it
    should be rounded up to it. As special cases, if {x} is 0, returns
    0; if {|x|} is too large, returns {±INF}. */

double pswr_round_dim_mm(double x_mm, double dpi);
  /* Interprets {x_mm} as a dimension in millimeters and rounds it to the nearest non-zero 
    even integer multiple of the printer's dot size, as specified by {dpi} (dots per inch). 
    Howeverm, if {x_mm} is exactly zero, the result is zero.  If {dpi} is zero, 
    returns {x_mm} itself.  */

#endif
