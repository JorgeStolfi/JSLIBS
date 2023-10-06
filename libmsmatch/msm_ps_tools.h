#ifndef msm_ps_tools_H
#define msm_ps_tools_H

/* Postscript plots of graphs and such. */
/* Last edited on 2023-10-01 19:43:49 by stolfi */

#define msm_ps_tools_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <stdint.h>

#include <vec.h>
#include <epswr.h>

#define msm_EPS_MARGIN_MM (1.0)
  /* EPS figure margin width (in mm). */

/* POSTCRIPT FIGURE STREAMS */

typedef struct msm_ps_tools_t msm_ps_tools_t;
  /* A handle to a Postscript graph-plotting stream.
  
    An {msm_ps_tools_t} is a handle to an {epswr_figure_t} object plus some
    internal state. The {epswr_figure_t} object in turn includes a
    handle to a file that is supposed to contain one Encapsulated
    Postscript image. Plot commands on this interface (or in the {epswr.h} 
    interface) usually write Postscript commands to that file.
    
    Procedures that create an {msm_ps_tools_t} take as parameters a
    {FILE*} {wr} and two strings {name} and {tag}. If {wr} is not
    NULL, the figure will be written to that file, and the {name} and
    {tag} are ignored.
    
    If {wr} is NULL, the procedure will open a file called
    "{name}{tag}.eps" and use that file instead of {wr}. In both
    cases, the procedure will write the appropriate EPS preamble to
    the file.

    Many plotting routines in this interface take their arguments in
    /Graph coordinates/ {x,y}. A {msm_ps_tools_t} object defines a
    mapping of those Graph coordinates to the {epswr} ``client''
    coordinates (here called /Epswr coordinates/) {h,v}, measured in
    millimeters from the lower left corner of the plot area.
    
    The client must call {msm_ps_tools_close} at the end, to write the necessary 
    postamble and close the underlying {FILE}. 
  */
 
msm_ps_tools_t *msm_ps_tools_new
  ( FILE *wr,
    char *name,
    char *tag,
    double hSize, 
    double vSize,
    double fontSize,
    int32_t maxXLabChars,
    int32_t maxYLabChars,
    double mrg
  );
  /* Creates a {msm_ps_tools_t} for an Encapsulated Postscript
    figure of width {hSize} and height {vSize}, plus an extra margin
    of width {mrg} all around. All these dimensions are in mm.
    
    The parameters {maxXLabChars,maxYLabChars} should be an upper bound
    on the number of characters in the X and Y plot scale labels.
    
    The {epswr} client and device coordinate systems and windows,
    as well as the Graph coordinates of this intervace, will be set
    initially to the rectangle {[0_hSize]×[0_vSize]] excluding the
    margin.
    
    The {fontSize} (in pt) will be used for labels and titles of 
    graphs written to this {msm_ps_tools_t}. */

msm_ps_tools_t *msm_ps_tools_new_graph
  ( FILE *wr,
    char *name,
    char *tag,
    double hGraphSize, 
    double vGraphSize,
    bool_t scaleL, bool_t titleL,
    bool_t scaleR, bool_t titleR,
    bool_t scaleB, bool_t titleB,
    bool_t scaleT, bool_t titleT,
    double fontSize,
    int32_t maxXLabChars,
    int32_t maxYLabChars,
    double mrg
  );
  /* Creates a {msm_ps_tools_t} for an Encapsulated Postscript
    figure that is to contain a graph of width {hGraphSize} and height
    {vGraphSize} (in mm).
    
    The plottable area will include space for the graph, plus extra
    space around it for graph scales at the left, right, bottom or top
    (as specified by the flags {scaleL,scaleR,scaleB,scaleT} and/or
    titles (as specified by the flags {titleL,titleR,titleB,titleT}. The
    width of these extra spaces is estimated from the {fontSize}
    parameter (the label font's nominal size, in pt) and
    {maxXLabChars,maxYLabChars} (upper bounds to the number of
    characters in any label of the X scale).
    
    In addition to these extra spaces, the figure will have an extra
    margin of {mrg} mm all around the plottable area.
    
    The Graph and Epswr reference rectangles will be initialized
    to the graph area only (excluding the scale and title areas). */
  
void msm_ps_tools_close(msm_ps_tools_t *mps);
  /* Terminates any pictures that have been written to {mps}, flushes
    and closes the underlying file (unless it is {stdout} or
    {stderr}). Then frees all internal storage used by {mps}, including
    the associated {PSStream} object and {*mps} itself. */

void msm_ps_tools_get_plot_size(msm_ps_tools_t *mps, double *hSize, double *vSize);
  /* Stores in {*hSize} and {*vSize} the dimensions (in mm) of the usableplot
    area, as specified to {msm_ps_tools_init}, excluding the margin. */
   
/* COORDINATE MAPPING */
  
double msm_ps_tools_map_x(msm_ps_tools_t *mps, double x);
double msm_ps_tools_map_y(msm_ps_tools_t *mps, double y);
double msm_ps_tools_map_coord(msm_ps_tools_t *mps, epswr_axis_t axis, double coord);
void msm_ps_tools_map_coords(msm_ps_tools_t *mps, double x, double y, double *h, double *v);
  /* These procedures map Graph coordinates {x,y} to Epswr
     coordinates {h,v}, as defined by {mps}. */

double msm_ps_tools_unmap_h(msm_ps_tools_t *mps, double h);
double msm_ps_tools_unmap_v(msm_ps_tools_t *mps, double v);
double msm_ps_tools_unmap_coord(msm_ps_tools_t *mps, epswr_axis_t axis, double coord);
void msm_ps_tools_unmap_coords(msm_ps_tools_t *mps, double h, double v, double *x, double *y);
  /* These procedures map Epswr coordinates {h,v} to Graph
     coordinates {x,y}, as defined by {mps}. */
  
/* SETTING THE COORDINATE MAPPINGS

  The Graph-to-Epswr mapping of a {msm_ps_tools_t} object is defined by
  two rectangles stored in the object: the /Epswr reference window/
  {ERW = [hMin_hMax]×[vMin_vMax]} and the /Graph reference window/
  {GRW = [xMin_xMax]×[yMin_yMax]}. The {ERW} is distinct from 
  the {epswr} ``client'' and ``device'' windows.
  
  The Graph-to-Epswr mapping is such that the {GRW} is mapped to the
  {ERW}. More precisely, {x=xMin} is mapped to {h=hMin}, {x=xMax} to
  {h=hMax}. So, if {hMin > hMax} (or {xMin > xMax}) the Graph {x}
  axis is directed from right to left. The same holds for the {y} and
  {v} coordinates.
  
  Changing either of these rectangles affects the way that Graph
  coordinates are interpreted by subsequent function calls.

  The rectangles {GRW} and {ERW} are not used for clipping.
  Regardless of them, the usable plotting area is always
  the rectangle {[0_hSize]×[0_vSize]} in Epswr coordinates. */

void msm_ps_tools_set_epswr_ref_window
  ( msm_ps_tools_t *mps,
    double hMin, double hMax, 
    double vMin, double vMax
  );
  /* Sets the Epswr reference window to the specified
    rectangle (in mm from the lower left corner of the 
    plottable area).*/

void msm_ps_tools_set_graph_ref_window
  ( msm_ps_tools_t *mps,
    double xMin, double xMax, 
    double yMin, double yMax
  );
  /* Sets the Graph reference window to the specified
    rectangle (in Graph coordinates).*/

void msm_ps_tools_shrink_epswr_ref_window
  ( msm_ps_tools_t *mps,
    double lMrg, double rMrg, 
    double bMrg, double tMrg
  );
  /* Displaces the Epswr reference window boundary INWARDS by the
    specified amounts (in mm) on the left, right, bottom, and top sides,
    respectively. */

void msm_ps_tools_expand_graph_ref_window
  ( msm_ps_tools_t *mps,
    double lMrg, double rMrg, 
    double bMrg, double tMrg
  );
  /* Displaces the Graph reference window boundary OUTWARDS by the
    specified amounts on the left, right, bottom, and top sides,
    respectively. */

void msm_ps_tools_compute_data_range
  ( int32_t n, 
    int32_t stride, 
    double z[], 
    double *zMinP, 
    double *zMaxP
  );
  /* Computes a data range {[*zMinP _ *zMaxP]} appropriate for plotting the 
    graph of {n} given values, namely {z[i*stride]} for {i = 0..n-1]}.
    The resulting range will be non-empty and non-degenerate even if {n == 0}
    or all {z} values are equal. */

/* LOW-LEVEL DRAWING TOOLS 

  Unless specified otherwise, all procedures below take Graph
  coordinates. */

void msm_ps_tools_draw_segment(msm_ps_tools_t *mps, double xa, double ya, double xb, double yb);
  /* Draws a line segment from {(xa,ya)} to {(xb,yb)} with the current pen
    and draw color. */

void msm_ps_tools_draw_ref_axis(msm_ps_tools_t *mps, epswr_axis_t axis, double R, double G, double B);
  /* Draw the {x} or {y} coordinate axes in dark gray.
    The axis will span the current reference rectangle 
    (NOT the current plottable area). The axis will be drawn 
    with color {R,G,B}. */

void msm_ps_tools_draw_tic
  ( msm_ps_tools_t *mps, 
    epswr_axis_t axis, 
    double xt, 
    double yt,
    double ticSize,
    double ticAlign,
    double R, double G, double B,
    char *label,
    double labAlign
  );
  /* Draws a tic mark, perpendicular to the given {axis}, at the
    nominal position {xt,yt}, with the current pen and draw color.
    The parameters {ticSize,ticAlign} are
    passed to {epswr_tic}.
    
    If the string {label} is not NULL, its is written using the
    current font and text color. The label is centered on the tic
    position along the given {axis}, and aligned as specified by
    {labAlign}, relative to the tic, along the other axis. Thus, for
    example, {labAlign=1.2} would be appropriate for a tic along the
    left or bottom edges of the plot area, and {labAlign=-0.2} would
    do for tics along the top and right edges. The tic will be drawn 
    with color {R,G,B}, the text in black. */

void msm_ps_tools_draw_scale
  ( msm_ps_tools_t *mps, 
    epswr_axis_t axis, 
    double pos,
    double ticSize,
    double ticAlign,
    double ticMinDist, 
    double ticMinStep, 
    double R, double G, double B,
    char *fmt,
    double labAlign,
    double labMinDist,
    double labMinStep
  );
  /* Draws the tic marks for the specified {axis}, in black. 
    
    The ideal axis of the tics (the /ticline/) will be at coordinate
    {pos} on the other axis. The set of tic marks will span the
    reference rectangle. The spacing of the tics is chosen by the
    procedure, with the constraint that the spacing between tic marks
    (labeled or unlabeled) must be at least {ticMinDist} in millimeters, and
    {ticMinStep} in Graph coordinates. (At least one of these
    constraints must be positive.)
    
    If the parameter {fmt} is not NULL, a subset of the tics will be
    labeled with the corresponding Graph coordinates, in the format
    {fmt} (as per {printf}). Labeled tics are longer than unlabeled
    ones ({1.5*ticSize}, in mm). The spacing between labeled tics
    will be at least {labMinDist} in millimeters, and {labMinStep}
    in Graph coordinates.  
    
    The parameters {ticSize,ticAlign,labAlign} are explained
    under {msm_ps_tools_draw_tic}. The lines and tics will be drawn 
    with color {R,G,B}. */

void msm_choose_label_coords
  ( msm_ps_tools_t *mps,
    epswr_axis_t axis,
    double ztMin,
    double ztMax, 
    double ztStep,
    double minDist,
    double minStep,
    int32_t *labPerP,
    int32_t *labSkpP
  );;
  /* Chooses Graph coordinates for major (labeled) tics,
    given the spacing {ztStep} (in Graph coordinates) between
    minor tics, and the Graph coordinates {ztMin,ztMax} of the first and last
    minor tics.

    Specifically, the procedure chooses the number {labPer} (a
    positive integer) of ordinary tic steps for each major tic step,
    so that the label increment {labPer*ztStep} is a nice round value.
    It also computes the number {lapSkp} of minor tic steps between
    {ztMin} and the first major tic after it.
    
    If {minDist} is positive, it specifies the minimum spacing (in mm)
    between labeled tics. If {minStep} is positive, it specifies the
    minimum increment (in Graph coordinates) between consecutive
    labels.
    
    The results are returned in {*labPerP,*lasbSkpP}.  If the constraints
    imply that there are no major tics in the range {ztMin,ztMax},
    returns {*labPerP=0}. */ 
    
void msm_ps_tools_choose_tic_coords
  ( msm_ps_tools_t *mps, 
    epswr_axis_t axis, 
    double cMin, 
    double cMax, 
    double minDist, 
    double minStep,
    double *zMinP,
    double *zMaxP,
    double *zStepP
  );
  /* Chooses Graph coordinates for minor tics within the range
    {[cMin_cMax]} of Epswr coordinates.
    
    Specifically, the procedure chooses a nice round increment
    {zStep} in Graph coordinates between two consecutive tics, and
    computes the Graph coordinates {zMin,zMax}, multiples of
    {zStep}, of the first and last tic that fit in that range.
    
    The parameter {minDist}, if positive, and specifies the minimum
    distance (in mm) between tics. The parameter {minStep}, if
    positive, specifies the minimum Graph coordinate increment
    between tics. At least one of {minDist,minStep} must be
    positive.
    
    The results are returned in {*zMinP,*zMaxP,*zStepP}.
    If there are no nice coordinate values in the specified range,
    returns {*zStepP=0}, {*zMinP>*zMaxP}. */ 

void msm_psplot_choose_label_coords
  ( msm_ps_tools_t *mps,
    epswr_axis_t axis,
    double ztMin,
    double ztMax, 
    double ztStep,
    double minDist,
    double minStep,
    int32_t *labPerP,
    int32_t *labSkpP
  );
  /* Chooses Graph coordinates for major (labeled) tics,
    given the spacing {ztStep} (in Graph coordinates) between
    minor tics, and the Graph coordinates {ztMin,ztMax} of the first and last
    minor tics.

    Specifically, the procedure chooses the number {labPer} (a
    positive integer) of ordinary tic steps for each major tic step,
    so that the label increment {labPer*ztStep} is a nice round value.
    It also computes the number {lapSkp} of minor tic steps between
    {ztMin} and the first major tic after it.
    
    If {minDist} is positive, it specifies the minimum spacing (in mm)
    between labeled tics. If {minStep} is positive, it specifies the
    minimum increment (in Graph coordinates) between consecutive
    labels.
    
    The results are returned in {*labPerP,*lasbSkpP}.  If the constraints
    imply that there are no major tics in the range {ztMin,ztMax},
    returns {*labPerP=0}. */ 

void msm_ps_tools_draw_ref_frame(msm_ps_tools_t *mps, double R, double G, double B);
  /* Draw a thin black frame around current reference rectangle
   (NOT the current plottable area). The frame will be drawn 
    with color {R,G,B}. */

void msm_ps_tools_draw_plot_frame(msm_ps_tools_t *mps, double R, double G, double B);
  /* Draw a thin black frame around current plottable area. The frame will be drawn 
    with color {R,G,B}. */

void msm_ps_tools_draw_y_polyline
  ( msm_ps_tools_t *mps,
    double xMin, 
    double xMax,
    double y[],
    int32_t n
  );
  /* Draws a polygonal line with {n} vertices and {n-1} line segments.
    The vertex X coordinates are univormly spaced in {[xMin_xMax]}.
    The Y coordinates are {y[0..n-1]}. */

void msm_ps_tools_draw_y_dots
  ( msm_ps_tools_t *mps,
    double xMin, 
    double xMax,
    double y[],
    int32_t n,
    double rad,
    bool_t fill,
    bool_t draw
  );
  /* Draws {n} dots of radius {rad}. The X coordinates of the
    dot centers are uniformly spaced in {[xMin_xMax]}. The Y
    coordinates are {y[0..n-1]}.
    
    The radius {rad} is in mm, while the centers are in Graph
    coordinates. The parameters {fill} and {draw} are as in
    {epswr_dot}. */

/* HIGH-LEVEL DRAWING TOOLS 

  Unless specified otherwise, all procedures below take Graph
  coordinates. */

void msm_ps_tools_draw_graphs
  ( msm_ps_tools_t *mps,
    int32_t nc,
    int32_t nd,
    double x[],
    double start,
    double step,
    double y[],
    double yMin,
    double yMax
  );
  /* Draws {nc} graphs on {mps}, each with a different color. 
  
    The graphs are defined by {nd} points. Point number {i} has
    abscissa {x[i]} and ordinate {y[nd*c + i]} for {i} in {0..nd-1}
    and {c} in {0..nc-1}.
    
    If {x} is NULL, assumes {x[i] == start + i*step} for all {i}. If {x}
    is not null, {skip} and {step} are ignored.
    
    Also draws the axes, tic marks, and a surrounding frame. The
    vertical plot scale is adjusted to include {[yMin _ yMax]}, plus
    tics, tic labels etc.. */

void msm_ps_tools_draw_histogram
  ( msm_ps_tools_t *mps,
    int32_t nd,
    double x[],
    double y[],
    double yMin,
    double yMax
  );
  /* Draws a histogram on {mps}.
    
    Histogram bar number {i} spans the interval {[x[i] _ x[i+1]]} and
    has height {y[i]}, for {i} in {0..nd-1}. Assumes that the sequence {x[0..nd]}
    is increasing. If {x} is NULL, assumes {x[i] == i - 0.5} for all {i}.
    
    Also draws the axes, tic marks, and a surrounding frame. The
    vertical plot scale is adjusted to include {[yMin _ yMax]}, plus
    tics, tic labels etc.. */

/* FOR LOW LEVEL HACKING */

epswr_figure_t *msm_ps_tools_get_eps_figure(msm_ps_tools_t *mps);
  /* Returns the {epswr_figure_t} object underlying {mps}. 
  
    This can be used to draw graphics objects not supported by this
    interface. Any coordinate arguments must be first mapped b {double
    msm_ps_tools_map_{x,y,coord,coords}} before being used as parameters
    to the {epswr} functions.
    
    The {epswr} scale-setting procedures ({epswr_set_device_window}, 
    {epswr_set_client_window}, etc) should not be used, or chaos will 
    result. */

#endif
