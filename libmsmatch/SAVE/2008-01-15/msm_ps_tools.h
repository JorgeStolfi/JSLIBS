#ifndef msm_ps_tools_H
#define msm_ps_tools_H

/* Postscript plots of graphs and such. */
/* Last edited on 2008-01-11 12:09:58 by stolfi */

#define msm_ps_tools_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>

#include <vec.h>
#include <pswr.h>

#define msm_PT_PER_MM (72.0/25.4)
  /* Postscript points (pt) per millimeter (mm). */ 

#define msm_EPS_MARGIN_MM (1.0)
  /* EPS figure margin width (in mm). */

/* POSTCRIPT FIGURE STREAMS */

typedef struct msm_ps_tools_t msm_ps_tools_t;
  /* A handle to a Postscript graph-plotting stream.
  
    An {msm_ps_tools_t} is a handle to a {FILE} descriptor plus some
    internal state. The file is supposed to contain one Encapsulated
    Postscript image.
    
    Procedures that create an {msm_ps_tools_t} take as parameters a
    {FILE*} {wr} and two strings {name} and {tag}. If {wr} is not
    NULL, the figure will be written to that file, and the {name} and
    {tag} are ignored.
    
    If {wr} is NULL, the procedure will open a file called
    "{name}{tag}.eps" and use that file instead of {wr}. In both
    cases, the procedure will write the appropriate EPS preamble to
    the file.
    
    All dimensional parameters are in mm, except where said otherwise.
    
    The client must call {msm_ps_tools_close} at the end, to write the necessary 
    postamble and close the underlying {FILE}. 
  */
 
msm_ps_tools_t *msm_ps_tools_new
  ( FILE *wr,
    char *name,
    char *tag,
    double hSize, 
    double vSize,
    double mrg
  );
  /* Creates a {msm_ps_tools_t} for an Encapsulated Postscript
    figure of width {hSize} and height {vSize}, plus an extra margin
    of width {mrg} all around. All dimensions are in mm.
    
    The client's coordinate system and plotable area will be set
    initially to the rectangle {[0_hSize]×[0_vSize]] excluding the
    margin. */

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
    int maxLabChars,
    double mrg
  );
  /* Creates a {msm_ps_tools_t} for an Encapsulated Postscript
    figure that is to contain a graph of width {hGraphSize} and height
    {vGraphSize}.
    
    The plottable area will include space for the graph, plus space
    around it for graph scales at the left, right, bottom or top (as
    specified by the flags {scaleL,scaleR,scaleB,scaleT} and/or titles
    (as specified by the flags {titleL,titleR,titleB,titleT}. The
    figure will have an extra margin of {mrg} mm all around the
    plottable area.  The width of the bands is estimated 
    from the {fontSize} parameter (the label font's nominal
    size, in pt) and {maxLabChars} (the max number of characters 
    in any label of the scale). 
    
    The client and device reference rectangles will be initialized
    to the graph area only (excluding the scale and title areas). */
  
void msm_ps_tools_close(msm_ps_tools_t *mps);
  /* Terminates any pictures that have been written to {mps}, flushes
    and closes the underlying file (unless it is {stdout} or
    {stderr}). Then frees all internal storage used by {mps}, including
    the associated {PSStream} object and {*mps} itself. */

PSStream *msm_ps_tools_get_ps_stream(msm_ps_tools_t *mps);
  /* Returns the {PSStream} object associated with the {msm_ps_tools_t}
    object {mps}.  The {pswr.h} routines can be used to draw on that
    stream. For those routines, the client coordinates are
    measured (in mm) from the lower left corner of the usable
    plotting area. */

void msm_ps_tools_get_plot_size(msm_ps_tools_t *mps, double *hSize, double *vSize);
  /* Stores in {*hSize} and {*vSize} the dimensions (in mm) of the usableplot
    area, as specified to {msm_ps_tools_init}, excluding the margin. */
   
/* COORDINATE MAPPING
  
  Many plotting routines in this interface take their arguments in
  /client coordinates/ {x,y}. A {msm_ps_tools_t} object defines
  a mapping of those client coordinates to /device
  coordinates/ {h,v}, measured in millimeters from the lower left
  corner of the plot area. */
  
double msm_ps_tools_map_x(msm_ps_tools_t *mps, double x);
double msm_ps_tools_map_y(msm_ps_tools_t *mps, double y);
double msm_ps_tools_map_coord(msm_ps_tools_t *mps, pswr_axis_t axis, double coord);
void msm_ps_tools_map_coords(msm_ps_tools_t *mps, double x, double y, double *h, double *v);
  /* These procedures map client coordinates {x,y} to device
     coordinates {h,v}, as defined by {mps}. */

double msm_ps_tools_unmap_h(msm_ps_tools_t *mps, double h);
double msm_ps_tools_unmap_v(msm_ps_tools_t *mps, double v);
double msm_ps_tools_unmap_coord(msm_ps_tools_t *mps, pswr_axis_t axis, double coord);
void msm_ps_tools_unmap_coords(msm_ps_tools_t *mps, double h, double v, double *x, double *y);
  /* These procedures map device coordinates {h,v} to client
     coordinates {x,y}, as defined by {mps}. */
  
/* SETTING THE COORDINATE MAPPINGS

  The client-to-device mapping of a {msm_ps_tools_t} object is defined by
  two rectangles stored in the object: the /device reference window/
  {DRW = [hMin_hMax]×[vMin_vMax]} and the /client reference window/
  {CRW = [xMin_xMax]×[yMin_yMax]}.
  
  The client-to-device mapping is such that the {CRW} is mapped to the
  {DRW}. More precisely, {x=xMin} is mapped to {h=hMin}, {x=xMax} to
  {h=hMax}. So, if {hMin > hMax} (or {xMin > xMax}) the client {x}
  axis is directed from right to left. The same holds for the {y} and
  {v} coordinates.
  
  Changing either of these rectangles affects the way that client
  coordinates are interpreted by subsequent function calls.

  The rectangles {CRW} and {DRW} are not used for clipping.
  Regardless of them, the usable plotting area is always
  the rectangle {[0_hSize]×[0_vSize]} in device coordinates. */

void msm_ps_tools_set_device_ref_window
  ( msm_ps_tools_t *mps,
    double hMin, double hMax, 
    double vMin, double vMax
  );
  /* Sets the device reference window to the specified
    rectangle (in mm from the lower left corner of the 
    plottable area).*/

void msm_ps_tools_set_client_ref_window
  ( msm_ps_tools_t *mps,
    double xMin, double xMax, 
    double yMin, double yMax
  );
  /* Sets the client reference window to the specified
    rectangle (in client coordinates).*/

void msm_ps_tools_shrink_device_ref_window
  ( msm_ps_tools_t *mps,
    double lMrg, double rMrg, 
    double bMrg, double tMrg
  );
  /* Displaces the device reference window boundary INWARDS by the
    specified amounts on the left, right, bottom, and top sides,
    respectively. */

void msm_ps_tools_expand_client_ref_window
  ( msm_ps_tools_t *mps,
    double lMrg, double rMrg, 
    double bMrg, double tMrg
  );
  /* Displaces the client reference window boundary OUTWARDS by the
    specified amounts on the left, right, bottom, and top sides,
    respectively. */

void msm_ps_tools_compute_data_range(int n, int stride, double z[], double *zMinP, double *zMaxP);
  /* Computes a data range {[*zMinP _ *zMaxP]} appropriate for plotting the 
    graph of {n} given values, namely {z[i*stride]} for {i = 0..n-1]}.
    The resulting range will be non-empty and non-degenerate even if {n == 0}
    or all {z} values are equal. */

/* LOW-LEVEL DRAWING TOOLS 

  Unless specified otherwise, all procedures below take client
  coordinates. */

void msm_ps_tools_draw_segment(msm_ps_tools_t *mps, double xa, double ya, double xb, double yb);
  /* Draws a line segment from {(xa,ya)} to {(xb,yb)} with the current pen
    and draw color. */

void msm_ps_tools_draw_ref_axis(msm_ps_tools_t *mps, pswr_axis_t axis);
  /* Draw the {x} or {y} coordinate axes in dark gray.
    The axis will span the current reference rectangle 
    (NOT the current plottable area). */

void msm_ps_tools_draw_tic
  ( msm_ps_tools_t *mps, 
    pswr_axis_t axis, 
    double xt, 
    double yt,
    double ticSize,
    double ticAlign,
    char *label,
    double labAlign
  );
  /* Draws a tic mark, perpendicular to the given {axis}, at the
    nominal position {xt,yt}, with the current pen and draw color.
    The parameters {ticSize,ticAlign} are
    passed to {pswr_tic}.
    
    If the string {label} is not NULL, its is written using the
    current font and text color. The label is centered on the tic
    position along the given {axis}, and aligned as specified by
    {labAlign}, relative to the tic, along the other axis. Thus, for
    example, {labAlign=1.2} would be appropriate for a tic along the
    left or bottom edges of the plot area, and {labAlign=-0.2} would
    do for tics along the top and right edges. */

void msm_ps_tools_draw_scale
  ( msm_ps_tools_t *mps, 
    pswr_axis_t axis, 
    double pos,
    double ticSize,
    double ticAlign,
    double ticMinDist, 
    double ticMinStep, 
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
    {ticMinStep} in client coordinates. (At least one of these
    constraints must be positive.)
    
    If the parameter {fmt} is not NULL, a subset of the tics will be
    labeled with the corresponding client coordinates, in the format
    {fmt} (as per {printf}). Labeled tics are longer than unlabeled
    ones ({1.5*ticSize}, in mm). The spacing between labeled tics
    will be at least {labMinDist} in millimeters, and {labMinStep}
    in client coordinates.  
    
    The parameters {ticSize,ticAlign,labAlign} are explained
    under {msm_ps_tools_draw_tic}. */

void msm_choose_label_coords
  ( msm_ps_tools_t *mps,
    pswr_axis_t axis,
    double ztMin,
    double ztMax, 
    double ztStep,
    double minDist,
    double minStep,
    int *labPerP,
    int *labSkpP
  );;
  /* Chooses client coordinates for major (labeled) tics,
    given the spacing {ztStep} (in client coordinates) between
    minor tics, and the client coordinates {ztMin,ztMax} of the first and last
    minor tics.

    Specifically, the procedure chooses the number {labPer} (a
    positive integer) of ordinary tic steps for each major tic step,
    so that the label increment {labPer*ztStep} is a nice round value.
    It also computes the number {lapSkp} of minor tic steps between
    {ztMin} and the first major tic after it.
    
    If {minDist} is positive, it specifies the minimum spacing (in mm)
    between labeled tics. If {minStep} is positive, it specifies the
    minimum increment (in client coordinates) between consecutive
    labels.
    
    The results are returned in {*labPerP,*lasbSkpP}.  If the constraints
    imply that there are no major tics in the range {ztMin,ztMax},
    returns {*labPerP=0}. */ 
    
void msm_ps_tools_choose_tic_coords
  ( msm_ps_tools_t *mps, 
    pswr_axis_t axis, 
    double cMin, 
    double cMax, 
    double minDist, 
    double minStep,
    double *zMinP,
    double *zMaxP,
    double *zStepP
  );
  /* Chooses client coordinates for minor tics within the range
    {[cMin_cMax]} of device coordinates.
    
    Specifically, the procedure chooses a nice round increment
    {zStep} in client coordinates between two consecutive tics, and
    computes the client coordinates {zMin,zMax}, multiples of
    {zStep}, of the first and last tic that fit in that range.
    
    The parameter {minDist}, if positive, and specifies the minimum
    distance (in mm) between tics. The parameter {minStep}, if
    positive, specifies the minimum client coordinate increment
    between tics. At least one of {minDist,minStep} must be
    positive.
    
    The results are returned in {*zMinP,*zMaxP,*zStepP}.
    If there are no nice coordinate values in the specified range,
    returns {*zStepP=0}, {*zMinP>*zMaxP}. */ 

void msm_psplot_choose_label_coords
  ( msm_ps_tools_t *mps,
    pswr_axis_t axis,
    double ztMin,
    double ztMax, 
    double ztStep,
    double minDist,
    double minStep,
    int *labPerP,
    int *labSkpP
  );
  /* Chooses client coordinates for major (labeled) tics,
    given the spacing {ztStep} (in client coordinates) between
    minor tics, and the client coordinates {ztMin,ztMax} of the first and last
    minor tics.

    Specifically, the procedure chooses the number {labPer} (a
    positive integer) of ordinary tic steps for each major tic step,
    so that the label increment {labPer*ztStep} is a nice round value.
    It also computes the number {lapSkp} of minor tic steps between
    {ztMin} and the first major tic after it.
    
    If {minDist} is positive, it specifies the minimum spacing (in mm)
    between labeled tics. If {minStep} is positive, it specifies the
    minimum increment (in client coordinates) between consecutive
    labels.
    
    The results are returned in {*labPerP,*lasbSkpP}.  If the constraints
    imply that there are no major tics in the range {ztMin,ztMax},
    returns {*labPerP=0}. */ 

void msm_ps_tools_draw_ref_frame(msm_ps_tools_t *mps);
  /* Draw a thin black frame around current reference rectangle
   (NOT the current plottable area). */

void msm_ps_tools_draw_plot_frame(msm_ps_tools_t *mps);
  /* Draw a thin black frame around current plottable area. */

void msm_ps_tools_draw_y_polyline
  ( msm_ps_tools_t *mps,
    double xMin, 
    double xMax,
    double y[],
    int n
  );
  /* Draws a polygonal line with {n} vertices and {n-1} line segments.
    The vertex X coordinates are univormly spaced in {[xMin_xMax]}.
    The Y coordinates are {y[0..n-1]}. */

void msm_ps_tools_draw_y_dots
  ( msm_ps_tools_t *mps,
    double xMin, 
    double xMax,
    double y[],
    int n,
    double rad,
    bool_t fill,
    bool_t draw
  );
  /* Draws {n} dots of radius {rad}. The X coordinates of the
    dot centers are uniformly spaced in {[xMin_xMax]}. The Y
    coordinates are {y[0..n-1]}.
    
    The radius {rad} is in mm, while the centers are in client
    coordinates. The parameters {fill} and {draw} are as in
    {pswr_dot}. */

/* HIGH-LEVEL DRAWING TOOLS 

  Unless specified otherwise, all procedures below take client
  coordinates. */

void msm_ps_tools_draw_graphs
  ( msm_ps_tools_t *mps,
    int nc,
    int nd,
    bool_t circ,
    double x[],
    double y[],
    double yMin,
    double yMax
  );
  /* Draws {nc} graphs on {mps}, each with a different color. 
  
    The graphs are defined by {nd} points. Point number {i} has
    abscissa {x[i]} and ordinate {y[nd*c + i]} for {i} in {0..nd-1}
    and {c} in {0..nc-1}.  If {x} is NULL, assumes {x[i] == i}
    for all {i}.
    
    If {circ} is true the graphs are assumed to be
    periodic, and an extra segment is drawn at each end of the plot to
    show that.  In that case the vector {x} must have {nd+1} elements,
    and the period is assumed to be {x[nd] - x[0]}.
    
    Also draws the axes, tic marks, and a surrounding frame. The
    vertical plot scale is adjusted to include {[yMin _ yMax]}, plus
    tics, tic labels etc.. */

void msm_ps_tools_draw_histogram
  ( msm_ps_tools_t *mps,
    int nd,
    double x[],
    double y[],
    double yMin,
    double yMax
  );
  /* Draws a histogram on {mps}.
    
    Histogram bar number {i} spans the interval {[x[i] _ x[i+1]]} and
    has height {y[i]}, for {i} in {0..nd-1}. Assumes that {x[0..nd]}
    are increasing. If {x} is NULL, assumes {x[i] == i - 0.5} for all {i}.
    
    Also draws the axes, tic marks, and a surrounding frame. The
    vertical plot scale is adjusted to include {[yMin _ yMax]}, plus
    tics, tic labels etc.. */

#endif
