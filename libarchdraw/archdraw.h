#ifndef archdraw_H
#define archdraw_H

/* Primitives for architectural drawings */
/* Last edited on 2023-02-20 18:15:24 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <values.h>
#include <stdint.h>

#include <epswr.h>
#include <r3.h>
#include <r2.h>
#include <vec.h>
#include <frgb.h>
#include <bool.h>

/* REFERENCE POINTS */

typedef struct adrw_point_t 
  { r3_t p;      /* Point coordinates (cm). */
    char *lab;   /* Point label (possibly NULL). */
    int32_t krot;    /* Direction angle of label, degrees ccw from {+X}. */
    bool_t used; /* {TRUE} if point has been used by a command from this library. */
  } adrw_point_t;
 
bool_t is_defined(r3_t *p);
  /* TRUE if {*p} is a well-defined point. */

/* POINT LISTS */
 
vec_typedef(adrw_point_vec_t, adrw_point_vec, adrw_point_t); 

void adrw_append_point
  ( char *lab, 
    int32_t ip, 
    int32_t jpx, 
    int32_t jpy, 
    int32_t jpz, 
    double dX, 
    double dY, 
    double dZ, 
    adrw_point_vec_t *P,
    int32_t *nP
  );
  /* Defines {P->el[ip]} as point {(P[jpx].x,P[jpy].y,P[jpz].z)+(dX,dY,dZ)}, with label {lab}.
    If {ip} is {*np} or more, also marks all points in {*np..ip-1} as undefined,
    and sets {*np} to {ip+1}. 
    
    If {S.P[ip]} was defined previously, the labels must be identical
    and the coordinates must agree, allowing for some accumulated
    floating point rounding error.  In any case the new definition is 
    ignored.
    
    If {jpx} is negative, {P[jpx].x} is assumed to be 0; that is, {dX} is assumed
    to be an absolute coordinate.  The same rule holds for the other two axes.
    
    All points created are marked with {used=FALSE}. */

adrw_point_t adrw_make_derived_point(adrw_point_t *ap, int32_t k, double dx, double dy, double dz, int32_t krot);
  /* Creates an {adrw_point_t} record that is the point {*ap} 
     displaced by {(dx,dy,dz)}.
    
    If the displacement is {(0,0,0)}, the result is just a copy of
    {*ap}, with the same label text and rotation angle. Otherwise the
    label will be "{ap.lab}.{k}" and the rotation angle will be
    {krot}. */

void show_holes(adrw_point_vec_t *P, int32_t ini, int32_t fin);
  /* Lists all points in {P} in the range {ini..fin} that are still undefined. */
    
bool_t adrw_polygon_is_closed(adrw_point_vec_t *P);
  /* TRUE if polygon {P.e[0..P.ne-1]} is closed in the XY projection. */

void adrw_compute_area_and_length(adrw_point_vec_t *P, double *areaP, double *lengthP);
  /* Computes the 3D length and the XY-projected area of the polygon {P}.
    If the polygon is not closed, the area is zero. To be closed, the 
    points must be the same */

/* UNIT TYPES AND PLOTTING STYLES */ 

typedef uint32_t adrw_space_type_t;

typedef struct adrw_unit_style_t
  { frgb_t fill_rgb;   /* Area fill color. */
    frgb_t draw_rgb;   /* Perimeter stroke color. */
    double pen_width;  /* Stroke line width (true mm). */
    frgb_t dots_rgb;   /* Dot fill color. */
    double dot_radius; /* Dot radius (true mm). */
  } adrw_unit_style_t ;
  /* Style parameters to be used when filing/drawing a room. */

adrw_unit_style_t *adrw_make_unit_style
  ( frgb_t *fill_rgb,
    frgb_t *draw_rgb,
    double pen_width,
    frgb_t *dots_rgb,
    double dot_radius
  );
  /* Allocates a {adrw_unit_style_t} record and fills it with the given parameters.
    If {fill_rgb} is NULL or negative, or polygon is open, assumes no fill.
    If {draw_rgb} is NULL or negative, or {pen_width} is zero, assumes no draw.
    If {dots_rgb} is NULL or negative, dots are not filled.
    If {dot_radius} is non-positive, does not draw dots. */

/* BUILDING UNITS */ 

typedef struct adrw_unit_t
  { char *label;               /* Short label of unit. */
    char *descr;               /* Description of unit (full name, occupants, detailed use, etc.). */
    adrw_point_vec_t pt;       /* Points that define the outline. */
    double round;              /* Corner rounding radius. */
    int32_t kfloor;                /* Floor index (0 = ground). */
    adrw_space_type_t type;    /* Usage type class (for tabulations, etc.). */
    adrw_unit_style_t *style;  /* Style to use when drawing the unit. */
    /* Computed fields: */
    double modules;            /* Number of modules, for tabulation. */
    double area;               /* Area of XY projection (0 if not closed) (cm^2). */
    double length;             /* Length ofpolygonal line or perimeter of closed polygon in 3D  (cm). */
  } adrw_unit_t;
  /* An {adrw_unit_t} is a polygonal region, normally on a single 
    horizontal plane.
    
    In an architectural drawing, usually an office or other room, but
    may also be a hallway, pillar, vertical passage, exterior perimeter
    of a building, etc. It may also be an open polygonal line, e.g. a
    pipe or cable; or a sequence of dots.
    
    The points {pt.e[0..pt.ne-1]} should be the corners. For a closed
    area, the last vertex must be repeated. If {round} is positive, they
    will be rounded off with arcs of that radius.  This is valid
    only if the points lie on the same horizontal plane. */

adrw_unit_t *adrw_make_poly
  ( char *label, 
    char *descr, 
    double modules,
    adrw_point_vec_t *P, 
    int32_t v[], 
    double round,
    int32_t kfloor, 
    adrw_space_type_t type, 
    adrw_unit_style_t *style
  );
  /* Creates a unit that is a polygon or polygonal line
    whose corners are the previously defined points {P.e[v[0..n-1]]},
    where {n} is the first index such that {v[n] < 0}.  For a 
    closed area, specify {v[n-1] == v[0]}. */

adrw_unit_t *adrw_make_box
  ( char *label, 
    char *descr, 
    adrw_point_vec_t *P, 
    int32_t ip, 
    double ctrx,
    double ctry,
    double ctrz,
    double wdx,
    double wdy,
    double round,
    int32_t kfloor, 
    adrw_space_type_t type, 
    adrw_unit_style_t *style
  );
  /* Creates a unit that is a horizontal rectangle, given its center
    and sides {wdx,wdy}. The center coordinates are 
    {P.e[ip]+(ctrx,ctry,ctrz)} if {ip >= 0}, or just {(ctrx,ctry,ctrz)} if {ip}
    is negative. The number of modules is assumed to be 0. */

adrw_unit_t *adrw_make_dot
  ( char *label, 
    char *descr, 
    adrw_point_vec_t *P, 
    int32_t ip, 
    double ctrx,
    double ctry,
    double ctrz,
    int32_t kfloor, 
    adrw_space_type_t type, 
    adrw_unit_style_t *style
  );
  /* Creates a unit that is a single point. The center coordinates are 
    {P.e[ip]+(ctrx,ctry,ctrz)} if {ip >= 0}, or just {(ctrx,ctry,ctrz)} if {ip}
    is negative. The area, length, and number of modules are assumed to be 0. */

vec_typedef(adrw_unit_vec_t, adrw_unit_vec, adrw_unit_t *);
  
typedef struct adrw_building_t
  { int32_t NU;                /* Number of units (incl. pillars, external outline, pipes, etc.) */
    adrw_unit_vec_t unit;  /* Valid rooms are {unit.e[0..NR-1]}. */
  } adrw_building_t;
  /* A {adrw_building_t} is a set of floorplans for a building.
    It is basically a list of units, possibly in different floors 
    or blocks. */

adrw_building_t *adrw_make_building(void);
  /* Creates a building, initally with zero units. */

void adrw_append_unit(adrw_building_t *B, adrw_unit_t *rm);
  /* Appends the unit {rm} to the building planset {B}. */

void adrw_append_seats
  ( adrw_building_t *B, 
    adrw_point_vec_t *P, 
    int32_t v00,
    int32_t v11,
    double szx,
    double szy,
    int32_t kfloor, 
    adrw_space_type_t type, 
    adrw_unit_style_t *style
  );
  /* Appends a series of units representing seats. Each unit is a
    square of side {fmin(szx,szy)} with the given {style}, nested
    inside an invisible /nominal rectangle/. These nominal rectangles
    tile an axis-aligned /seating area/ defined by the opposite corners
    {P.e[v00]} and {P.e[v11]}. The corner {v00} must be at the back of
    the seating area. The nominal rectangles measure at least {szx ×
    szy}, but are expanded as needed to completely cover the seating
    area. All rectangles are counted as 0 modules. */
  
/* PRINTOUT */

void adrw_print_points(FILE *wr, adrw_point_vec_t *P);

void adrw_print_unit
  ( FILE *wr,
    adrw_unit_t *rm,
    char *type_tag[],
    bool_t print_modules,
    bool_t print_area,
    bool_t print_length,
    bool_t TeX
  );
  /* Prints one line of a building-units table.
    The {type_tag} array is used to convert unit types to strings. 
    If any of {print_modules,print_area,print_length} is false,
    omith the corresponing field (and its separator).
    Areas are printed in {m^2}, lengths in {m}. */

void adrw_compute_building_stats
  ( adrw_building_t *B,
    int32_t ntypes,
    double units[],
    double modules[],
    double area[],
    double length[]
  );
  /* Computes the total units, modules, area, and length of units by unit type.
    It any aray is NULL, the corresponding stats are not computed. */

void adrw_print_building
  ( FILE *wr,
    adrw_building_t *B,
    bool_t select[],
    char *type_tag[],
    bool_t print_modules,
    bool_t print_area,
    bool_t print_length,
    bool_t TeX
  );
  /* Prints a table with all units of the building, sorted by type and
    label. If {select} is NULL, prints all units. If {select} is not
    NULL, prints a unit {rm} only if {select[rm.type]} is TRUE. If any
    of {print_modules,print_area,print_length} is false, omith the
    corresponing column of the table. The {type_tag} array is used to
    convert unit types to strings. */

/* PLOTTING */

epswr_figure_t *adrw_new_figure
  ( char *dir,
    char *prefix,
    char *suffix,
    double xmin, 
    double xmax, 
    double ymin, 
    double ymax, 
    int32_t ox, 
    int32_t nx, 
    int32_t oy, 
    int32_t ny, 
    char *title
  );
  /* Starts a new EPS figure, writing to a file 
    "{dir}/{prefix}_{name}_{NNNNN}_{suffix}.eps". 
   
    The whole plot is assumed to span a rectangle {[xmin _ xmax} × [ymin
    _ ymax]} in Client coordinates. This rectangle is divided into {nx}
    by {ny} sub-rectangles. Then the scale is set so that the
    sub-rectangle in column {ox} from left and row {oy} from bottom fits
    into the available plotting area.
    
    The {title} as well as the ranges are printed at the bottom of the
    figure.
    
    Any of the {dir}, {prefix}, {name}, or {suffix} arguments may be
    {NULL} in which case those parts of the name and the associated "/"
    or "_" are omitted. */

void adrw_plot_unit
  ( epswr_figure_t *epsf,
    adrw_unit_t *rm,
    bool_t show_dots
  );
  /* Plots the unit {rm} to {epsf}.  If {show_dots} is true, forces display of vertices in black. */

void adrw_plot_building
  ( epswr_figure_t *epsf,
    adrw_building_t *B,
    bool_t show_dots
  );
  /* Plots the building {B} to {epsf}.  If {show_dots} is true, forces display of vertices in black. */

void adrw_plot_point
  ( epswr_figure_t *epsf, 
    char *lab, 
    int32_t ip, 
    double x, 
    double y, 
    double rot, 
    double hAlign, 
    double vAlign
  );
  /* Plots a dot at coordinates {(x,y)}.  The dot will be labeled {lab} 
    (if {lab != NULL}) or "P{ip}" (if {lab == NULL}).  
    
    The parameters {rot,hAlign,vAlign} define the rotation and placement of 
    the label relative to {(x,y)}, as in {epswr_label}. */

void adrw_plot_type_legend
  ( char *dir,
    char *prefix,
    char *suffix,
    int32_t key_type[], 
    int32_t ntypes,
    int32_t ncols, 
    char *type_tag[], 
    adrw_unit_style_t *style[]
  );
  /* Creates an EPS file called "{dir}/{prefix}_{suffix}.eps" containing
    the color key for the usage types {tp[0..nt-1]}, where {nt} is the
    first index such that {tp[nt] == -1}. The key will have {ncols}
    columns, pairing the color in {style[type[ip]]} with the tag in
    {type_tag[type[ip]]}.
    
    Each entry in the key is 8 mm tall and 24 mm wide, including spacing
    and margins. */

void adrw_plot_histogram_bar
  ( char *dir,
    char *prefix,
    char *name,
    char *suffix,
    adrw_space_type_t type,
    adrw_unit_style_t *st,
    double val,
    double vmax
  );
  /* Writes an EPS file called "{dir}/{prefix}_{name}_{suffix}.eps" that contains a single 
    horizontal histogram bar, with length proportional to {val/valmax}.

    The bar will be filled with the color specified in {st} and 
    its border will be stroked in black. The figure will be 160 
    mm wide and 4 mm tall. The bar, including the stroked border. */

#endif
