/* stmap - tools for reading, plotting, and manipulating street maps */
/* Last edited on 2022-10-20 05:58:29 by stolfi */

#ifndef stmap_H
#define stmap_H

#include <math.h>
#include <stdint.h>

#include <quad.h>
#include <r2.h>
#include <vec.h>

#include <stdio.h>
#include <pswr.h>

typedef struct Interval { double lo, hi; } Interval;

typedef struct RGBColor { float R, G, B; } RGBColor;

typedef r2_t Point;

typedef struct EdgeData
  { float cost[2];  /* Traversal cost for each direction. */
    int32_t id;         /* Undirected edge ID.  */
  } EdgeData;
  /* 
    Data for an undirected edge.  If {a} is a quad_arc_t on the edge,
    then {cost[quad_sym_bit(a)]} the traversal cost for that edge in the
    direction of {a}. */

typedef struct VertexData
  { Point p;      /* Vertex coords (m). */
    int32_t id;       /* Vertex ID. */
    int32_t deg;      /* Count of incident undirected edges. */
  } VertexData;

typedef struct Map
  { int32_t nv;             /* Vertex count. */
    int32_t ne;             /* Edge count. */
    VertexData **vd;    /* Vertex data records. */
    EdgeData **ed;      /* Edge data records. */
    quad_arc_t *out;    /* {out[vi]} is some {quad_arc_t} out of vertex number {vi}. */
    quad_arc_t *along;  /* {along[2*ei+s]} is edge {ei} taken in direction {s}. */
  } Map;
  /* A street map is an undirected graph drawn on the plane. The
    quad_arcs {out[vi]} and {along[2*ei+s]} belong to a quad_edge that
    describes the map's topology. For any integer {ai}, we have 
    {quad_sym_bit(along[ai]) = (ai % 2)}. */

Map *st_map_read(FILE *f);
  /* Reads a street map from the file {f}. Plot map format:

    ----------------------------
      p NAME NV NE
      v 0 XXX YYY NIN NOUT
      v 1 XXX YYY NIN NOUT
      ...
      v NV-1 XXX YYY NIN NOUT
      a ORGV DSTV BLOCKED
      a ORGV DSTV BLOCKED
      ...
      a ORGV DSTV BLOCKED
    ----------------------------

  where NV and NE are the vertex and edge counts, XXX and YYY are
  the vertex coordinates (in meters), NIN and NOUT are the vertex
  in- and out-degrees, ORGV and DSTV are the indices of the 
  edge's origin and destination vertices, and BLOCKED is 1
  of the street is non-transitable, 0 otherwise.

  The default edge traversal cost (in both senses) is +oo if the
  edge is BLOCKED, otherwise it is the Euclidean distance between
  its endpoints. */

/* 
  DISTANCES
  
  This module assumes that the incremental cost of traversing an edge
  {a} is {quad_ldata(a)->cost[quad_sym_bit(a)]}. The cost of a path is by
  definition the sum of the costs of its arcs.
  
  The /{u}-cost/ of a vertex {v} (resp. an arc {a}) is the minimum
  cost of any path that starts at {u} and ends with {v} (resp. {a}).
  The /{u}-cost/ of an undirected edge {e} is the minimum of the
  {u}-costs of its two directed arcs.

  A vertex {v} (resp. arc {a}, edge {e}) is /{d}-reachable/ from a
  vertex {u} if its {u}-cost is at most {d}.
  
  The set of all vertices and edges that are {d}-reachable from a
  vertex {u} is the /{d}-ball of {u}/.  */

void st_map_compute_costs
  ( Map *m, 
    int32_t u, 
    float dMax,
    int32_t *r,
    int32_t *nr, 
    float *d, 
    quad_arc_t *e,
    float *c
  );
  /* Returns in {nr} and {r[0..nr]} the count and the IDs of the
    vertices that are {dMax}-reachable from {u}. Also, for every
    {dMax}-reachable vertex {v}, sets {d[v]} to the {u}-cost of {v},
    and saves in {e[v]} the last arc of an optimum path from {u} to
    {v} (or {NULL_REF} if {u == v}). Also, for each {dMax}-reachable
    arc {m->along[a]}, sets {c[a]} to its {u}-cost.
    
    The vectors {d,e,r} should have size {m.nv}, and the {c} vector
    should have size {2*m.ne}. IMPORTANT: the vectors {d}, {e} and {c}
    must have been initialized with {st_init_costs} below, and they
    should be re-initialized with {st_reset_costs} before each
    subsequent call. */

void st_map_init_costs(Map *m, float *d, quad_arc_t *e, float *c);
  /* Initializes the vectors {d}, {e}, and {c} as expected by
    {st_map_compute_costs}. Specifically, sets {d[v] = INFINITY},
    {e[v] = NULL_REF}, {c[ai] = INFINITY} for all {v} and {ai}. */
    
void st_map_reset_costs(Map *m, int32_t *r, int32_t nr, float *d, quad_arc_t *e, float *c);
  /* Re-initializes the vectors {d}, {e} and {c} after a call
    to {st_map_compute_costs}, in preparation for a new call. Only
    the vertices {r[0..nr-1]} and their incident edges are re-initialized. */

void st_compute_coverage
  ( Map *m, 
    int32_t *u, 
    float *dMax, 
    int32_t n, 
    int32_t *vcover,
    int32_t *ecover
  );
  /* Given a list of vertices {u[0..n-1]} (the /sites/), and
    respective cost bounds {dMax[0..n-1]}, computes for each vertex {v}
    the number {vcover[v]} of {dMax[i]}-balls centered at {u[i]} that 
    contain {v}; and, for each undirected edge {e}, the count
    {ecover[e]} of such balls that contain {e}. */

void st_increment_coverage
  ( Map* m,
    float dMax,
    int32_t *r, 
    int32_t nr, 
    float *d,
    float *c,
    int32_t *vcover, 
    int32_t *ecover
  );
  /* Increments {vcover[v]} and {ecover[e]} for every vertex 
    {v} and every undirected edge {e} in the {dMax}-ball of
    some vertex {u}.
   
    Assumes that {r[0..nr-1]} are the vertices in that ball. Also, for
    every vertex {v} and arc {a} in that ball, assumes that {d[v]} and
    {c[a]} are the respective {u}-costs, as would be returned from
    {st_map_compute_costs}. */

void st_map_plot
  ( PSStream *ps, 
    Map *m, 
    Interval xr, Interval yr, 
    float *vwidth, RGBColor *vcolor,
    float *ewidth, RGBColor *ecolor
  );
  /* Plots the map {m} on the Postscript file {f}. The client is
    responsible for writing {f}'s header and trailer.
    
    If the ranges {xr} and {yr} are neither empty nor infinte, the
    procedure assumes that the current clipping region is contained in
    the rectangle {xr × yr}, and threfore it can safely skip any
    vertices or edges that lie entirely outside that rectangle.
    
    Each undirected edge {ei} is drawn with line width {ewidth[ei]}
    (mm) and color {ecolor[ei]}. Each vertex {vi} with is plotted as a
    dot with diameter {vwidth[vi]} (mm) and color {vcolor[vi]}.
    Elements with zero or negative width/diameter are omitted. Any of
    the arrays {ewidth}, {ecolor}, {vwidth}, and {vcolor} can be NULL,
    in which case suitable defaults are used.  */

int32_t st_map_nearest_vertex(Map *m, Point p);
  /* Returns the ID of the vertex of {m} nearest to {p}. */

void st_map_get_bbox(Map *m, Interval *xr, Interval *yr);
  /* Returns in {*xr} and {*yr} the bounding box of all map vertices. */

/* DATA I/O */

void st_write_double_vec_t(char *name, double_vec_t *d);
  /* Writes to file {name} an arbitrary data vector {d}.
    If {name} is "-", writes to {stdout}. */
  
void st_fwrite_double_vec_t(FILE *wr, double_vec_t *d);
  /* Same as {st_write_double_vec_t}, but writes to the open file {wr}. */

double_vec_t st_read_double_vec_t(char *name);
  /* Reads from file {name} an arbitrary data vector {d}. The
    file must be in the format created by {st_write_double_vec_t}. If
    {name} is "-", reads from {stdin}. */

double_vec_t st_fread_double_vec_t(FILE *rd);
  /* Same as {st_read_double_vec_t}, but reads from the open file {rd}. */

/* INTERVAL HACKS */

void st_interval_widen(Interval *r, double margin);
  /* Widens {*r} by the specified {margin} on both sides */

void st_adjust_rect_shape(Interval *xr, Interval *yr, double tx, double ty);
  /* 
    Widens either {*xr} or {*yr}, as needed, to ensure that 
    their dimensions are in the ratio {tx:ty}. */

#endif
