#ifndef hxg_paint_H
#define hxg_paint_H

/* Routines to paint figures on an hexagonal grid. */
/* Last edited on 2025-01-01 01:27:27 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <r2.h>
#include <i2.h>
#include <interval.h>

#include <hxg_canvas.h>

/* TIE BREAKING

  In all the procedures below, ties are broken by implicitly shifting
  the figure down by an infinitesimal epsilon, and left by a larger
  infinitesimal delta. */

typedef void hxg_paint_op_t(uint32_t ix, uint32_t iy, uint32_t k);
  /* A painting procedure {Pr} usually takes a pixel operation {op} of
    this type, and will apply it to all pixels inside the specified figure.
    The arguments are the integer column and row indices of the 
    pixel and the linearized pixel index {k = ix + nx*iy}, where
    {nx} is the number of pixels per row. */ 

void hxg_paint_canvas(hxg_canvas_t *cvs, hxg_paint_op_t *op);
  /* applies {op} to the whole canvas {cvs}. */

void hxg_paint_polygon
  ( hxg_canvas_t *cvs,   /* Where to paint. */
    r2_vec_t *p,          /* Vertices of polygon, ccw around interior. */
    hxg_paint_op_t *op   /* The pixel painting operation. */
  );
  /* Applies {op} to all canvas pixels that lie inside the polygon {p}.
    
    The polygon is defined by the ordered list of its vertices; an edge
    is assumed between every two consecutive vertices. The list
    must be explicitly closed; that is, the last vertex 
    must be identical to the first one.
    
    The interior is defined by the even/odd rule. To represent
    holes and multiple components, just concatenate their vertex lists,
    each explicitly closed as above, separated by a dummy point with 
    infinite coordinates {# = (INF,INF)}. That is, for a polygon 
    with two islands with vertices {A,B,C,D} and {M,N,O}, and a hole
    {U,V,W}, the vector {p} should contain {A,B,C,D,A,#,M,N,O,M,#,U,V,W}. */
  
void hxg_paint_rect_stroke
  ( hxg_canvas_t *cvs,  /* Where to paint. */
    r2_t *a,            /* Beginning of stroke axis. */
    r2_t *b,            /* End of stroke axis. */
    double r,           /* Half-width of stroke. */
    hxg_paint_op_t *op  /* The pixel painting operation. */
  );
  /* Applies {op} to all canvas pixels that lie inside the rectangle
    whose axis is the segment {a,b} and whose half-width (extent
    perpendicular to that axis, in each direction) is {r}. Equivalent
    to painting a polygon with the appropriate vertices. */

void hxg_paint_circle
  ( hxg_canvas_t *cvs,  /* Where to paint. */
    r2_t *c,            /* Center of circle. */
    double r,           /* Radius of circle. */
    hxg_paint_op_t *op  /* The pixel painting operation. */
  );
  /* Paints all canvas pixels that lie inside the circle
    with center {c} and radius {r}. */
  
void hxg_paint_sausage
  ( hxg_canvas_t *cvs,  /* Where to paint. */
    r2_vec_t *p,         /* Vertices of sausage(s). */
    double r,           /* Half-width of sausage. */
    hxg_paint_op_t *op  /* The pixel painting operation. */
  );
  /* Paints all canvas pixels that are at most {r} away from
    any point of a set of polygonal lines, specified by their vertices
    {p[0..n-1]}.
    
    The lines need not be closed, and must be separated from each
    other by the dummy vertex {# = (INF,INF)}. It is equivalent to
    painting a circle of radius {r} at every finite vertex, and a
    rectangular stroke of half-width {r} between every two consecitive
    finite vertices.  
    
    Beware that, unlike most other routines from this module, this
    procedure may call {op} two or more times for the same pixel! */

#endif
