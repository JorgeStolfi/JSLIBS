#ifndef hxg_canvas_H
#define hxg_canvas_H

/* Pixel arrays with hexagonal geometry. */
/* Last edited on 2025-01-01 03:04:09 by stolfi */

#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <interval.h>

#include <cpk_basic.h>

#define hxg_canvas_MAX_SIZE 2048
/* Maximum width or height of an hexagonal canvas. */

typedef struct hxg_canvas_t 
  { r2_t org;          /* Position of lower left grid point. */
    r2_t step;         /* Grid steps. */
    uint32_t size[2];  /* Number of grid points along each axis. */    
  } hxg_canvas_t;
  /* A {hxg_canvas_t} is an array of "pixels", associated with points
    of an hexagonal grid. The bottom left point of the grid is {org};
    points are spaced {step[0]} along each horizontal row, and each
    row is spaced {step[1]} vertically from the previous one. Both
    steps are positive. Even- and odd-numbered rows are displaced
    horizontally by {-0.25*step[0]} and {+0.25*step[0]}, respectivelly.
    Each row has {size[0]} points, and there are {size[1]} rows.
    
    Storage area for the pixels should be provided separately by
    clients. The standard way is to provide a vector {pix[0..N-1]}
    where {N = size[0]*size[1]}, where the pixels are stored by
    scanlines -- i.e. pixel {ix,iy} is stored at {pix[ix +
    size[0]*iy]}. */
    
hxg_canvas_t* hxg_canvas_new(uint32_t nx, uint32_t ny, double xStep, bool_t centered);
  /* Creates a new canvas with {ny} rows of {nx} pixels each. Both
    dismensions must be at least 2. The horizontal step (which must be
    positive) will be {xStep} and the vertical step will be
    {xSet*sqrt(3)/2}, so that the grid cells will be equilateral
    triangles. 
    
    If {centered} is false, the origin wil be the midpoint between
    pixels 0 and 1 in row 0. If {cenetered} is true, origin is set so
    that the bounding pox of all the pixels, considering the odd-even
    shifts, is symmetrical about {(0,0)}. */

r2_t hxg_canvas_pixel_pos(hxg_canvas_t *cvs, i2_t p);
  /* Returns the nominal coordinates of the point number {p.c[0]} on
    row {p.c[1]} of the grid of canvas {cvs}. */

void hxg_canvas_bbox(hxg_canvas_t *cvs, interval_t B[]);
  /* For {k} in {0..1}, sets {B[k]} to the range of coordinates of the 
    pixels of {cvs} along coordinate axis {k}, considering the
    even/odd row shift). */

/* 
  ENUMERATION OF THE INFINITE GRID
  
  The following procedures consider the infinite hex grid underlying 
  the canvas {cvs}, ignoring its bounds. */

i2_vec_t hxg_canvas_half_disk_arcs(hxg_canvas_t *cvs, i2_t c, double r);      
  /* Returns a list of integer vectors {D[k]} such that the pixel with
    indices {p + D[k]} lies at distance less than {r} from the pixel
    with indices {c}. Lists all pairs {D[k] = (ex,ey)} such that {ey >
    0}, or {ey = 0} and {ex > 0}. Note that the set {D} varies
    depending on whether {c} lies on an  even or odd scanline. */

double hxg_canvas_sum_over_circle
  ( hxg_canvas_t *cvs, 
    r2_t c, 
    double r, 
    double f(int32_t ix, int32_t iy, r2_t *p),
    r2_t *p
  );
  /* Computes the sum of {f(ix,iy,p)} over all grid pixels {ix,iy} that lie
    inside the circle with center at the pixel {c} and radius {r}. */

#endif
