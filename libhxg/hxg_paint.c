/* See hxg_paint.h */
/* Last edited on 2025-01-01 01:40:16 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <r2.h>
#include <i2.h>
#include <interval.h>
#include <affirm.h>

#include <cpk_basic.h>
#include <hxg_canvas.h>
#include <hxg_paint.h>

/* INTERNAL PROTOTYPES */

int32_t cmpy (r2_t *a, r2_t *b);
  /* Compares points {a} and {b} in Y-coordinate, then X: */

void hxg_scanline_range
  ( double yMin,  /* Minimum Y of object. */
    double yMax,  /* Maximum Y of object. */
    double oy,    /* Ordinate of scanline 0. */
    double sy,    /* Spacing between scanlines. */
    uint32_t ny,       /* Number of scanlines. */
    uint32_t *iMin,    /* OUT: Lowest scanline. */
    uint32_t *iMax     /* OUT: Highest scanline. */
  );
  /* Computes the range {iMin..iMax} of indices of 
    all scanlines in {0..ny-1} whose ordinates intercept
    the interval {[yMin _ yMax]}.  The result may be 
    an empty interval, that is, {iMin > yiMax};
    otherwise we have {0 <= iMin <= iMax <= ny-1}. */

void hxg_add_edge_crossings
  ( r2_t *u,               /* Edge start. */
    r2_t *v,               /* Edge end. */
    double oy,             /* Ordinate of scanline 0. */
    double sy,             /* Spacing between scanlines. */
    uint32_t ny,                /* Number of scanlines. */
    uint32_t ncross[],          /* Count of rdgr crossings per scanline. */
    double_vec_t xcross[]  /* Crossing abscissas per scanline. */
  );
  /* Computes all crossings of edge {(u,v)} with the scanlines of the
    canvas. For every scanline {i} that is crossed, stores the
    crossing abscissa into {xcross[i][n]}, where {n = ncross[i]}; and
    increments {ncross[i]}. The list {xcross[i]} is allocated if {n}
    is zero, otherwise it is extended as needed. */

/* IMPLEMENTATIONS */

int32_t cmpy (r2_t *a, r2_t *b)
  { 
    if (Y(*a) < Y(*b)) 
      { return -1; }
    else if (Y(*a) > Y(*b))
      { return +1; }
    else if (X(*a) < X(*b)) 
      { return -1; }
    else if (X(*a) > X(*b))
      { return +1; }
    else
      { return 0; }
  }

void hxg_paint_canvas(hxg_canvas_t *cvs, hxg_paint_op_t *op)
  { uint32_t nx = cvs->size[0]; /* Num of points per canvas row. */
    uint32_t ny = cvs->size[1]; /* Num of rows in canvas. */
    for (uint32_t iy = 0; iy < ny; iy++)
      { uint32_t k = iy*nx; 
        for (uint32_t ix = 0; ix < nx; ix++) { op(ix, iy, k); k++; }
      }
  }

void hxg_paint_polygon
  ( hxg_canvas_t *cvs,  /* Where to paint. */
    r2_vec_t *p,        /* Vertices of polygon, ccw around interior. */
    hxg_paint_op_t *op  /* Pixel operator. */
  )
  { /* Get hold of params: */
    uint32_t np = p->ne;      /* Number of vertices in list. */
    uint32_t nx = cvs->size[0]; /* Num of points per canvas row. */
    uint32_t ny = cvs->size[1]; /* Num of rows in canvas. */
    double ox = X(cvs->org), oy = Y(cvs->org);    /* Lower left corner of grid. */
    double sx = X(cvs->step), sy = Y(cvs->step);  /* Steps of canvas's grid. */
    /* Find range {iyMin..iyMax} of scanlines that are affected: */
    interval_t B[2];
    r2_bbox(p->ne, p->e, B, TRUE);
    uint32_t iyMin, iyMax;
    hxg_scanline_range(LO(B[1]), HI(B[1]), oy, sy, ny, &iyMin, &iyMax);
    /* Make a list of all edge crossings for each scanline: */
    uint32_t ncross[ny];           /* Number of crossings in each scanline. */
    double_vec_t xcross[ny];  /* Abscissas of crossings in each scanline. */
    for (uint32_t iy = iyMin; iy <= iyMax; iy++) { ncross[iy] = 0; }
    /* Scan the polygon's loops: */
    { uint32_t ip = 0;
      while (ip < np)
        { /* Scan loop starting at {p[ip]}, including final {#} if any: */
          r2_t *pf = &(p->e[ip]);
          demand(r2_is_finite(pf), "unexpected infinite vertex");
          ip++;
          r2_t *pa = pf;
          while (ip < np)
            { r2_t *pb = &(p->e[ip]);
              if (! r2_is_finite(pb)) { ip++; break; }
              hxg_add_edge_crossings(pa, pb, oy, sy, ny, ncross, xcross);
              pa = pb; 
              ip++;
            }
          /* Close polygon if not closed: */
          demand(r2_eq(pa, pf), "non-closed polygon loop"); 
        }
    }
    /* Now fill each scanline according to its crossings: */
    { for (uint32_t iy = iyMin; iy <= iyMax; iy++)
        { /* Get the crossings of this scanline: */
          uint32_t nc = ncross[iy];
          if (nc != 0)
            { assert(nc % 2 == 0); /* Jordan's curve theorem. */
              double *cx = xcross[iy].e;
              /* Sort crossings (insertion sort should be fast enough): */
              for (int32_t ic = 1; ic < nc; ic++)
                { double cxi = cx[ic];
                  int32_t jc = ic; 
                  while ((jc > 0) && (cx[jc-1] > cxi))
                    { cx[jc] = cx[jc-1]; jc--; }
                  cx[jc] = cxi;
                }
              /* Process crossing by even-odd rule: */
              for (int32_t ic = 0; ic < nc; ic += 2)
                { /* Get next "inside" interval: */
                  double xMin = cx[ic], xMax = cx[ic+1]; 
                  /* Get range {ixMin..ixMax} of pixels contained in interval: */
                  /* We assume that the polygon is shifted left by delta: */
                  int32_t ixMin = (int32_t)ceil((xMin - ox)/sx - (iy & 1)/2.0);
                  int32_t ixMax = (int32_t)ceil((xMax - ox)/sx - (iy & 1)/2.0) - 1;
                  /* Intersect {ixMin..ixMax} with {0..nx-1}: */
                  if (ixMin < 0) { ixMin = 0; }
                  if (ixMax > nx-1) { ixMax = (int32_t)nx-1; }
                  if (ixMin <= ixMax)
                    { /* Paint pixels: */
                      uint32_t k = iy*nx + (uint32_t)ixMin;
                      for (uint32_t ix = (uint32_t)ixMin; ix <= (uint32_t)ixMax; ix++) { op(ix, iy, k); k++; }
                    }
                }
              /* Free list of crossings: */
              free(xcross[iy].e);
            }
        }
    }
  }       
  
void hxg_scanline_range
  ( double yMin,  /* Minimum Y of object. */
    double yMax,  /* Maximum Y of object. */
    double oy,    /* Ordinate of scanline 0. */
    double sy,    /* Spacing between scanlines. */
    uint32_t ny,       /* Number of scanlines. */
    uint32_t *iyMin,   /* OUT: Lowest scanline. */
    uint32_t *iyMax    /* OUT: Highest scanline. */
  )
  { /* Compute range {iyMin..iyMax} of scanlines crossed by edge. */
    int32_t jyMin = (int32_t)ceil((yMin - oy)/sy);
    int32_t jyMax = (int32_t)ceil((yMax - oy)/sy)-1;
    /* Intersect {iyMin..iyMax} with {0..ny-1}: */
    if (jyMin < 0) { jyMin = 0; }
    if (jyMax > ny-1) { jyMax = (int32_t)ny-1; }
    if (jyMin > jyMax) { jyMin = 1; jyMax = 0; }
    (*iyMin) = (uint32_t)jyMin; 
    (*iyMax) = (uint32_t)jyMax;
  }

void hxg_add_edge_crossings
  ( r2_t *u,               /* Edge start. */
    r2_t *v,               /* Edge end. */
    double oy,             /* Ordinate of scanline 0. */
    double sy,             /* Spacing between scanlines. */
    uint32_t ny,           /* Number of scanlines. */
    uint32_t ncross[],     /* Count of rdgr crossings per scanline. */
    double_vec_t xcross[]  /* Crossing abscissas per scanline. */
  )
  { /* Make sure {u} is the lower endpoint: */
    if (Y(*u) > Y(*v)) { r2_t *t = u; u = v; v = t; }
    /* Compute range {iyMin..iyMax} of scanlines crossed by edge. */
    uint32_t iyMin, iyMax;
    hxg_scanline_range(Y(*u), Y(*v), oy, sy, ny, &iyMin, &iyMax);
    if (iyMin <= iyMax)
      { /* Get edge's slope {Sx/Sy}: */
        double slope = (X(*v) - X(*u))/(Y(*v) - Y(*u)); 
        /* Add crossings to the crossings lists: */
        for (uint32_t iy = iyMin; iy <= iyMax; iy++)
          { double y = oy + iy*sy;
            double x = X(*u) + slope*(y - Y(*u));
            uint32_t *nc = &(ncross[iy]);
            double_vec_t *cr = &(xcross[iy]);
            if (*nc == 0) 
              { *cr =  double_vec_new(6); }
            else
              { double_vec_expand(cr, (vec_index_t)*nc); }
            cr->e[*nc] = x;
            (*nc)++;
          }
      }
  }

void hxg_paint_circle
  ( hxg_canvas_t *cvs, /* Where to paint. */
    r2_t *c,         /* Center of circle. */
    double r,        /* Radius of circle. */
    hxg_paint_op_t *op        /* Pixel operator. */
  )
  { /* Get hold of params: */
    uint32_t nx = cvs->size[0]; /* Num of points per canvas row. */
    uint32_t ny = cvs->size[1]; /* Num of rows in canvas. */
    double ox = X(cvs->org), oy = Y(cvs->org);    /* Lower left corner of grid. */
    double sx = X(cvs->step), sy = Y(cvs->step);  /* Steps of canvas's grid. */
    double cx = X(*c), cy = Y(*c);
    double r2 = r*r;
    /* Compute range {iyMin..iyMax} of scanlines crossed by circle. */
    double yMin = cy - r;
    double yMax = cy + r;
    uint32_t iyMin, iyMax;
    hxg_scanline_range(yMin, yMax, oy, sy, ny, &iyMin, &iyMax);
    if (iyMin > iyMax) { return; }
    /* fprintf(stderr, "  c = (%f,%f) r = %f iy = %d..%d", cx, cy, r, iyMin, iyMax); */
    /* Enumerate relevant scanlines: */
    for (uint32_t iy = iyMin; iy <= iyMax; iy++)
      { /* Displacement of scanline {iy} from center. */
        double dy = oy + iy*sy - cy;
        /* Find half-extend {hx} of circle on scanline {iy}: */
        double dy2 = dy*dy;
        double hx = (dy2 > r2 ? 0.0 : sqrt(r2 - dy2));
        /* Find span {xMin,xMax} of circle on scanline {iy}: */
        double xMin = cx - hx, xMax = cx + hx; 
        /* Get range {ixMin..ixMax} of pixels contained in interval: */
        /* We assume that the circle is infinitesimally shifted to the left: */
        int32_t jxMin = (int32_t)ceil((xMin - ox)/sx - (iy & 1)/2.0);
        int32_t jxMax = (int32_t)ceil((xMax - ox)/sx - (iy & 1)/2.0) - 1;
        /* Intersect {jxMin..jxMax} with {0..nx-1}: */
        if (jxMin < 0) { jxMin = 0; }
        if (jxMax > nx-1) { jxMax = (int32_t)nx-1; }
        if (jxMin > jxMax) { jxMin = 1; jxMax = 0; }
        uint32_t ixMin = (uint32_t)jxMin;
        uint32_t ixMax = (uint32_t)jxMax;
        /* Paint pixels: */
        uint32_t k = iy*nx + ixMin; 
        for (uint32_t ix = ixMin; ix <= ixMax; ix++) { op(ix, iy, k); k++;  }
      }
  }

void hxg_paint_rect_stroke
  ( hxg_canvas_t *cvs, /* Where to paint. */
    r2_t *a,           /* Vertex of triangle. */
    r2_t *b,           /* Vertex of triangle. */
    double r,          /* Half-width of stroke. */
    hxg_paint_op_t *op        /* Pixel operator. */
  )
  { /* Zero-area figures contain no pixels: */
    if (r == 0) { return; }
    /* Get the displacement vector from {a} to {b}: */
    r2_t d; r2_sub(b, a, &d);
    /* Zero-area figures contain no pixels: */
    if ((X(d) == 0) && (Y(d) == 0)) { return; }
    /* Get its direction: */
    r2_dir(&d, &d);
    /* Get a vector {v} with length {r} perpendicular to {d}: */
    r2_t v = (r2_t){{-r*Y(d), r*X(d)}};
    /* Build the polygon: */
    r2_vec_t p = r2_vec_new(5);
    r2_add(a, &v, &(p.e[0]));
    r2_add(b, &v, &(p.e[1]));
    r2_sub(b, &v, &(p.e[2]));
    r2_sub(a, &v, &(p.e[3]));
    /* Close it: */
    p.e[4] = p.e[0];
    /* Plot the polygon: */
    hxg_paint_polygon(cvs, &p, op);
    free(p.e);
  }

void hxg_paint_sausage
  ( hxg_canvas_t *cvs,   /* Where to paint. */
    r2_vec_t *p,         /* Vertices of sausage(s). */
    double r,            /* Half-width of sausage. */
    hxg_paint_op_t *op   /* Pixel operator. */
  )
  { uint32_t np = p->ne;
    uint32_t ip;
    r2_t *pa = NULL;
    for (ip = 0; ip < np; ip++)
      { r2_t *pb = &(p->e[ip]);
        if (r2_is_finite(pb))
          { if (pa != NULL) 
              { hxg_paint_rect_stroke(cvs, pa, pb, r, op); }
            hxg_paint_circle(cvs, pb, r, op);
            pa = pb;
          }
        else
          { pa = NULL; }
      }
  }
  
