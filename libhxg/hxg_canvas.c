/* See hxg_canvas.h */
/* Last edited on 2025-01-01 03:06:31 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <r2.h>
#include <interval.h>
#include <i2.h>

#include <hxg_canvas.h>
#include <cpk_graph.h>
#include <cpk_basic.h>
  
hxg_canvas_t* hxg_canvas_new(uint32_t nx, uint32_t ny, double xStep, bool_t centered)
  { demand((nx >= 2) && (nx <= hxg_canvas_MAX_SIZE), "invalid grid dimension {nx}");
    demand((ny >= 2) && (ny <= hxg_canvas_MAX_SIZE), "invalid grid dimension {ny}");
    demand(xStep >= 1.0e-100, "invalid grid spacing {xStep}");
    hxg_canvas_t *cvs = talloc(1, hxg_canvas_t);
    cvs->size[0] = nx; cvs->size[1] = ny;
    cvs->step = (r2_t){{ xStep, xStep*sqrt(3)/2 }};
    if (centered)
      { cvs->org = (r2_t){{ -0.5*nx*cvs->step.c[0], -0.5*ny*cvs->step.c[1] }}; }
    else
      { cvs->org = (r2_t){{ 0.0, 0.0 }}; }
    return cvs;
  }

r2_t hxg_canvas_pixel_pos(hxg_canvas_t *cvs, i2_t p)
  { int32_t jx = X(p), iy = Y(p); 
    double sx = X(cvs->step), sy = Y(cvs->step);
    double x = X(cvs->org) + (jx + 0.5*(double)(iy & 1) - 0.25)*sx;
    double y = Y(cvs->org) + iy*sy;
    return (r2_t){{x,y}};
  }
  
void hxg_canvas_bbox(hxg_canvas_t *cvs, interval_t B[])
  { double xlo = cvs->org.c[0] - cvs->step.c[0]*0.25;
    double xhi = cvs->org.c[0] + cvs->step.c[0]*(cvs->size[0] - 1 + 0.25);
    
    double ylo = cvs->org.c[1];
    double yhi = cvs->org.c[1] + cvs->step.c[1]*(cvs->size[1] - 1);
    
    B[0] = (interval_t){{ xlo, xhi }};
    B[1] = (interval_t){{ ylo, yhi }};
  }

i2_vec_t hxg_canvas_half_disk_arcs(hxg_canvas_t *cvs, i2_t c, double r)
  {
    /* Get hold of params: */
    demand(r > 0, "invalid radius {r}");
    double sx = X(cvs->step), sy = Y(cvs->step);  /* Steps of canvas's grid. */
    demand((sx > 0) & (sy > 0), "invalid canvas grid steps}");
    int32_t ciy = Y(c);
    double rad2 = r*r;
    /* Compute largest displacement {eymax} in Y: */  
    uint32_t eymax = (uint32_t)floor(r/sy);
    /* Conservative estimate of template size: */
    uint32_t mD = ((uint32_t)floor((3.15/2)*r/sx)+1)*(eymax+1);
    /* Allocate the template and fill it: */
    i2_vec_t D = i2_vec_new(mD);
    uint32_t nD = 0;
    for (int32_t ey = 0; ey <= eymax; ey++)
      { int32_t iy = ciy + ey;
        double dy = ey*sy;
        /* Compute circle half-extent {hx} on this scanline: */
        double dy2 = dy*dy;
        double hx = (dy2 > rad2 ? 0.0 : sqrt(rad2 - dy2));
        /* Get index range {exmin..exmax} of pixels inside interval: */
        double fx = ((iy & 1) - (ciy & 1))/2.0; /* Rel. scanline shift */
        int32_t exmin = (int32_t)floor(-hx/sx - fx)+1;
        int32_t exmax = (int32_t)ceil(hx/sx - fx)-1;
        /* On scanline {ciy}, consider only forward pixels: */
        if (ey == 0) { exmin = 1; }
        /* Enumerate interior pixels on this scanline: */
        for (int32_t ex = exmin; ex <= exmax; ex++)
          { i2_t e = (i2_t){{ex,ey}};
            i2_vec_expand(&D, (vec_index_t)nD); /* Just in case... */
            D.e[nD] = e; nD++;
          }
      }
    i2_vec_trim(&D, nD);
    return D;
  }

double hxg_canvas_sum_over_circle
  ( hxg_canvas_t *cvs, 
    r2_t c, 
    double r, 
    double f(int32_t ix, int32_t iy, r2_t *p),
    r2_t *p
  )
  { 
    /* Get hold of params: */
    double ox = X(cvs->org), oy = Y(cvs->org);    /* Lower left corner of grid. */
    double sx = X(cvs->step), sy = Y(cvs->step);  /* Steps of canvas's grid. */
    double cx = X(c), cy = Y(c);
    double r2 = r*r;
    /* Compute range {iymin..iymax} of scanlines for {dMAx} circle at {ox,oy}. */
    int32_t iymin = (int32_t)ceil((cy - r - oy)/sy);
    int32_t iymax = (int32_t)ceil((cy + r - oy)/sy)-1;
    /* Enumerate relevant scanlines: */
    double fsum = 0;
    for (int32_t iy = iymin; iy <= iymax; iy++)
      { /* Displacement of scanline {iy} from center. */
        double dy = oy + iy*sy - cy;
        /* Find half-extent {hx} of circle on scanline {iy}: */
        double dy2 = dy*dy;
        double hx = (dy2 > r2 ? 0.0 : sqrt(r2 - dy2));
        /* Get range {ixmin..ixmax} of circle pixels on this scanline: */
        /* We assume that the circle is infinitesimally shifted to the left: */
        int32_t ixmin = (int32_t)ceil((cx - hx - ox)/sx - (iy & 1)/2.0);
        int32_t ixmax = (int32_t)ceil((cx + hx - ox)/sx - (iy & 1)/2.0) - 1;
        for (int32_t ix = ixmin; ix <= ixmax; ix++) 
          { fsum += f(ix,iy,p); }
      }
    return fsum;
  }
