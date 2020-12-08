#define PROG_NAME "testshade"
#define PROG_DESC "test of {pswr.h} and {pswr_shade.h}"
#define PROG_VERS "1.0"

/* Last edited on 2011-06-06 17:49:38 by stolfi */

#define testshade_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

/* Created by J. Stolfi, UNICAMP sometime before 2003-09-30. */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <affirm.h>

#include <pswr.h>
#include <pswr_aux.h>
#include <pswr_shade_tri.h>
#include <pswr_shade_quad.h>
/* #include <pswr_color.h> */

#define OUT_PREFIX "out/"
  /* Prefix for output file names. */

#define FIG_XSZ (12.0)
#define FIG_YSZ (12.0)
  /* Subfigure dimensions in client units. */

#define MIN_EPS_HMG (4.0)
#define MIN_EPS_VMG (4.0)
  /* Minimum margins for EPS figures (pt). */

#define MIN_PS_HMG (72.0)
#define MIN_PS_VMG (72.0)
  /* Minimum margins for PS documents (pt). */

int main (int argc, char **argv);
void SetPlotWindow(PSStream *ps, double xsz, double ysz, int ix, int iy, int nx, int ny);
void DoEPSTests(int nx, int ny);
void DoPSTests(int nx, int ny);
void DoPaintings(PSStream *ps, int nx, int ny);
void PlotShadedPoly
  ( PSStream *ps,
    int ix,
    int iy,
    int nx, 
    int ny,
    bool_t quad,
    int ns
  );
  /* Draws a shading test plot on file {ps}.  The page
    is divided conceptually into an array of {nx × ny} sub-figures.
    The plot is done within the sub-figure in column
    {ix} and row {iy} (counting from 0 at bottom left).
    
    The test plot is a triangle if {quad} is FALSE, and
    a quadrilateral if {quad} is TRUE.  The parameter {ns}
    is the subdivision order; {ns == 0} means a single solid color,
    {ns == -1} means Gouraud shading. */ 

int main (int argc, char **argv)
  { DoEPSTests(2,4);
    DoPSTests(2,4);
    return 0;
  }
  
void SetPlotWindow(PSStream *ps, double xsz, double ysz, int ix, int iy, int nx, int ny)
  {
    /* Compute minimum margins in pt: */
    double min_hmg = (ps->eps ? MIN_EPS_HMG : MIN_PS_HMG);
    double min_vmg = (ps->eps ? MIN_EPS_VMG : MIN_PS_VMG);
    /* Get the total page/figure size in pt: */
    double tot_hsz, tot_vsz;
    pswr_get_canvas_size(ps, &tot_hsz, &tot_vsz);
    /* Compute max width and height available for all subfigures in pt: */
    double max_hsz = tot_hsz - 2*min_hmg;
    double max_vsz = tot_vsz - 2*min_vmg;
    /* Compute scale factor (pt per client unit): */
    double hscale = max_hsz/(nx*xsz);
    double vscale = max_vsz/(ny*ysz);
    double scale = fmin(hscale, vscale);
    /* Subfigure dimensions in points: */
    double fig_hsz = scale*xsz;
    double fig_vsz = scale*ysz;
    /* Recompute true margins in points: */
    double hmg = (ps->hCanvasSize - nx*fig_hsz)/2;
    double vmg = (ps->vCanvasSize - ny*fig_vsz)/2;
    /* Compute subfigure coords in canvas: */
    double fig_hlo = hmg + ix*fig_hsz;
    double fig_hhi = fig_hlo + fig_hsz;
    double fig_vlo = vmg + iy*fig_vsz;
    double fig_vhi = fig_vlo + fig_vsz;
    /* Set the window: */
    pswr_set_window(ps, 0, xsz, 0, ysz, fig_hlo, fig_hhi, fig_vlo, fig_vhi);
  }

void DoEPSTests(int nx, int ny)
  { /* Choose scale: */
    double hscale = 6.0*72/(nx*FIG_XSZ);
    double vscale = 6.0*72/(ny*FIG_XSZ);
    double scale = fmin(hscale, vscale);
    /* EPS figure size proportional to total image size. */
    double eps_hsz = ceil(nx*scale*FIG_XSZ + 2*MIN_EPS_HMG);
    double eps_vsz = ceil(ny*scale*FIG_YSZ + 2*MIN_EPS_VMG);
    PSStream *ps = pswr_new_stream(OUT_PREFIX, NULL, TRUE, "doc", NULL, FALSE, eps_hsz, eps_vsz);
    pswr_new_canvas(ps, "one");
    DoPaintings(ps, nx, ny);
    pswr_close_stream(ps);
  }
      
void DoPSTests(int nx, int ny)
  { /* EPS figure size will be set later. */
    PSStream *ps = pswr_new_stream(OUT_PREFIX, NULL, FALSE, "doc", "letter", FALSE, 0, 0);
    pswr_new_canvas(ps, "one");
    DoPaintings(ps, nx, ny);
    pswr_close_stream(ps);
  }

double FA(double x, double y);
double FB(double x, double y);

void DoPaintings(PSStream *ps, int nx, int ny)
  { 
    PlotShadedPoly(ps, 0, 0, nx, ny, FALSE,  0);

    PlotShadedPoly(ps, 0, 1, nx, ny, FALSE,  2);

    PlotShadedPoly(ps, 0, 2, nx, ny, FALSE, 15);

    PlotShadedPoly(ps, 0, 3, nx, ny, FALSE, -1);

    PlotShadedPoly(ps, 1, 0, nx, ny, TRUE,   0);

    PlotShadedPoly(ps, 1, 1, nx, ny, TRUE,   2);

    PlotShadedPoly(ps, 1, 2, nx, ny, TRUE,  15);

    PlotShadedPoly(ps, 1, 3, nx, ny, TRUE,  -1);
  }

void PlotShadedPoly
  ( PSStream *ps,
    int ix,
    int iy,
    int nx, 
    int ny,
    bool_t quad,
    int ns
  )
  { assert((ix >= 0) && (ix < nx));
    assert((iy >= 0) && (iy < ny));
    SetPlotWindow(ps, FIG_XSZ, FIG_YSZ, ix, iy, nx, ny);
    
    double wx = FIG_XSZ;
    double wy = FIG_YSZ;
    
    /* Usable area {[0 _ wx]×[0 _ wy]} */
    
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.40,  0.0, 0.0);
    /* int ticlo = 0; */
    /* int tichi = (int)floor(FIG_XY_SZ); */
    /* pswr_tics(ps, HOR, ticlo, tichi, tichi - ticlo, NULL, 1.0, 0.0); */
    /* pswr_tics(ps, VER, ticlo, tichi, tichi - ticlo, NULL, 1.0, 0.0); */
    pswr_rectangle(ps, -0.25, +0.25+wx, -0.25, +0.25+wy, FALSE, TRUE);
    
    pswr_set_pen(ps, 0.000, 0.000, 0.333,  0.10,  0.0, 0.0);
    
    char *cmt = NULL;
    asprintf(&cmt, "Shading %s (ns = %d)", (char *[2]){ "triangle", "quadrilateral"}[quad], ns);
    pswr_comment(ps, cmt);
    free(cmt);
        
    if (quad)
      { double x00 = 0.10*wx,  y00 = 0.10*wy;
        double x01 = 0.10*wx,  y01 = 0.60*wy;
        double x10 = 0.90*wx,  y10 = 0.20*wy;
        double x11 = 0.70*wx,  y11 = 0.50*wy;
        
        pswr_shade_quadrilateral
          ( ps,
            x00, y00, 1.000, 0.300, 0.050,
            x01, y01, 0.300, 1.000, 0.050,
            x10, y10, 0.700, 0.400, 0.250,
            x11, y11, 1.000, 0.300, 1.000,
            ns
          );
            
        pswr_quadrilateral(ps, x00, y00, x01, y01, x10, y10, x11, y11, FALSE, TRUE);
      }
    else
      { double x1 = 0.10*wx,  y1 = 0.10*wy;
        double x2 = 0.80*wx,  y2 = 0.60*wy;
        double x3 = 0.40*wx,  y3 = 0.90*wy;
        
        pswr_shade_triangle
          ( ps,
            x1, y1, 1.000, 0.300, 0.050,
            x2, y2, 0.300, 1.000, 0.050,
            x3, y3, 0.300, 0.400, 0.250,
            ns
          );
        
        pswr_triangle(ps, x1, y1, x2, y2, x3, y3, FALSE, TRUE);
      }
  }

