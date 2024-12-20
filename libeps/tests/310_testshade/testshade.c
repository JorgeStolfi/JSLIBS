#define PROG_NAME "testshade"
#define PROG_DESC "test of {epswr.h} and {epswr_shade.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-05 10:15:14 by stolfi */

#define testshade_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

/* Created by J. Stolfi, UNICAMP sometime before 2003-09-30. */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>

#include <epswr.h>
#include <epswr_dev.h>
#include <epswr_shade_tri.h>
#include <epswr_shade_quad.h>
/* #include <epswr_color.h> */

#define FIG_XSZ (12.0)
#define FIG_YSZ (12.0)
  /* Subfigure dimensions in Client units. */

#define MIN_EPS_HMG (4.0)
#define MIN_EPS_VMG (4.0)
  /* Minimum margins for EPS figures (pt). */

int32_t main (int32_t argc, char **argv);
void SetPlotWindow(epswr_figure_t *epsf, double xsz, double ysz, int32_t ix, int32_t iy, int32_t nx, int32_t ny);
void DoEPSTests(int32_t nx, int32_t ny);
void DoPSTests(int32_t nx, int32_t ny);
void DoPaintings(epswr_figure_t *epsf, int32_t nx, int32_t ny);
void PlotShadedPoly
  ( epswr_figure_t *epsf,
    int32_t ix,
    int32_t iy,
    int32_t nx, 
    int32_t ny,
    bool_t quad,
    int32_t ns
  );
  /* Draws a shading test plot on file {epsf}.  The page
    is divided conceptually into an array of {nx × ny} sub-figures.
    The plot is done within the sub-figure in column
    {ix} and row {iy} (counting from 0 at bottom left).
    
    The test plot is a triangle if {quad} is FALSE, and
    a quadrilateral if {quad} is TRUE.  The parameter {ns}
    is the subdivision order; {ns == 0} means a single solid color,
    {ns == -1} means Gouraud shading. */ 

int32_t main (int32_t argc, char **argv)
  { DoEPSTests(2,4);
    return 0;
  }
  
void SetPlotWindow(epswr_figure_t *epsf, double xsz, double ysz, int32_t ix, int32_t iy, int32_t nx, int32_t ny)
  {
    /* Compute minimum margins in pt: */
    double min_hmg = MIN_EPS_HMG;
    double min_vmg = MIN_EPS_VMG;

    /* Get the total page/figure size in pt: */
    double tot_hsz, tot_vsz;
    epswr_get_figure_size(epsf, &tot_hsz, &tot_vsz);

    /* Compute max width and height available for all subfigures in pt: */
    double max_hsz = tot_hsz - 2*min_hmg;
    double max_vsz = tot_vsz - 2*min_vmg;

    /* Compute scale factor (pt per Client unit): */
    double hscale = max_hsz/(nx*xsz);
    double vscale = max_vsz/(ny*ysz);
    double scale = fmin(hscale, vscale);

    /* Subfigure dimensions in points: */
    double fig_hsz = scale*xsz;
    double fig_vsz = scale*ysz;

    /* Recompute true margins in points: */
    double hmg = (tot_hsz - nx*fig_hsz)/2;
    double vmg = (tot_vsz - ny*fig_vsz)/2;

    /* Compute subfigure coords in figure: */
    double fig_hlo = hmg + ix*fig_hsz;
    double fig_hhi = fig_hlo + fig_hsz;
    double fig_vlo = vmg + iy*fig_vsz;
    double fig_vhi = fig_vlo + fig_vsz;

    /* Set the window: */
    epswr_set_window(epsf, fig_hlo, fig_hhi, fig_vlo, fig_vhi, FALSE, 0, xsz, 0, ysz);
  }

void DoEPSTests(int32_t nx, int32_t ny)
  { /* Choose scale: */
    double hscale = 6.0*72/(nx*FIG_XSZ);
    double vscale = 6.0*72/(ny*FIG_XSZ);
    double scale = fmin(hscale, vscale);
    /* EPS figure size proportional to total image size. */
    double eps_hsz = ceil(nx*scale*FIG_XSZ + 2*MIN_EPS_HMG);
    double eps_vsz = ceil(ny*scale*FIG_YSZ + 2*MIN_EPS_VMG);
    char *fname = jsprintf("out/fig_%03d_%03d.eps", nx,  ny);
    FILE *wr = open_write(fname, TRUE);
    bool_t verbose = TRUE;
    epswr_figure_t *epsf = epswr_new_figure
      (wr, eps_hsz, eps_vsz, MIN_EPS_HMG, MIN_EPS_HMG, MIN_EPS_VMG, MIN_EPS_VMG, verbose);
    DoPaintings(epsf, nx, ny);
    epswr_end_figure(epsf);
    free(fname);
  }
      
double FA(double x, double y);
double FB(double x, double y);

void DoPaintings(epswr_figure_t *epsf, int32_t nx, int32_t ny)
  { 
    PlotShadedPoly(epsf, 0, 0, nx, ny, FALSE,  0);

    PlotShadedPoly(epsf, 0, 1, nx, ny, FALSE,  2);

    PlotShadedPoly(epsf, 0, 2, nx, ny, FALSE, 15);

    PlotShadedPoly(epsf, 0, 3, nx, ny, FALSE, -1);

    PlotShadedPoly(epsf, 1, 0, nx, ny, TRUE,   0);

    PlotShadedPoly(epsf, 1, 1, nx, ny, TRUE,   2);

    PlotShadedPoly(epsf, 1, 2, nx, ny, TRUE,  15);

    PlotShadedPoly(epsf, 1, 3, nx, ny, TRUE,  -1);
  }

void PlotShadedPoly
  ( epswr_figure_t *epsf,
    int32_t ix,
    int32_t iy,
    int32_t nx, 
    int32_t ny,
    bool_t quad,
    int32_t ns
  )
  { assert((ix >= 0) && (ix < nx));
    assert((iy >= 0) && (iy < ny));
    SetPlotWindow(epsf, FIG_XSZ, FIG_YSZ, ix, iy, nx, ny);
    
    double wx = FIG_XSZ;
    double wy = FIG_YSZ;
    
    /* Usable area {[0 _ wx]×[0 _ wy]} */
    
    epswr_set_pen(epsf, 0.000, 0.000, 0.000,  0.40,  0.0, 0.0);
    /* int32_t ticlo = 0; */
    /* int32_t tichi = (int32_t)floor(FIG_XY_SZ); */
    /* epswr_tics(epsf, epswr_axis_HOR, ticlo, tichi, tichi - ticlo, NULL, 1.0, 0.0); */
    /* epswr_tics(epsf, epswr_axis_VER, ticlo, tichi, tichi - ticlo, NULL, 1.0, 0.0); */
    epswr_rectangle(epsf, -0.25,+0.25+wx, -0.25,+0.25+wy, FALSE, TRUE);
    
    epswr_set_pen(epsf, 0.000, 0.000, 0.333,  0.10,  0.0, 0.0);
    
    char *cmt = jsprintf("Shading %s (ns = %d)", (char *[2]){ "triangle", "quadrilateral"}[quad], ns);
    epswr_comment(epsf, cmt);
    free(cmt);
        
    if (quad)
      { double x00 = 0.10*wx,  y00 = 0.10*wy;
        double x01 = 0.10*wx,  y01 = 0.60*wy;
        double x10 = 0.90*wx,  y10 = 0.20*wy;
        double x11 = 0.70*wx,  y11 = 0.50*wy;
        
        epswr_shade_quadrilateral
          ( epsf,
            x00, y00, 1.000, 0.300, 0.050,
            x01, y01, 0.300, 1.000, 0.050,
            x10, y10, 0.700, 0.400, 0.250,
            x11, y11, 1.000, 0.300, 1.000,
            ns
          );
            
        epswr_quadrilateral(epsf, x00, y00, x01, y01, x10, y10, x11, y11, FALSE, TRUE);
      }
    else
      { double x1 = 0.10*wx,  y1 = 0.10*wy;
        double x2 = 0.80*wx,  y2 = 0.60*wy;
        double x3 = 0.40*wx,  y3 = 0.90*wy;
        
        epswr_shade_triangle
          ( epsf,
            x1, y1, 1.000, 0.300, 0.050,
            x2, y2, 0.300, 1.000, 0.050,
            x3, y3, 0.300, 0.400, 0.250,
            ns
          );
        
        epswr_triangle(epsf, x1, y1, x2, y2, x3, y3, FALSE, TRUE);
      }
  }

