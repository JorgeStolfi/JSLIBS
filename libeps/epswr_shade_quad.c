/* See epswr_shade_quad.h */
/* Last edited on 2024-06-22 18:52:42 by stolfi */

#define _GNU_SOURCE
#include <assert.h>
#include <stdint.h>

#include <bool.h>
#include <affirm.h>

#include <epswr.h>
#include <epswr_def.h>
#include <epswr_dev.h>
#include <epswr_vis.h>

#include <epswr_shade_quad.h>

void epswr_shade_quadrilateral_subd
  ( epswr_figure_t *epsf,
    double x00, double y00, double R00, double G00, double B00,
    double x01, double y01, double R01, double G01, double B01,
    double x10, double y10, double R10, double G10, double B10,
    double x11, double y11, double R11, double G11, double B11,
    int32_t ns         /* Number of subdivisions. */
  );

void epswr_shade_quadrilateral_gouraud
  ( epswr_figure_t *epsf,
    double x00, double y00, double R00, double G00, double B00,
    double x01, double y01, double R01, double G01, double B01,
    double x10, double y10, double R10, double G10, double B10,
    double x11, double y11, double R11, double G11, double B11
  );

void epswr_shade_quadrilateral
  ( epswr_figure_t *epsf,
    double x00, double y00, double R00, double G00, double B00,
    double x01, double y01, double R01, double G01, double B01,
    double x10, double y10, double R10, double G10, double B10,
    double x11, double y11, double R11, double G11, double B11,
    int32_t ns         /* Number of subdivisions. */
  )
  { if (ns < 0)
      { /* Gouraud-shaded quadrilateral: */
        epswr_shade_quadrilateral_gouraud
          ( epsf,
            x00, y00, R00, G00, B00,
            x01, y01, R01, G01, B01,
            x10, y10, R10, G10, B10,
            x11, y11, R11, G11, B11
          );
      }
    else if (ns == 0)
      { /* Solid quadrilateral: */
        double Rt = (R00 + R01 + R10 + R11)/4;
        double Gt = (G00 + G01 + G10 + G11)/4;
        double Bt = (B00 + B01 + B10 + B11)/4;
        epswr_set_fill_color(epsf, Rt, Gt, Bt);
        epswr_quadrilateral(epsf,  x00, y00,  x01, y01,  x10, y10,  x11, y11,  TRUE, FALSE);
      }
    else
      { /* Shading by subdivision: */
        epswr_shade_quadrilateral_subd
          ( epsf,
            x00, y00, R00, G00, B00,
            x01, y01, R01, G01, B01,
            x10, y10, R10, G10, B10,
            x11, y11, R11, G11, B11,
            ns
          );
      }
  }

void epswr_shade_quadrilateral_subd
  ( epswr_figure_t *epsf,
    double x00, double y00, double R00, double G00, double B00,
    double x01, double y01, double R01, double G01, double B01,
    double x10, double y10, double R10, double G10, double B10,
    double x11, double y11, double R11, double G11, double B11,
    int32_t ns         /* Number of subdivisions. */
  )
  { 
    assert(ns > 0);
    int32_t nt = ns + 1; /* Number of intervals. */
    int32_t nv = nt + 1; /* Number of corners. */
    double x[nv], y[nv];  /* Coords of previous corner row. */
    
    int32_t i0, i1;
    /* Initialize the first row of corner coords: */
    for (i0 = 0; i0 <= nt; i0++)
      { double u0 = ((double)i0)/((double)nt);
        x[i0] = (1 - u0)*x00 + u0*x10;
        y[i0] = (1 - u0)*y00 + u0*y10;
      }
      
    /* Scan the remaining rows: */
    for (i1 = 1; i1 <= nt; i1++)
      { /* Relative position of corner row: */
        double u1 = ((double)i1)/((double)nt);
        
        /* Relative position of tile row: */
        double vt = ((double)i1-0.5)/((double)nt);
        
        /* Interpolate positions and colors along the side {x00--x01}: */
        double xt0 = (1 - u1)*x00 + u1*x01;
        double yt0 = (1 - u1)*y00 + u1*y01;
        
        double R0t = (1 - vt)*R00 + vt*R01;
        double G0t = (1 - vt)*G00 + vt*G01;
        double B0t = (1 - vt)*B00 + vt*B01;
        
        /* Interpolate positions and colors along the side {x10--x11}: */
        double xt1 = (1 - u1)*x10 + u1*x11;
        double yt1 = (1 - u1)*y10 + u1*y11;
        
        double R1t = (1 - vt)*R10 + vt*R11;
        double G1t = (1 - vt)*G10 + vt*G11;
        double B1t = (1 - vt)*B10 + vt*B11;
        
        /* Current corner on row {i1}: */
        double xt11 = xt0;
        double yt11 = yt0;
        
        /* Current corner on row {i1-1}: */
        double xt10 = x[0];
        double yt10 = y[0];
        
        for (i0 = 1; i0 <= nt; i0++)
          { /* Relative position of corner column: */
            double u0 = ((double)i0)/((double)nt);
            
            /* Relative position of tile column: */
            double ut = ((double)i0-0.5)/((double)nt);
            
            /* Update previous corner on row {i1}: */
            double xt01 = xt11;
            double yt01 = yt11;
            
            /* Update previous corner on row {i1-1}: */
            double xt00 = xt10;
            double yt00 = yt10;
            
            /* Recover current corner of previous row {i1-1}: */
            xt10 = x[i0];
            yt10 = y[i0];
        
            /* Interpolate current corner on row {i1} along the line {xt0--xt1}: */
            xt11 = (1 - u0)*xt0 + u0*xt1;
            yt11 = (1 - u0)*yt0 + u0*yt1;
            
            /* Interpolate color of tile: */
            double Rt = (1 - ut)*R0t + ut*R1t;
            double Gt = (1 - ut)*G0t + ut*G1t;
            double Bt = (1 - ut)*B0t + ut*B1t;
            
            /* Paint tile: */
            epswr_set_fill_color(epsf, Rt, Gt, Bt);
            epswr_quadrilateral(epsf, xt00, yt00, xt01, yt01, xt10, yt10, xt11, yt11, TRUE, FALSE);
            
            /* Save current corner for next row: */
            x[i0-1] = xt01;
            y[i0-1] = yt01;
          }
          
        /* Save last point: */
        x[nt] = xt11;
        y[nt] = yt11;
      }
  }

void epswr_shade_quadrilateral_gouraud
  ( epswr_figure_t *epsf,
    double x00, double y00, double R00, double G00, double B00,
    double x01, double y01, double R01, double G01, double B01,
    double x10, double y10, double R10, double G10, double B10,
    double x11, double y11, double R11, double G11, double B11
  )
  { if ((R00 < 0.0) || (R01 < 0.0) || (R10 < 0.0) || (R11 < 0.0)) { return; }
    
    /* Vertices in Postscript coordinates and clock order for visib. testing: */
    double psx[4], psy[4];
    epswr_x_to_h_coord(epsf, x00, &(psx[0])); 
    epswr_y_to_v_coord(epsf, y00, &(psy[0]));
    
    epswr_x_to_h_coord(epsf, x01, &(psx[1])); 
    epswr_y_to_v_coord(epsf, y01, &(psy[1]));
    
    epswr_x_to_h_coord(epsf, x11, &(psx[2])); 
    epswr_y_to_v_coord(epsf, y11, &(psy[2]));
    
    epswr_x_to_h_coord(epsf, x10, &(psx[3])); 
    epswr_y_to_v_coord(epsf, y10, &(psy[3]));
    
    if (! epswr_vis_polygon_is_invisible(epsf, psx, psy, 4))
      { /* Plot them: */
        FILE *wr = epsf->wr;
        /* Note that {gsquadf} wants row-by-row order, i.e [0] [1] [3] [2]: */
        fprintf(wr, "%6.1f %6.1f  %5.3f %5.3f %5.3f\n", psx[0], psy[0], R00, G00, B00);
        fprintf(wr, "%6.1f %6.1f  %5.3f %5.3f %5.3f\n", psx[1], psy[1], R01, G01, B01);
        fprintf(wr, "%6.1f %6.1f  %5.3f %5.3f %5.3f\n", psx[3], psy[3], R10, G10, B10);
        fprintf(wr, "%6.1f %6.1f  %5.3f %5.3f %5.3f",   psx[2], psy[2], R11, G11, B11);
        fprintf(wr, " gsquadf\n");
        fflush(wr);
      }
  }
