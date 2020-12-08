#define PROG_NAME "testiso"
#define PROG_DESC "test of {pswr.h} and {pswr_iso.h}"
#define PROG_VERS "1.0"

/* Last edited on 2019-04-09 19:46:12 by jstolfi */

#define testiso_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

/* Created by J. Stolfi, UNICAMP sometime before 2003-09-30. */

#define _GNU_SOURCE
#include <affirm.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <pswr.h>
#include <pswr_iso.h>
#include <pswr_color.h>

#define OUT_PREFIX "out/"
  /* Prefix for output file names. */
  
#define FIG_SIZE (140.0)
  /* Figure size in points. */

int main (int argc, char **argv);
void DoEPSTests(int mx, int my);
void DoPSTests(int mx, int my);
void DoPaintings(PSStream *ps, int mx, int my);
void PlotFunc2D
  ( PSStream *ps, 
    double func(double x, double y),
    int mx,
    int my,
    int kx,
    int ky,
    bool_t bands, 
    bool_t lines,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax,       /* Maximum isoline index. */
    double *R, double *G, double *B
  );

double FA(double x, double y);
double FB(double x, double y);
double FC(double x, double y);
double FX(double x, double y);
double FY(double x, double y);

int main (int argc, char **argv)
{ DoEPSTests(3,3);
    DoPSTests(3,3);
    return 0;
  }
  
void DoEPSTests(int mx, int my)
  {
    double hsize = FIG_SIZE*mx;
    double vsize = FIG_SIZE*my;
    PSStream *ps = pswr_new_stream(OUT_PREFIX, NULL, TRUE, "doc", NULL, FALSE, hsize+8, vsize+8);
    pswr_new_canvas(ps, NULL);
    double xsize = 12.0*mx;
    double ysize = 12.0*my;
    pswr_set_window
      ( ps, 
        0, xsize, 0, ysize,   
        4, hsize+4, 4, vsize+4
      );
    pswr_set_grid(ps, 12, 36);
    DoPaintings(ps, mx, my);
    pswr_close_stream(ps);
  }
      
void DoPSTests(int mx, int my)
  { 
    double hsize = FIG_SIZE*mx;
    double vsize = FIG_SIZE*my;
    double hmrg = (6.5*72.0 - hsize)/2;
    double vmrg = (9.0*72.0 - vsize)/2;
    PSStream *ps = pswr_new_stream(OUT_PREFIX, NULL, FALSE, "doc", "letter", FALSE, 0, 0);
    pswr_new_canvas(ps, "one");
    double xsize = 12.0*mx;
    double ysize = 12.0*my;
    pswr_set_window
      ( ps, 
        0, xsize, 0, ysize,   
        hmrg, hmrg+hsize,  vmrg, vmrg+vsize
      );
    pswr_set_grid(ps, 12, 36);
    DoPaintings(ps, mx, my);
    pswr_close_stream(ps);
  }

void DoPaintings(PSStream *ps, int mx, int my)
  { 
    double vMin = -1.0, vMax = +1.0, vStep = 0.2, vStart = vStep/2;
    
    /* Indices of first and last isoline: */
    int kMin = pswr_sup_isoline(vStart, vStep, vMin);
    int kMax = pswr_inf_isoline(vStart, vStep, vMax);
    affirm(pswr_level(vStart, vStep, kMin) >= vMin, "duh");
    affirm(pswr_level(vStart, vStep, kMin-1) < vMin, "duh");
    affirm(pswr_level(vStart, vStep, kMax) <= vMax, "duh");
    affirm(pswr_level(vStart, vStep, kMax+1) > vMax, "duh");

    /* Create standard color table (blue-white-red): */
    int N; /* Number of color bands: */
    double *Ra, *Ga, *Ba;
    pswr_make_color_table
      ( vStart, vStep, kMin, kMax, 
        0.000, 0.667, 1.000,
        1.000, 1.000, 1.000,
        1.000, 0.333, 0.000,
        &N,
        &Ra, &Ga, &Ba
      );
    affirm(N == kMax - kMin + 2, "bad N");
    
    /* Create a modified color table with some invisible colors: */
    double Rb[N], Gb[N], Bb[N];
    int i;  /* Index of color band. */
    for (i = 0; i < N; i++)
      { if ((i == 0) || (i == N-1))
          { Rb[i] = -1.000; Gb[i] = -1.000; Bb[i] = -1.000; }
        else
          { double Y = 0.299*Ra[i] + 0.587*Ga[i] + 0.114*Ba[i];
            double Rc = Ra[i]-Y, Gc = Ga[i]-Y, Bc = Ba[i]-Y; 
            double cost = (Rc - Gc)/sqrt(2);
            double sint = (2*Bc - Rc - Gc)/sqrt(6);
            double r = ((double)i)/((double)N-1);
            double rot = 2*M_PI*r;
            double cosw = cos(rot);
            double sinw = sin(rot);
            double sinu = sint*cosw + sinw*cost;
            double cosu = cost*cosw - sint*sinw;
            double Rd = + cosu/sqrt(2) - sinu/sqrt(6);
            double Gd = - cosu/sqrt(2) - sinu/sqrt(6);
            double Bd = 2*sinu/sqrt(6);
            Rb[i] = Rd+0.8*Y; Gb[i] = Gd+0.8*Y; Bb[i] = Bd+0.8*Y;
          }
      }
  
    pswr_set_pen(ps, 0.800, 0.700, 0.600,  0.10,  0.0, 0.0);
    pswr_grid_lines(ps);

    pswr_comment(ps, "Isolines only:");
    PlotFunc2D(ps, FB, mx,my,  0,0, FALSE, TRUE, vStart,vStep, kMin,kMax, Ra,Ga,Ba);

    pswr_comment(ps, "Color bands only:");
    PlotFunc2D(ps, FB, mx,my,  0,1, TRUE, FALSE, vStart,vStep, kMin,kMax, Ra,Ga,Ba);

    pswr_comment(ps, "Color bands and isolines:");
    PlotFunc2D(ps, FB, mx,my,  0,2, TRUE, TRUE,  vStart,vStep, kMin,kMax, Rb,Gb,Bb);

    pswr_comment(ps, "Isolines only:");
    PlotFunc2D(ps, FA, mx,my,  1,0, FALSE, TRUE, vStart,vStep, kMin,kMax, Ra,Ga,Ba);

    pswr_comment(ps, "Color bands with invisible colors:");
    PlotFunc2D(ps, FA, mx,my,  1,1, TRUE, FALSE, vStart,vStep, kMin,kMax, Ra,Ga,Ba);

    pswr_comment(ps, "Color bands with invisible colors and isolines:");
    PlotFunc2D(ps, FA, mx,my,  1,2, TRUE, TRUE,  vStart,vStep, kMin,kMax, Rb,Gb,Bb);

    pswr_comment(ps, "Color bands and isolines (affine function):");
    PlotFunc2D(ps, FC, mx,my,  2,0, TRUE, TRUE,  vStart,vStep, kMin,kMax, Rb,Gb,Bb);

    pswr_comment(ps, "Color bands and isolines (X coordinate):");
    PlotFunc2D(ps, FX, mx,my,  2,1, TRUE, TRUE,  vStart,vStep, kMin,kMax, Rb,Gb,Bb);

    pswr_comment(ps, "Isolines only (Y coordinate):");
    PlotFunc2D(ps, FY, mx,my,  2,2, FALSE, TRUE, vStart,vStep, kMin,kMax, Rb,Gb,Bb);
  }

void PlotFunc2D
  ( PSStream *ps,
    double func(double x, double y),
    int mx,
    int my,
    int kx,
    int ky,
    bool_t bands, 
    bool_t lines,
    double vStart,  /* Synchronize levels with this value. */
    double vStep,   /* Spacing between levels. */
    int kMin,       /* Minimum isoline index. */
    int kMax,       /* Maximum isoline index. */
    double *R, double *G, double *B
  )
  { /* Usable area {[0 _ wx]×[0 _ wy]} */
    double xc = kx*12.0 + 1.0;
    double yc = ky*12.0 + 1.0;
    double xw = 12.0 - 2*1.0;
    double yw = 12.0 - 2*1.0;
    
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
    pswr_rectangle(ps, -0.25+xc, +0.25+xw+xc, -0.25+yc, +0.25+yw+yc, FALSE, TRUE);
      
    int nx = 20;
    int ny = 20;
    
    pswr_set_pen(ps, 0.000, 0.000, 0.333,  0.10,  0.0, 0.0);
    int pass;
    for (pass = 0; pass < 2; pass++)
      { if ((pass == 0) && (! bands)) { continue; }
        if ((pass == 1) && (! lines)) { continue; }
        int ix, iy;
        for (ix = 0; ix < nx; ix++)
          { for (iy = 0; iy < ny; iy++)
              { double x0 = ((double)ix)/((double)nx);
                double x1 = ((double)ix+1)/((double)nx);
                double y0 = ((double)iy)/((double)ny);
                double y1 = ((double)iy+1)/((double)ny);
                double f00 = func(x0,y0);
                double f01 = func(x0,y1);
                double f10 = func(x1,y0);
                double f11 = func(x1,y1);
                if (bands && (pass == 0))
                  { pswr_bands_in_quadrilateral
                      ( ps,
                        xw*x0+xc, yw*y0+yc, f00,
                        xw*x0+xc, yw*y1+yc, f01, 
                        xw*x1+xc, yw*y0+yc, f10, 
                        xw*x1+xc, yw*y1+yc, f11,
                        vStart, vStep, kMin, kMax,
                        R, G, B
                      );
                  }
                if (lines && (pass == 1))
                  { pswr_isolines_in_quadrilateral
                      ( ps,
                        xw*x0+xc, yw*y0+yc, f00,
                        xw*x0+xc, yw*y1+yc, f01, 
                        xw*x1+xc, yw*y0+yc, f10, 
                        xw*x1+xc, yw*y1+yc, f11,
                        vStart, vStep, kMin, kMax
                      );
                  }
              }
          }
      }
  }

double FA(double x, double y)
  { x -= 0.5; y -= 0.5;
    double r2 = sqrt(x*x + y*y + 0.001);
    return sin(r2*4*M_PI) + 4.0*x*y;
  }

double FB(double x, double y)
  { double den = (0.1 - 0.5)*(0.1 - sqrt(0.85));
    return (x - sqrt(0.85))*(y - 0.5)/den;
  }

double FC(double x, double y)
  { return (7*x + 5*y)/12; }

double FX(double x, double y)
  { return x; }

double FY(double x, double y)
  { return y; }

