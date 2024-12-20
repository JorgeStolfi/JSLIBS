#define PROG_NAME "testfig"
#define PROG_DESC "test of {epswr.h} basic figure creation and windowing ops"
#define PROG_VERS "1.0"

/* Created by J. Stolfi, UNICAMP sometime before 2003-10-11. */
/* Last edited on 2024-12-05 10:15:08 by stolfi */

#define testfig_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>
#include <jsfile.h>
#include <jsprintf.h>

#include <epswr.h>
#include <epswr_dev.h>

int32_t main (int32_t argc, char **argv);

void DoTests(void);

void DrawPictures(epswr_figure_t *epsf, double lightDir, int32_t nh, int32_t nv);
  /* The {lightDir} is the elevation of the light source: 0.0 for equatorial, 1.0 for polar. */

void DrawPicture(epswr_figure_t *epsf, double cubeRot, double lightDir);
  /* The {cubeRot} is the rotation angle of the cube: 0.0 to 1.0 for 90 degrees. */

void DrawCubeFace
  ( epswr_figure_t *epsf, 
    int32_t ax, int32_t bx, int32_t cx, /* The normal axis and two parallel axes. */
    int32_t fc,                 /* Coordinate along axis {ax} ({+1} or {-1}). */
    double ct, double st,   /* Cosine and sine of rotation angle. */
    double cs, double ss    /* Cosine and sine of light source elevation. */
  );
  
void TestWindowOps(epswr_figure_t *epsf, int32_t kk);
  /* Tries to redefine the current window to the same Client and Device
    rectangle, using various windowing ops depending on {kk} */

void CheckWindow
  ( char *fname,
    char *coords, 
    double aMin0, double aMax0, double bMin0, double bMax0,
    double aMin1, double aMax1, double bMin1, double bMax1
  );
  /* Compares the rectangle {[aMin0 _ aMax0] × [bMin0 _ bMax0]}
    with {[aMin1 _ aMax1] × [bMin1 _ bMax1]}, and aborts if 
    not 9almost) equal. The {fname} is the name
    of the function likely to have failed, and {coords} identifies
    the coordinate system ("hv" or "xy"). */

void CheckValue(char *fname, char coord, char* side, double v0, double v1);
  /* Compares {v0} with {v1}, bombs if not (almost) equal.  The {fname} is the name
    of the function likely to have failed; {coord} is a {char} that 
    identifies the coordinate ('h', 'v', 'x', or 'y'), and {side} is "Max" or "Min. */
    
/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { DoTests();
    return 0;
  }
  
void DoTests(void)
  { 
    /* Subfigure dimensions (pt): */
    int32_t nh = 5; double hFigSize = 160.0;
    int32_t nv = 4; double vFigSize = 120.0;
    
    /* Figrue dimensions (pt) for {epswr_new_figure,epswr_new_named_figure}: */
    double hPlotSize = nh * hFigSize;
    double vPlotSize = nv * vFigSize;
    double leftMargin = 8.0;
    double rightMargin = 36.0;
    double botMargin = 72.0;
    double topMargin = 12.0;
    
    /* Temporary client window for {epswr_new_captioned_figure}: */
    double xw = hPlotSize/epswr_pt_per_mm;
    double yw = vPlotSize/epswr_pt_per_mm;
    double xmin = -xw/2, xmax = +xw/2;
    double ymin = -yw/2, ymax = +yw/2;
    int32_t capLines = 3;
    int32_t fontHeight = 10;

    int32_t nFig = 3;
    bool_t verbose = TRUE;
    for (int32_t iFig = 0;  iFig < nFig; iFig++)
      { epswr_figure_t *epsf = NULL;
        if ((iFig % 3) == 0)
          { /* Test {epswr_new_figure}, with given open file: */
            char *fileName = jsprintf("out/test_fig_%05d_A.eps", iFig);
            FILE *wr = open_write(fileName, TRUE);
            free(fileName);

            /* Create the figure object: */
            epsf = epswr_new_figure
              ( wr, hPlotSize, vPlotSize, 
                leftMargin, rightMargin, botMargin, topMargin, 
                verbose    
              );
          }
        else if ((iFig % 3) == 1)
          { /* Test {epswr_new_named_figure}: */
            epsf = epswr_new_named_figure
              ( "out", "test", "fig", iFig, "B",
                hPlotSize, vPlotSize, 
                leftMargin, rightMargin, botMargin, topMargin, 
                verbose    
              );
          }
        else if ((iFig % 3) == 2)
          { /* Test {epswr_new_captioned_figure}: */
            epsf = epswr_new_captioned_figure
              ( "out", "test", "fig", iFig, "C",
                xmin,xmax, ymin,ymax, hPlotSize, 1.2*vPlotSize,
                capLines, fontHeight,
                verbose    
              );
          }
        else
          { assert(FALSE); }
        
        /* Plot the cube in various positions: */
        double lightDir = ((double)iFig + 0.5)/((double)nFig);
        DrawPictures(epsf, lightDir, nh, nv);
        epswr_end_figure(epsf);
     }
  }
      
void DrawPictures(epswr_figure_t *epsf, double lightDir, int32_t nh, int32_t nv)
  { 
    double hMin, hMax, vMin, vMax;
    epswr_get_device_window(epsf, &hMin, &hMax, &vMin, &vMax);
    int32_t nSub = nh * nv;
    int32_t iSub = 0;
    for (int32_t iv = nv-1; iv >= 0; iv--)
      { for (int32_t ih = 0;  ih < nh; ih++)
          { /* Leave the bottom row incomplete: */
            if ((iv == 0) && (ih == (nh + 1)/2)) { return; }
            /* Set the plot window to the desired subfigure: */
            epswr_set_device_window_to_grid_cell(epsf, hMin, hMax, ih, nh, vMin, vMax, iv, nv);
            epswr_shrink_device_window(epsf, 8.0, 8.0, 18.0, 8.0);
            if (iSub < 20) { TestWindowOps(epsf, iSub); }
            if (iSub%2 == 0) 
              { epswr_set_pen(epsf, 0.000, 0.000, 0.000,  0.15,  0.0,0.0);
                epswr_frame(epsf);
              }
            double cubeRot = ((double)iSub + 0.5)/((double)nSub);
            DrawPicture(epsf, lightDir, cubeRot);
            iSub++;
          }
      }
  }
  
void TestWindowOps(epswr_figure_t *epsf, int32_t kk)
  { 
    /* Save current Device and Client windows: */
    double hMin0, hMax0, vMin0, vMax0;
    epswr_get_device_window(epsf, &hMin0, &hMax0, &vMin0, &vMax0);
    double xMin0, xMax0, yMin0, yMax0;
    epswr_get_client_window(epsf, &xMin0, &xMax0, &yMin0, &yMax0);
    
    /* Change the windows in various ways: */
    double xMin1, xMax1, yMin1, yMax1;
    double hMin1, hMax1, vMin1, vMax1;
    int32_t mod = 7;
    int32_t kkm = (kk %mod); 
    if (kkm == 0)
      { /* Test {epswr_shrink_device_window}: */
        char *fname = "epswr_shrink_device_window";
        fprintf(stderr, "--------------------\ntesting %s ...\n", fname);
        /* Shrinking amounts: */
        double dhMin = 2.0, dhMax = 5.0, dvMin = 7.0, dvMax = 9.0;
        epswr_shrink_device_window(epsf, dhMin, dhMax, dvMin, dvMax);
        epswr_get_device_window(epsf, &hMin1, &hMax1, &vMin1, &vMax1);
        CheckWindow
          ( fname, "hv", 
            hMin1, hMax1, vMin1, vMax1, 
            hMin0+dhMin, hMax0-dhMax, vMin0+dvMin, vMax0-dvMax
          );
        epswr_shrink_device_window(epsf, -dhMin, -dhMax, -dvMin, -dvMax);
        epswr_get_device_window(epsf, &hMin1, &hMax1, &vMin1, &vMax1);
        CheckWindow
          ( fname, "hv", 
            hMin1, hMax1, vMin1, vMax1, 
            hMin0, hMax0, vMin0, vMax0
          );
      }
    else if ((kkm == 1) || (kkm == 2))
      { /* Test {set_device_window}: */
        bool_t relative = ((kkm % 2) == 0);
        char *fname = (relative ? "epswr_set_device_window:rel" : "epswr_set_device_window:abs" );
        fprintf(stderr, "--------------------\ntesting %s ...\n", fname);
        /* New Device window (absolute) */
        double hMin2 = hMin0+2.0, hMax2 = hMax0-5.0, vMin2 = vMin0+7.0, vMax2 = vMax0-9.0;
        if (relative)
          { /* Set Device window relative to current one: */
            epswr_set_device_window(epsf, hMin2-hMin0, hMax2-hMin0, vMin2-vMin0, vMax2-vMin0, TRUE);
          }
        else
          { /* Set Device window absolute: */
            epswr_set_device_window(epsf, hMin2, hMax2, vMin2, vMax2, FALSE);
          }
        epswr_get_device_window(epsf, &hMin1, &hMax1, &vMin1, &vMax1);
        CheckWindow
          ( fname, "hv", 
            hMin1, hMax1, vMin1, vMax1, 
            hMin2, hMax2, vMin2, vMax2
          );
        epswr_get_client_window(epsf, &xMin1, &xMax1, &yMin1, &yMax1);
        CheckWindow
          ( fname, "xy",
            xMin1, xMax1, yMin1, yMax1, 
            0, hMax2-hMin2, 0, vMax2-vMin2
          );
      }
    else if (kkm == 3)
      { /* Test {set_client_window}: */
        char *fname = "epswr_set_client_window";
        fprintf(stderr, "--------------------\ntesting %s ...\n", fname);
        /* New Client window: */
        double xMin2 = 100, xMax2 = xMin2 + 10*(hMax0-hMin0);
        double yMin2 = 200, yMax2 = yMin2 + 10*(vMax0-vMin0);
        epswr_set_client_window(epsf, xMin2, xMax2, yMin2, yMax2);
        epswr_get_client_window(epsf, &xMin1, &xMax1, &yMin1, &yMax1);
        CheckWindow
          ( fname, "xy", 
            xMin1, xMax1, yMin1, yMax1, 
            xMin2, xMax2, yMin2, yMax2
          );
      }
    else if ((kkm == 4) || (kkm == 5))
      { /* Test {set_window}: */
        bool_t relative = ((kkm % 2) == 0);
        char *fname = (relative ? "epswr_set_window:rel" : "epswr_set_window:abs" );
        fprintf(stderr, "--------------------\ntesting %s ...\n", fname);
        /* New Device window (absolite): */
        double hMin2 = hMin0+2.0, hMax2 = hMax0-5.0, vMin2 = vMin0+7.0, vMax2 = vMax0-9.0;
        /* New Client window: */
        double xMin2 = 100, xMax2 = xMin2 + 10*(hMax2-hMin2);
        double yMin2 = 200, yMax2 = yMin2 + 10*(vMax2-vMin2);
        if (relative)
          { /* Set Device window relative to current one: */
            epswr_set_window
              ( epsf, 
                hMin2-hMin0, hMax2-hMin0, vMin2-vMin0, vMax2-vMin0, TRUE,
                xMin2, xMax2, yMin2, yMax2
              );
          }
        else
          { /* Set Device window absolute: */
            epswr_set_window
              ( epsf, 
                hMin2, hMax2, vMin2, vMax2, FALSE,
                xMin2, xMax2, yMin2, yMax2
              );
          }
        epswr_get_device_window(epsf, &hMin1, &hMax1, &vMin1, &vMax1);
        CheckWindow
          ( fname, "hv", 
            hMin1, hMax1, vMin1, vMax1, 
            hMin2, hMax2, vMin2, vMax2
          ); 
        epswr_get_client_window(epsf, &xMin1, &xMax1, &yMin1, &yMax1);
        CheckWindow
          ( fname, "xy", 
            xMin1, xMax1, yMin1, yMax1,
            xMin2, xMax2, yMin2, yMax2
          );
      }
    else
      { /* Do nothing: */
      }
    /* Restore the original windows: */
    epswr_set_device_window(epsf, hMin0, hMax0, vMin0, vMax0, FALSE);
    epswr_set_client_window(epsf, xMin0, xMax0, yMin0, yMax0);
  }
  
void CheckWindow
  ( char *fname,
    char *coords, 
    double aMin0, double aMax0, double bMin0, double bMax0,
    double aMin1, double aMax1, double bMin1, double bMax1
  )
  {
    CheckValue(fname, coords[0], "Min", aMin0, aMin1);
    CheckValue(fname, coords[0], "Max", aMax0, aMax1);
    CheckValue(fname, coords[1], "Min", bMin0, bMin1);
    CheckValue(fname, coords[1], "Max", bMax0, bMax1);
  }
  
void CheckValue(char *fname, char coord, char* side, double v0, double v1)
  {
    if (fabs(v0 - v1) > 1.0e-6)
      { fprintf(stderr, "** error in {%s}: %c%s", fname, coord, side);
        fprintf(stderr, " is %.8f, should be %.8f\n", v0, v1);
        affirm(FALSE, "test failed");
      }
  }

void DrawPicture(epswr_figure_t *epsf, double lightDir, double cubeRot)
  {
    /* Get the Device window: */
    double hMin, hMax, vMin, vMax;
    epswr_get_device_window(epsf, &hMin, &hMax, &vMin, &vMax);
    
    /* Paint and draw the cube: */
    double t = cubeRot*0.5*M_PI; /* Cube rotation angle in radians. */
    double ct = cos(t), st = sin(t);
    
    double s = lightDir*0.5*M_PI; /* Light source elevation angle in radians. */
    double cs = cos(s), ss = sin(s);
    
    double R = sqrt(3.0) + 0.2;
    epswr_set_client_window(epsf, -R, +R, -R, +R);
    
    /* Enumerate the eight faces of the cube: */
    for (int32_t ax = 0;  ax < 3; ax++)
      { int32_t bx = (ax + 1) % 3, cx = (ax + 2) % 3;
        for (int32_t fc = -1; fc <= +1; fc += 2)
          { DrawCubeFace(epsf, ax, bx, cx, fc, ct, st, cs, ss); }
      }

    /* Write the caption: */
    double capFontSize = 12.0;
    epswr_set_text_font(epsf, "Courier", capFontSize);
    double vTopText = vMin - 4.0;
    double vBotText = vTopText - 2*capFontSize;
    epswr_dev_set_text_geometry(epsf, hMin, hMax, vBotText, vTopText, 0.0);
    char *xRot = jsprintf("# %5.3f #", cubeRot);
    epswr_set_fill_color(epsf, 0.000,0.000,1.000);
    epswr_text(epsf, xRot, FALSE, 0.5, TRUE, FALSE);
    free(xRot);
  }
            
void DrawCubeFace
  ( epswr_figure_t *epsf, 
    int32_t ax, int32_t bx, int32_t cx, /* The normal axis and two parallel axes. */
    int32_t fc,                 /* Coordinate along axis {ax} ({+1} or {-1}). */
    double ct, double st,   /* Cosine and sine of rotation angle. */
    double cs, double ss    /* Cosine and sine of light source elevation. */
  )
  {
    auto void norm(double u[]);            /* Normalizes {u[0..2]} to unit length.*/
    auto double dot(double u[], double v[]); /* Inner product of {u[0..2]} and {v[0..2]}.*/
    
    void norm(double u[])
      { double s2 = 0.0;
        for (int32_t k = 0;  k < 3; k++) { s2 += u[k]*u[k]; }
        double m = sqrt(s2);
        for (int32_t k = 0;  k < 3; k++) { u[k] /= m; }
      }
  
    double dot(double u[], double v[])
      { double s = 0.0;
        for (int32_t k = 0;  k < 3; k++) { s += u[k]*v[k]; }
        return s;
      }
  
    /* Camera and lighting direction vectors: */
    double obs[3] = {  +3.0,  +2.0,  +1.0 }; norm(obs); /* Observer. */
    double xpv[3] = {  -2.0,  +3.0,  00.0 }; norm(xpv); /* Horizontal camera axis. */
    double ypv[3] = {  -3.0,  -2.0, +13.0 }; norm(ypv); /* Vertical camera axis. */
    double ldv[3] = { +5.0*cs, +12.0*cs, +13*ss }; norm(ldv); /* Light source. */
    
    /* Object points: */
    double fn[3]; /* Face normal. */
    double xp[4], yp[4];  /* Plot coordinates of vertices. */
    
    /* Clear face normal: */
    for (int32_t ar = 0;  ar < 3; ar++) { fn[ar] = 0.0; }
    
    /* Enumerate vertices of face {v[ax] == fc} in cyclic order, compute face normal: */
    double v[3];
    int32_t r = 0;
    for (int32_t j = -1; j <= +1; j += 2)
      { for (int32_t k = -1; k <= +1; k += 2)
          { v[ax] = fc;
            v[bx] = j;
            v[cx] = k*j; /* Hack to get the right order. */
            
            /* Rotate point {v} by angle {arg(ct,st)} around all axes: */
            for (int32_t ar = 0;  ar < 3; ar++)
              { int32_t br = (ar + 1) % 3, cr = (ar + 2) % 3;
                double xt =  ct*v[br] + st*v[cr];
                double yt = -st*v[br] + ct*v[cr];
                v[br] = xt; v[cr] = yt;
              }
            
            /* Accumulate into face normal: */
            for (int32_t ar = 0;  ar < 3; ar++) { fn[ar] += v[ar]; }
            
            /* Project and store in {xp[r],yp[r]}: */
            
            xp[r] = dot(v, xpv);
            yp[r] = dot(v, ypv);
            r++;
          }
      }
    norm(fn);
    
    /* Check visibility and plot: */
    double cv = dot(fn, obs);
    if (cv > 0.0)
      { /* Compute illumination: */
        double cg = dot(fn, ldv);
        double shade = 0.5 + 0.5*(cg <= 0.0 ? 0.0 : cg);
        epswr_set_pen(epsf, 0.000, 0.000, 0.000,  0.10,  0.0,0.0);
        epswr_set_fill_color(epsf, shade*1.000, shade*0.850, shade*0.500);
        epswr_polygon(epsf, TRUE, xp, yp, 4, TRUE, TRUE, FALSE);
        double xpfn = dot(fn, xpv);
        double ypfn = dot(fn, ypv);
        epswr_set_pen(epsf, shade*1.000, shade*0.000, shade*0.000,  0.20,  0.0,0.0);
        epswr_segment(epsf, xpfn, ypfn, 2*xpfn, 2*ypfn);
      }
  }
