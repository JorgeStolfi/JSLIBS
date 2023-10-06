/* See epswr.h */
/* Last edited on 2023-10-01 20:03:02 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsstring.h>
#include <jsfile.h>
#include <jstime.h>
#include <rn.h>

#include <epswr.h>
#include <epswr_def.h>
#include <epswr_dev.h>
#include <epswr_font_list.h>
#include <epswr_vis.h>

/* INTERNAL PROTOTYPES */

#define INF INFINITY
  /* IEEE plus infinity. */

void epswr_write_file_preamble(epswr_figure_t *eps);
  /*   */

void epswr_write_file_postamble(epswr_figure_t *eps);
  /* Writes the file's postamble, including the fonts list. */

void epswr_font_list_initialize(int32_t *nfontsP, char ***fontsP);
  /* Initializes the font list with the "Courier" font.
    The list must be {NULL}. */

void epswr_font_list_free(int32_t *nfontsP, char ***fontsP);
  /* Frees all storage used by the font list and sets it to NULL. */
    
void epswr_check_param
  ( const char *name, 
    double z, 
    double zMin, 
    double zMax
  );
  /* Checks whether {z} is in the range {[zMin __ zMax]}; 
    aborts with error if not. */

/* START OF NEW FIGURE */

epswr_figure_t *epswr_new_figure
  ( FILE *wr,
    double hPlotSize,    /* Plot window width (in pt). */
    double vPlotSize,    /* Plot window height (in pt). */
    double leftMargin,   /* Extra margin at left (pt). */
    double rightMargin,  /* Extra margin at right (pt). */
    double botMargin,    /* Extra margin at bottom (in pt). */
    double topMargin,    /* Extra margin at top (in pt). */
    bool_t verbose       /* TRUE to print diagnostics. */
  )
  {
    double hSize = ceil(leftMargin + hPlotSize + rightMargin);
    double vSize = ceil(botMargin + vPlotSize + topMargin);
    
    epswr_figure_t *eps = epswr_dev_new_figure(wr, hSize, vSize, verbose);
    
    /* Set the initial plot window and reset the grid: */
    double hMin = leftMargin;
    double hMax = hMin + hPlotSize;
    double vMin = botMargin;
    double vMax = vMin + vPlotSize;
    epswr_dev_set_window(eps, hMin, hMax, vMin, vMax, FALSE); 
    epswr_def_reset_client_window(eps);

    /* Set the text window to the plot window: */
    epswr_dev_set_text_geometry(eps, hMin, hMax, vMin, vMax, 0.0); 
    
    fflush(eps->wr);
 
    return eps;
  }
  
epswr_figure_t *epswr_new_named_figure
  ( char *dir, 
    char *prefix,
    char *name,
    int32_t seq, 
    char *suffix,
    double hPlotSize,    /* Initial plot window width (in pt). */
    double vPlotSize,    /* Initial plot window height (in pt). */
    double leftMargin,   /* Extra margin at left (pt). */
    double rightMargin,  /* Extra margin at right (pt). */
    double botMargin,    /* Extra margin at bottom (in pt). */
    double topMargin,    /* Extra margin at top (in pt). */
    bool_t verbose       /* TRUE to print diagnostics. */
  )
  {
    if (dir == NULL) { dir = ""; }
    if (prefix == NULL) { prefix = ""; }
    if (name == NULL) { name = ""; }
    if (suffix == NULL) { suffix = ""; }
    FILE *wr;
    if (strlen(dir) + strlen(prefix) + strlen(name) + strlen(suffix) == 0)
      { /* All name parts are omitted: */
        fprintf(stderr, "writing EPS figure to {stdout}\n");
        wr = stdout;
      }
    else
      { char *dir_s = (dir[0] == 0 ? "" : "/");
        char *prefix_u = ((prefix[0] == 0) || (name[0] == 0) ? "" : "_");
        char *fname = NULL;
        if (seq >= 0) 
          { char *name_u = ((prefix[0] == 0) && (name[0] == 0) ? "" : "_");
            char *seq_u = (suffix[0] == 0 ? "" : "_");
            asprintf(&fname, "%s%s%s%s%s%s%05d%s%s.eps", dir, dir_s, prefix, prefix_u, name, name_u, seq, seq_u, suffix); }
        else
          { char *name_u = (((prefix[0] == 0) && (name[0] == 0)) || (suffix[0] == 0) ? "" : "_");
            asprintf(&fname, "%s%s%s%s%s%s%s.eps", dir, dir_s, prefix, prefix_u, name, name_u, suffix);
          }
        wr = open_write(fname, TRUE);
        free(fname);
      }
    epswr_figure_t *eps = epswr_new_figure
      ( wr, hPlotSize, vPlotSize,
        leftMargin, rightMargin, botMargin, topMargin,
        verbose
      );
    return eps;
  }
  
epswr_figure_t *epswr_new_captioned_figure
  ( char *dir,
    char *prefix,
    char *name,
    int32_t seq, 
    char *suffix,
    double xMin,         
    double xMax,         
    double yMin,
    double yMax,
    double hMaxSize,     /* Max plot window width (in pt). */
    double vMaxSize,     /* Max plot window height (in pt). */
    int32_t capLines,    /* Number of lines to reserve for caption. */
    double fontHeight,   /* Font height for caption text. */
    bool_t verbose       /* TRUE to print diagnostics. */
  )
  {
    double mrg = 4.0; /* Margin width in {pt}. */
    double mrg_bot = mrg + (capLines == 0 ? 0 : capLines*fontHeight + mrg);
    
    /* Determine the plot area size {hSize,vSize} (pt): */
    double wx = fabs(xMax-xMin);
    double wy = fabs(yMax-yMin);
    demand((wx > 1.0e-200) && (wy > 1.0e-200), "client window too small");
    double scale = +INF;
    if (isfinite(hMaxSize))
      { demand(hMaxSize > 1.0e-200, "invalid {hMaxSize}");
        scale = fmin(scale, hMaxSize/wx);
      }
    if (isfinite(vMaxSize))
      { demand(vMaxSize > 1.0e-200, "invalid {vMaxSize}");
        scale = fmin(scale, vMaxSize/wy);
      }
    demand(isfinite(scale), "must specify at least one {hMaxSize,vMaxSize}");
    double hSize = scale*wx;
    double vSize = scale*wy;
    /* Sanity check: */
    demand((hSize <= epswr_MAX_SIZE) && (vSize <= epswr_MAX_SIZE), "figure too big");
    
    /* Create the figure:  */
    epswr_figure_t *eps = epswr_new_named_figure
      ( dir, prefix, name, seq, suffix,
        hSize, vSize, mrg, mrg, mrg_bot, mrg,
        verbose
      ); 
      
    /* Client window: */
    epswr_set_client_window(eps, xMin,xMax, yMin,yMax);
    
    /* Caption area: */
    epswr_set_text_geometry(eps, FALSE, 0,hSize, mrg-mrg_bot, -mrg, 0.0);
    epswr_set_text_font(eps, "Courier", fontHeight);
    epswr_set_fill_color(eps, 0,0,0);
    return eps;
  }

void epswr_end_figure(epswr_figure_t *eps)
  { epswr_dev_end_figure(eps); }

void epswr_get_figure_size
  ( epswr_figure_t *eps,  /* Picture stream. */
    double *hSizeP,  /* OUT: Width of figure (in pt). */
    double *vSizeP   /* OUT: Height of figure (in pt). */
  )  
  { epswr_dev_get_figure_size(eps, hSizeP, vSizeP); }

/* PLOTTING WINDOW */

void epswr_set_device_window
  ( epswr_figure_t *eps,
    double hMin, double hMax,
    double vMin, double vMax,
    bool_t relative
  )
  { epswr_dev_set_window(eps, hMin, hMax, vMin, vMax, relative);
    epswr_def_reset_client_window(eps);
  }

void epswr_get_device_window
  ( epswr_figure_t *eps,
    double *hMinP, double *hMaxP,
    double *vMinP, double *vMaxP
  )
  { epswr_dev_get_window(eps, hMinP, hMaxP, vMinP, vMaxP); }

void epswr_shrink_device_window
  ( epswr_figure_t *eps, 
    double dhMin, double dhMax, 
    double dvMin, double dvMax
  )
  { epswr_dev_shrink_window(eps, dhMin, dhMax, dvMin, dvMax);
    epswr_def_reset_client_window(eps);
  }
    
void epswr_set_device_window_to_grid_cell
  ( epswr_figure_t *eps, 
    double hMin, double hMax, int32_t col, int32_t cols,
    double vMin, double vMax, int32_t row, int32_t rows
  )
  { epswr_dev_set_window_to_grid_cell(eps, hMin, hMax, col, cols, vMin, vMax, row, rows);
    epswr_def_reset_client_window(eps);
  }

void epswr_set_client_window
  ( epswr_figure_t *eps,
    double xMin, double xMax,
    double yMin, double yMax
  )
  { 
    if (eps->verbose)  
      { epswr_def_report_window(eps, "epswr_set_client_window", "Client", "requested", xMin, xMax, xMin, xMax); }

    demand(xMin < xMax, "Client X scale must be positive");
    demand(yMin < yMax, "Client Y scale must be positive");
    
    /* Shrink the plot window to ensure equal scales on both axes: */
    double hMin = eps->hMin, hMax = eps->hMax, hCtr = (hMax + hMin)/2;
    double vMin = eps->vMin, vMax = eps->vMax, vCtr = (vMax + vMin)/2;
    
    double hScale = (hMax - hMin)/(xMax - xMin); /* Device units per Client unit. */
    double vScale = (vMax - vMin)/(yMax - yMin); /* Device units per Client unit. */
    
    if (fabs(hScale) < fabs(vScale))
      { /* Reduce {fabs(vScale)} by shrinking {vMax,vMin}: */
        double vShrink = fabs(hScale)/fabs(vScale);
        vMin = vCtr + vShrink*(vMin - vCtr);
        vMax = vCtr + vShrink*(vMax - vCtr);
      }
    else if (fabs(vScale) < fabs(hScale))
      { /* Reduce {fabs(hScale)} by shrinking {hMax,hMin}: */
        double hShrink = fabs(vScale)/fabs(hScale);
        hMin = hCtr + hShrink*(hMin - hCtr);
        hMax = hCtr + hShrink*(hMax - hCtr);
      }

    if ((hMin != eps->hMin) || (hMax != eps->hMax) || (vMin != eps->vMin) || (vMax != eps->vMax))
      { epswr_dev_set_window(eps, hMin, hMax, vMin, vMax, FALSE); }

    /* Paranoia: check for equal scales, allowing for some rounding errors: */
    hScale = (hMax - hMin)/(xMax - xMin);
    vScale = (vMax - vMin)/(yMax - yMin);
    { double scalesum = fabs(hScale) + fabs(vScale);
      double scalediff = fabs(fabs(hScale) - fabs(vScale));
      demand(scalediff/scalesum < 0.0001, "unequal scales");
    }

    /* Save the window parameters in the {epswr_figure_t} record: */
    epswr_def_set_client_window(eps, xMin, xMax, yMin, yMax, "epswr_set_client_window");
    /* Write the window setup commands to the Postcript file: */
    FILE *wr = eps->wr;

    epswr_dev_write_window_set_cmds(wr, hMin, hMax, vMin, vMax);
  }

void epswr_get_client_window
  ( epswr_figure_t *eps,
    double *xMinP, double *xMaxP,
    double *yMinP, double *yMaxP
  )
  { 
    epswr_def_get_client_window(eps, xMinP, xMaxP, yMinP, yMaxP);
  }
  
void epswr_set_window
  ( epswr_figure_t *eps,
    
    double hMin, double hMax,
    double vMin, double vMax,
    bool_t relative,

    double xMin, double xMax,
    double yMin, double yMax
  )
  {
    epswr_dev_set_window(eps, hMin, hMax, vMin, vMax, relative);
    epswr_set_client_window(eps, xMin, xMax, yMin, yMax);
  }

/* DRAWING COMMANDS */

void epswr_set_pen
  ( epswr_figure_t *eps,
    double R, double G, double B,
    double width,
    double dashLength,
    double dashSpace
  )
  { double pswidth = width * epswr_pt_per_mm;
    double psdashLength = dashLength * epswr_pt_per_mm;
    double psdashSpace = dashSpace * epswr_pt_per_mm;
    epswr_dev_set_pen(eps, R, G, B, pswidth, psdashLength, psdashSpace);
  }

void epswr_segment
  ( epswr_figure_t *eps,
    double xa, double ya,
    double xb, double yb
  )
  { double psxa; epswr_x_to_h_coord(eps, xa, &(psxa));
    double psya; epswr_y_to_v_coord(eps, ya, &(psya));
    double psxb; epswr_x_to_h_coord(eps, xb, &(psxb));
    double psyb; epswr_y_to_v_coord(eps, yb, &(psyb));
    epswr_dev_segment(eps, psxa, psya, psxb, psyb);
  }

void epswr_curve
  ( epswr_figure_t *eps,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double xd, double yd
  )
  { double psxa; epswr_x_to_h_coord(eps, xa, &(psxa));
    double psya; epswr_y_to_v_coord(eps, ya, &(psya));
    double psxb; epswr_x_to_h_coord(eps, xb, &(psxb));
    double psyb; epswr_y_to_v_coord(eps, yb, &(psyb));
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psxd; epswr_x_to_h_coord(eps, xd, &(psxd));
    double psyd; epswr_y_to_v_coord(eps, yd, &(psyd));
    epswr_dev_curve(eps, psxa, psya, psxb, psyb, psxc, psyc, psxd, psyd);
  }

void epswr_coord_line(epswr_figure_t *eps, epswr_axis_t axis, double pos)
  { double pspos;
    if (axis == epswr_axis_HOR)
      { epswr_x_to_h_coord(eps, pos, &(pspos)); }
    else if (axis == epswr_axis_VER)
      { epswr_y_to_v_coord(eps, pos, &(pspos)); }
    else
      { affirm(FALSE, "invalid axis"); pspos = 0.0; }
    epswr_dev_coord_line(eps, axis, pspos);
  }

void epswr_coord_lines(epswr_figure_t *eps, epswr_axis_t axis, double start, double step)
  { double psstart, psstep;
    if (axis == epswr_axis_HOR)
      { epswr_x_to_h_coord(eps, start, &(psstart));
        epswr_x_to_h_dist(eps, step, &(psstep));
      }
    else if (axis == epswr_axis_VER)
      { epswr_y_to_v_coord(eps, start, &(psstart));
        epswr_y_to_v_dist(eps, step, &(psstep));
      }
    else
      { affirm(FALSE, "invalid axis"); 
        psstart = 0.0; psstep = 0.0;
      }
    epswr_dev_coord_lines(eps, axis, psstart, psstep);
  }

void epswr_axis(epswr_figure_t *eps, epswr_axis_t axis, double pos, double lo, double hi)
  { 
    if (axis == epswr_axis_HOR)
      { epswr_segment(eps, lo, pos, hi, pos); }
    else if (axis == epswr_axis_VER)
      { epswr_segment(eps, pos, lo, pos, hi); }
    else
      { affirm(FALSE, "invalid axis"); }
  }

void epswr_frame(epswr_figure_t *eps)
  { epswr_dev_frame(eps); }

/* CLOSED FIGURES */

void epswr_set_fill_color(epswr_figure_t *eps, double R, double G, double B)
  { epswr_dev_set_fill_color(eps, R, G, B); }

void epswr_rectangle
  ( epswr_figure_t *eps,
    double xlo, double xhi,
    double ylo, double yhi,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxlo; epswr_x_to_h_coord(eps, xlo, &(psxlo));
    double psxhi; epswr_x_to_h_coord(eps, xhi, &(psxhi));
    double psylo; epswr_y_to_v_coord(eps, ylo, &(psylo));
    double psyhi; epswr_y_to_v_coord(eps, yhi, &(psyhi));
    epswr_dev_rectangle(eps, psxlo, psxhi, psylo, psyhi, fill, draw);
  }

void epswr_triangle
  ( epswr_figure_t *eps,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxa; epswr_x_to_h_coord(eps, xa, &(psxa));
    double psya; epswr_y_to_v_coord(eps, ya, &(psya));
    double psxb; epswr_x_to_h_coord(eps, xb, &(psxb));
    double psyb; epswr_y_to_v_coord(eps, yb, &(psyb));
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    epswr_dev_triangle(eps, psxa, psya, psxb, psyb, psxc, psyc, fill, draw);
  }

void epswr_quadrilateral
  ( epswr_figure_t *eps,
    double x00, double y00,
    double x01, double y01,
    double x10, double y10,
    double x11, double y11,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double x[4], y[4];
    x[0] = x00; y[0] = y00;
    x[1] = x01; y[1] = y01;
    x[2] = x11; y[2] = y11;
    x[3] = x10; y[3] = y10;
    epswr_polygon(eps, TRUE, x, y, 4, fill, draw, TRUE);
  }

void epswr_parallelogram
  ( epswr_figure_t *eps,
    double xc, double yc,
    double xu, double yu,
    double xv, double yv,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double x[4], y[4];
    x[0] = xc - xu - xv; y[0] = yc - yu - yv;
    x[1] = xc - xu + xv; y[1] = yc - yu + yv;
    x[2] = xc + xu + xv; y[2] = yc + yu + yv;
    x[3] = xc + xu - xv; y[3] = yc + yu - yv;
    epswr_polygon(eps, TRUE, x, y, 4, fill, draw, TRUE);
  }

void epswr_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double x[], double y[],
    int32_t n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double *psx = rn_alloc(n);
    double *psy = rn_alloc(n);
    /* Map points to Device coordinates: */
    for (int32_t i=0; i<n; i++)
      { epswr_x_to_h_coord(eps, x[i], &(psx[i]));
        epswr_y_to_v_coord(eps, y[i], &(psy[i]));
      }
    epswr_dev_polygon(eps, closed, psx, psy, n, fill, draw, evenOdd);
    free(psx); free(psy);
  }

void epswr_rounded_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double x[], double y[],
    int32_t n,
    double rad,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    /* Corners: */
    double *psx = rn_alloc(n);
    double *psy = rn_alloc(n);
    /* Map points and radius to Device coordinates: */
    for (int32_t i=0; i<n; i++)
      { epswr_x_to_h_coord(eps, x[i], &(psx[i]));
        epswr_y_to_v_coord(eps, y[i], &(psy[i]));
      }
    double psrad; epswr_y_to_v_dist(eps, rad, &(psrad)); /* Assumes equal scales. */
    epswr_dev_rounded_polygon(eps, closed, psx, psy, n, psrad, fill, draw, evenOdd);
    free(psx); free(psy);
  }

void epswr_bezier_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double x[], double y[],
    int32_t n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    int32_t np = 4*n; /* Number of points. */
    double *psx = rn_alloc(np);
    double *psy = rn_alloc(np);
    /* Map points to Device coordinates: */
    for (int32_t i = 0; i < np; i++)
      { epswr_x_to_h_coord(eps, x[i], &(psx[i]));
        epswr_y_to_v_coord(eps, y[i], &(psy[i]));
      }
    epswr_dev_bezier_polygon(eps, closed, psx, psy, n, fill, draw, evenOdd);
    free(psx); free(psy);
  }

void epswr_circle
  ( epswr_figure_t *eps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psrad; epswr_y_to_v_dist(eps, rad, &(psrad));
    epswr_dev_circle(eps, psxc, psyc, psrad, fill, draw);
  }

void epswr_lune
  ( epswr_figure_t *eps,
    double xc, double yc, double rad, double tilt,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psrad; epswr_y_to_v_dist(eps, rad, &(psrad));
    /* {epswr_dev_lune} wants the tilt in radians: */
    double pstilt = tilt * 180.0 / M_PI;
    epswr_dev_lune(eps, psxc, psyc, psrad, pstilt, fill, draw);
  }

void epswr_slice
  ( epswr_figure_t *eps,
    double xc, double yc, double rad, 
    double start, double stop,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psrad; epswr_y_to_v_dist(eps, rad, &(psrad));
    /* {epswr_dev_slice} wants the angles in degrees too: */
    epswr_dev_slice(eps, psxc, psyc, psrad, start, stop, fill, draw);
  }

/* PLOT MARKS */

void epswr_dot
  ( epswr_figure_t *eps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psrad = rad * epswr_pt_per_mm;
    epswr_dev_dot(eps, psxc, psyc, psrad, fill, draw);
  }

void epswr_tic
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double xc, double yc, 
    double ticSize,
    double align 
  )
  { double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psticSize = ticSize*epswr_pt_per_mm;
    epswr_dev_tic(eps, axis, psxc, psyc, psticSize, align);
  }   
 
void epswr_cross
  ( epswr_figure_t *eps, 
    double xc, double yc, double rad, bool_t diag,
    bool_t draw
  )
  { if (! draw) { return; }
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psrad = rad * epswr_pt_per_mm;
    epswr_dev_cross(eps, psxc, psyc, psrad, diag, draw);
  }

void epswr_asterisk
  ( epswr_figure_t *eps, 
    double xc, double yc, double rad,
    bool_t draw
  )
  { if (!draw) { return; }
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psrad = rad * epswr_pt_per_mm;
    epswr_dev_asterisk(eps, psxc, psyc, psrad, draw);
  }

void epswr_square
  ( epswr_figure_t *eps,
    double xc, double yc, double rad,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psrad = rad * epswr_pt_per_mm;
    epswr_dev_square(eps, psxc, psyc, psrad, fill, draw);
  }

void epswr_diamond
  ( epswr_figure_t *eps, 
    double xc, double yc,
    double xRad, double yRad,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxc; epswr_x_to_h_coord(eps, xc, &(psxc));
    double psyc; epswr_y_to_v_coord(eps, yc, &(psyc));
    double psxRad = xRad * epswr_pt_per_mm;
    double psyRad = yRad * epswr_pt_per_mm;
    epswr_dev_diamond(eps, psxc, psyc, psxRad, psyRad, fill, draw);
  }   

void epswr_arrowhead 
  ( epswr_figure_t *eps,
    double xa, double ya, double xb, double yb,
    double width, double length, 
    double fraction,
    bool_t fill, bool_t draw
  )
  {
    if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }

    double psxa; epswr_x_to_h_coord(eps, xa, &(psxa));
    double psya; epswr_y_to_v_coord(eps, ya, &(psya));
    double psxb; epswr_x_to_h_coord(eps, xb, &(psxb));
    double psyb; epswr_y_to_v_coord(eps, yb, &(psyb));
    
    double pswidth = width*epswr_pt_per_mm;
    double pslength = length * epswr_pt_per_mm;
    epswr_dev_arrowhead(eps, psxa, psya, psxb, psyb, pswidth, pslength, fraction, fill, draw);
  } 

void epswr_grid_lines(epswr_figure_t *eps, int32_t cols, int32_t rows)
  { epswr_dev_grid_lines(eps, cols, rows); }

void epswr_grid_cell
  ( epswr_figure_t *eps, 
    int32_t col, int32_t cols,
    int32_t row, int32_t rows,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    epswr_dev_grid_cell(eps, col, cols, row, rows, fill, draw);
  }
    
void epswr_set_label_font(epswr_figure_t *eps, const char *font, double size)
  { epswr_dev_set_label_font(eps, font, size); }

void epswr_label
  ( epswr_figure_t *eps, 
    const char *text,
    const char *strut,
    double x, double y, 
    double rot, 
    bool_t clipped,
    double hAlign, double vAlign,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psx; epswr_x_to_h_coord(eps, x, &(psx));
    double psy; epswr_y_to_v_coord(eps, y, &(psy));
    epswr_dev_label(eps, text, strut, psx, psy, rot, clipped, hAlign, vAlign, fill, draw);
  }
  
void epswr_set_text_geometry
  ( epswr_figure_t *eps, 
    bool_t client,
    double xhMin, double xhMax, 
    double yvMin, double yvMax,
    double rot
  )
  {
    double hMin, hMax, vMin, vMax;
    if (client)
      { epswr_x_to_h_coord(eps, xhMin, &(hMin));
        epswr_x_to_h_coord(eps, xhMax, &(hMax));
        epswr_y_to_v_coord(eps, yvMin, &(vMin));
        epswr_y_to_v_coord(eps, yvMax, &(vMax));
      }
    else
      { hMin = xhMin + eps->hMin; 
        hMax = xhMax + eps->hMin; 
        vMin = yvMin + eps->vMin; 
        vMax = yvMax + eps->vMin;
      }
    epswr_dev_set_text_geometry(eps, hMin, hMax, vMin, vMax, rot);
  }

void epswr_set_text_font(epswr_figure_t *eps, const char *font, double size)
  { epswr_dev_set_text_font(eps, font, size); }

void epswr_text
  ( epswr_figure_t *eps, 
    const char *text, 
    bool_t clipped,
    double hAlign,
    bool_t fill, bool_t draw
  )
  {
    if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    epswr_dev_text(eps, text, clipped, hAlign, fill, draw);
  }
 
/* CLIENT TO DEVICE COORDINATE CONVERSION */ 

void epswr_x_to_h_coord(epswr_figure_t *eps, double x, double *hP)
  { double scale = (eps->hMax - eps->hMin)/(eps->xMax - eps->xMin);
    (*hP) = eps->hMin + scale * (x - eps->xMin);
  }

void epswr_h_to_x_coord(epswr_figure_t *eps, double h, double *xP)
  { double scale = (eps->xMax - eps->xMin)/(eps->hMax - eps->hMin);
    (*xP) = eps->xMin + scale * (h - eps->hMin);
  }
 
void epswr_y_to_v_coord(epswr_figure_t *eps, double y, double *vP)
  { double scale = (eps->vMax - eps->vMin)/(eps->yMax - eps->yMin);
    (*vP) = eps->vMin + scale * (y - eps->yMin);
  }

void epswr_v_to_y_coord(epswr_figure_t *eps, double v, double *yP)
  { double scale = (eps->yMax - eps->yMin)/(eps->vMax - eps->vMin);
    (*yP) = eps->yMin + scale * (v - eps->vMin);
  }

void epswr_x_to_h_dist(epswr_figure_t *eps, double dx, double *dhP)
  { double scale = (eps->hMax - eps->hMin)/(eps->xMax - eps->xMin);
    (*dhP) = scale * dx;
  }

void epswr_h_to_x_dist(epswr_figure_t *eps, double dh, double *dxP)
  { double scale = (eps->xMax - eps->xMin)/(eps->hMax - eps->hMin);
    (*dxP) = scale * dh;
  }

void epswr_y_to_v_dist(epswr_figure_t *eps, double dy, double *dvP)
  { double scale = (eps->vMax - eps->vMin)/(eps->yMax - eps->yMin);
    (*dvP) = scale * dy;
  }
 
void epswr_v_to_y_dist(epswr_figure_t *eps, double dv, double *dyP)
  { double scale = (eps->yMax - eps->yMin)/(eps->vMax - eps->vMin);
    (*dyP) = scale * dv;
  }

void epswr_xy_to_hv_dist(epswr_figure_t *eps, double dxy, double *dhvP)
  { double h_from_x_scale = fabs((eps->hMax - eps->hMin)/(eps->xMax - eps->xMin));
    double v_from_y_scale = fabs((eps->vMax - eps->vMin)/(eps->yMax - eps->yMin));
    double scale = sqrt(h_from_x_scale*v_from_y_scale);
    (*dhvP) = scale * dxy;
  }
    
void epswr_hv_to_xy_dist(epswr_figure_t *eps, double dhv, double *dxyP)
  { double x_from_h_scale = fabs((eps->xMax - eps->xMin)/(eps->hMax - eps->hMin));
    double y_from_v_scale = fabs((eps->yMax - eps->yMin)/(eps->vMax - eps->vMin));
    double scale = sqrt(x_from_h_scale*y_from_v_scale);
    (*dxyP) = scale * dhv;
  }

void epswr_set_verbose(epswr_figure_t *eps, const bool_t verbose)
  { eps->verbose = verbose; }

void epswr_comment(epswr_figure_t *eps, const char *title)
  { epswr_dev_comment(eps, title); }
 
void epswr_show_stack(epswr_figure_t *eps, int32_t code)
  { epswr_dev_show_stack(eps, code); }
 
void epswr_flush (epswr_figure_t *eps)
  { FILE *wr = eps->wr;
    affirm(wr != NULL, "no wr");
    fflush(wr);
  }

/* UTILITIES */

void epswr_check_param(const char *name, double z, double zMin, double zMax)
  { if ((z < zMin) || (z > zMax))
      { fprintf(stderr, "** error: %s = %8.3f not in [%8.3f __ %8.3f]\n", name, z, zMin, zMax);
        assert(FALSE);
      }
  }

double epswr_round_to_nice(double x)
  { /* Special cases for 0 and infinity: */
    if ((x == 0) || (x == INF)) { return x; }
    /* Remove the sign of {x} and save it in {s}: */
    double s = (x < 0 ? -1 : +1);
    x = fabs(x);
    /* Compensate for rounding errors in the algorithm below: */
    x = 0.9999999 * x;
    /* If {x} is too small, don't bother: */
    if (x < 1.0e-300) { return s*1.0e-300; }
    /* Finds the smallest power {z} of 10 such that {z > x}: */
    double z = 1.0;
    while (z/10.0 > x) { z /= 10.0; }
    while (z <= x ) { z *= 10.0; }
    /* Check for overflow: */
    if (z == INF) { return s * INF; }
    /* Finds a nice multiple of that power of 10: */
    if (x < 0.20 * z)
      { return 0.20 * s * z; }
    else if (x < 0.25 * z)
      { return 0.25 * s * z; }
    else if (x < 0.5 * z)
      { return 0.5 * s * z; }
    else
      { return s * z; }
  }

double epswr_round_dim_mm(double x_mm, double dpi)
  {
    if (dpi == 0.0)
      { return x_mm; }
    else if (x_mm == 0.0)
      { return 0.0; }
    else
      { double dp_per_mm = (dpi/25.4)/2.0;   /* Dot pairs per  mm. */
        double x_dp_raw = x_mm * dp_per_mm;  /* Dimension in dot pairs, unrounded. */
        double x_dp = floor(x_dp_raw + 0.5); /* Dimension in dot pairs, rounded; */
        if (fabs(x_dp) == 0) { x_dp = (x_mm < 0 ? -1.0 : 1.0); };  /* Ensure min size. */
        return x_dp/dp_per_mm;
      }
  }
