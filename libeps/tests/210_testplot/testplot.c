#define PROG_NAME "testplot"
#define PROG_DESC "test of {epswr.h} plotting ops"
#define PROG_VERS "1.0"
/* Last edited on 2023-02-21 12:18:32 by stolfi */

#define testplot_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

/* Created by J. Stolfi, UNICAMP sometime before 2003-09-30. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include <jsfile.h>
#include <bool.h>
#include <affirm.h>
#include <jsstring.h>

#include <epswr.h>
#include <epswr_dim.h>

int32_t main (int32_t argc, char **argv);

void DoTests(void);

void DrawThings(epswr_figure_t *epsf);
void DrawLabels(epswr_figure_t *epsf, double xc, double yc);
void DrawTexts(epswr_figure_t *epsf, double xc, double yc);
void DrawLines(epswr_figure_t *epsf, double xc, double yc, bool_t arrowheads);
void DrawFigures(epswr_figure_t *epsf, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed);
void DrawFrame(epswr_figure_t *epsf, double xlo, double xhi, double ylo, double yhi);
void DrawPolygon(epswr_figure_t *epsf, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed);
void DrawRoundedPolygon(epswr_figure_t *epsf, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed);
void DrawBezierPolygon(epswr_figure_t *epsf, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed);
void DrawRoundedFrame(epswr_figure_t *epsf, double xlo, double xhi, double ylo, double yhi);

int32_t main (int32_t argc, char **argv)
  { DoTests();
    return 0;
  }
  
void DoTests(void)
  {
    double hPlotSize = 800.0;
    double vPlotSize = 400.0;
    double hvMarg = 8.0;
    bool_t verbose = TRUE;

    /* Create the figure file: */
    char *fileName = NULL;
    asprintf(&fileName, "out/fig.eps");
    FILE *wr = open_write(fileName, TRUE);
    free(fileName);

    /* Create the figure object: */
    epswr_figure_t *epsf = epswr_new_figure
      ( wr, hPlotSize, vPlotSize, 
        hvMarg, hvMarg, hvMarg, hvMarg, 
        verbose    
      );
    epswr_set_client_window(epsf, -29.00, +29.00,  -16.50, +16.50);

    DrawThings(epsf);
    epswr_end_figure(epsf);
  }

void DrawThings(epswr_figure_t *epsf)
  { epswr_comment(epsf, "Thick solid red frame:");
    epswr_set_pen(epsf, 1.000, 0.000, 0.000,  0.40,  0.0, 0.0);
    epswr_frame(epsf);
    
    epswr_comment(epsf, "Medium solid black coordinate lines:");
    epswr_set_pen(epsf, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
    epswr_coord_line(epsf, epswr_axis_HOR, 0.17);
    epswr_coord_line(epsf, epswr_axis_VER, 3.14);

    epswr_comment(epsf, "Medium solid light blue coordinate line grid:");
    epswr_set_pen(epsf, 0.500, 0.750, 1.000,  0.20,  0.0, 0.0);
    epswr_coord_lines(epsf, epswr_axis_HOR, 1.0, 2.0);
    epswr_coord_lines(epsf, epswr_axis_VER, 2.0, 4.0);

    epswr_comment(epsf, "Thin dashed light yellow gridlines:");
    epswr_set_pen(epsf, 0.800, 0.800, 0.300,  0.10,  2.0, 1.0);
    epswr_grid_lines(epsf, 44, 33);

    epswr_comment(epsf, "Labels in various positions:");
    DrawLabels(epsf, -27.0, -15.5);

    epswr_comment(epsf, "Text in various positions:");
    DrawTexts(epsf, -13.0, -15.5);

    epswr_comment(epsf, "Medium solid black segments:");
    epswr_set_pen(epsf, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
    DrawLines(epsf,  +1.0, -15.5, FALSE);

    epswr_comment(epsf, "Thicker blue segments with arrowheads:");
    epswr_set_pen(epsf, 0.000, 0.000, 1.000,  0.40,  0.0, 0.0);
    DrawLines(epsf, +15.0, -15.5, TRUE);

    epswr_comment(epsf, "Thin solid black figures, closed, yellow filled:");
    epswr_set_pen(epsf, 0.000, 0.000, 0.000,  0.10,  0.0, 0.0);
    epswr_set_fill_color(epsf, 1.000, 1.000, 0.000);
    DrawFigures(epsf, -27.0,  +0.5,  TRUE, TRUE, FALSE, TRUE);

    epswr_comment(epsf, "Medium solid red figures, open, unfilled:");
    epswr_set_pen(epsf, 1.000, 0.000, 0.000,  0.20,  0.0, 0.0);
    epswr_set_fill_color(epsf, -1.00, -1.00, -1.00);
    epswr_grid_cell(epsf,  5, 44, 2, 33,  TRUE, TRUE);
    DrawFigures(epsf, -14.0,  +0.5,  TRUE, TRUE, FALSE, FALSE);

    epswr_comment(epsf, "Unstroked figures, closed, pink e-o-filled:");
    epswr_set_pen(epsf, 0.000, 0.000, 0.000,  0.50,  0.0, 0.0);
    epswr_set_fill_color(epsf, 1.000, 0.800, 0.700);
    epswr_grid_cell(epsf,  7, 44, 2, 33,  TRUE, FALSE);
    DrawFigures(epsf,  +1.0,  +0.5,  TRUE, FALSE, TRUE, TRUE);
  }

void DrawLabels(epswr_figure_t *epsf, double xc, double yc)
  {
    /* Usable area [0 _ 12]×[0 _ 15] */
    DrawFrame(epsf, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);
    
    auto void do_lab
      ( char *fontname, double fontsize, 
        char *text, 
        char *strut, 
        double xd, double yd, 
        double rot, double xalign, double yalign, 
        double R, double G, double B
      );
    
    do_lab("Courier-Bold", 10.0, "lb-g-CB10", "g",   +6.50,   +9.50, +180.0, 0.0,0.0, 0.000,0.800,0.000);
    do_lab("Courier",       8.0, "lb-x-C8",   "x",   +1.00,   +6.00,    0.0, 0.0,0.0, 1.000,0.000,0.000);
    do_lab("Times-Roman",  12.0, "cb-X-TR12", "X",   +6.00,   +5.50,  +90.0, 0.5,0.0, 1.000,0.000,0.000);
    do_lab("Helvetica",     6.0, "rc-Xg-H6",  "Xg", +11.00,  +12.00,  -45.0, 1.0,0.5, 0.000,0.000,1.000);
    
    return;
    
    void do_lab
      ( char *fontname, double fontsize, 
        char *text, 
        char *strut,
        double xd, double yd, 
        double rot, double xalign, double yalign, 
        double R, double G, double B
      )
      { 
        epswr_set_label_font(epsf, fontname, fontsize);
    
        double h = fontsize/8.0; /* Estimated intrline displacement. */
        double ddx = +h*sin(rot*M_PI/180);
        double ddy = -h*cos(rot*M_PI/180);
        
        epswr_set_pen(epsf, R, G, B, 0.05,  0.0, 0.0);
        
        epswr_set_fill_color(epsf, 0.00,0.00,0.00);
        epswr_dot(epsf, xc+xd, yc+yd, 0.2,  TRUE, FALSE);

        epswr_set_fill_color(epsf, 1-R, 1-G, 1-B);
        epswr_label(epsf, txtcat("1-",text), strut, xc+xd, yc+yd,  rot, FALSE, xalign, yalign, TRUE, FALSE);
        
        xd += ddx; yd += ddy;
        
        epswr_set_fill_color(epsf, 0.00,0.00,0.00);
        epswr_dot(epsf, xc+xd, yc+yd, 0.2,  TRUE, FALSE);
      
        epswr_set_fill_color(epsf, 1-R, 1-G, 1-B);
        epswr_label(epsf, txtcat("2-",text), strut,  xc+xd, yc+yd,  rot, FALSE, xalign, yalign, TRUE, FALSE);
        
        xd += ddx; yd += ddy;
        
        epswr_set_fill_color(epsf, 0.00,0.00,0.00);
        epswr_dot(epsf, xc+xd, yc+yd, 0.2,  TRUE, FALSE);
        
        epswr_set_fill_color(epsf, 1-R, 1-G, 1-B);
        epswr_label(epsf, txtcat("3-",text), strut,  xc+xd, yc+yd,  rot, FALSE, xalign, yalign, FALSE, TRUE);
        
        xd += ddx; yd += ddy;
        
        epswr_set_fill_color(epsf, 0.00,0.00,0.00);
        epswr_dot(epsf, xc+xd, yc+yd, 0.2,  TRUE, FALSE);
      
        epswr_set_fill_color(epsf, 1-R, 1-G, 1-B);
        epswr_label(epsf, txtcat("4-",text), strut,  xc+xd, yc+yd,  rot, FALSE, xalign, yalign, TRUE, TRUE);
      }

  }

void DrawTexts(epswr_figure_t *epsf, double xc, double yc)
  {
    /* Usable area [0 _ 12]×[0 _ 15] */
    DrawFrame(epsf, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);
    
    auto void set_text
      ( char *fontname, double fontsize, 
        double xMin, double xMax,
        double yMin, double yMax,
        double rot, 
        double R, double G, double B
      );
    
    set_text("Courier",       7.0,    2.0,  6.0,  2.0,  5.0,  +30.0, 1.000,0.000,0.000);
    epswr_text(epsf, "Left 1\nLeftius 2", FALSE, 0.0, TRUE, FALSE);
    epswr_text(epsf, "xaman",             FALSE, 0.0, TRUE, FALSE);
    epswr_text(epsf, "rugga \nmuGGa",     FALSE, 1.0, TRUE, FALSE);

    set_text("Times-Roman",   8.0,    3.0,  9.0,  7.5, 10.7,  -30.0, 0.000,0.700,0.000);
    epswr_text(epsf, "Cent 1\nCentrum 2", FALSE, 0.5, TRUE, FALSE);
    epswr_text(epsf, "Centesim 3",        FALSE, 0.5, TRUE, FALSE);
    epswr_text(epsf, "Census 4\nCena 5",  FALSE, 0.5, TRUE, FALSE);

    set_text("Helvetica",     6.0,    7.5, 10.0,  11.5, 14.0,  000.00, 0.000,0.300,1.000);
    epswr_text(epsf, "Rite 1\nRighsky", FALSE, 1.0, TRUE, FALSE);
    epswr_text(epsf, "riGhsky ",          FALSE, 1.0, TRUE, FALSE);
    epswr_text(epsf, "ramen \nRAMEN",     FALSE, 1.0, TRUE, FALSE);
    
    return;
    
    void set_text
      ( char *fontname, double fontsize, 
        double xMin, double xMax,
        double yMin, double yMax,
        double rot, 
        double R, double G, double B
      )
      { 
        /* Set the text geometry, font, color: */
        epswr_set_text_geometry(epsf, TRUE, xMin+xc, xMax+xc, yMin+yc, yMax+yc, rot);
        epswr_set_text_font(epsf, fontname, fontsize);
        epswr_set_fill_color(epsf, R,G,B);
        
        /* Compute frame (unrotated) relative to center: */
        double dx = (xMax - xMin)/2; /* Half-width of frame. */
        double dy = (yMax - yMin)/2; /* Half-height of frame. */

        double xctr = (xMin + xMax)/2 + xc;
        double yctr = (yMin + yMax)/2 + yc;
        double xp[4], yp[4];
        
        /* Draw rotated frame: */
        double tmrg = 0.2;
        for (int32_t km = 0; km < 2; km++)
          { 
            double dm = km*tmrg;
            double ang = rot*M_PI/180; /* Rotation angle in radians. */
            double sa = sin(ang), ca = cos(ang);
            double dhx = +(dx + dm)*ca, dvx = +(dx + dm)*sa;
            double dhy = -(dy + dm)*sa, dvy = +(dy + dm)*ca;

            xp[0] = xctr - dhx - dhy;  yp[0] = yctr - dvx - dvy;
            xp[1] = xctr + dhx - dhy;  yp[1] = yctr + dvx - dvy;
            xp[2] = xctr + dhx + dhy;  yp[2] = yctr + dvx + dvy;
            xp[3] = xctr - dhx + dhy;  yp[3] = yctr - dvx + dvy;
        
            if (km == 0)
              { epswr_set_pen(epsf, 0.800,0.800,0.800, 0.15, 0,0); }
            else
              { epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0,0); }
            epswr_polygon(epsf, TRUE, xp, yp, 4,   FALSE,TRUE, TRUE);
          }
        return;
      }

  }

void DrawLines(epswr_figure_t *epsf, double xc, double yc, bool_t arrowheads)
  {
    /* Usable area [0 _ 12]×[0 _ 15] */
    DrawFrame(epsf, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);
    
    auto void do_seg(double xa, double ya, double b, double yb);
    void do_seg(double xa, double ya, double xb, double yb)
      { epswr_segment(epsf, xa+xc, ya+yc, xb+xc, yb+yc);
        epswr_dot(epsf, xa+xc, ya+yc, 1.0, FALSE, TRUE);
        epswr_dot(epsf, xb+xc, yb+yc, 1.0, FALSE, TRUE);
        if (arrowheads)
          { epswr_arrowhead(epsf, xa+xc, ya+yc, xb+xc, yb+yc, 2.0, 3.0, 0.85, TRUE, TRUE); }
      }
    
    do_seg(1.0, 1.0, 11.0, 3.0);
    do_seg(1.0, 3.0, 11.0, 1.0);
    
    auto void dodim(double xp, double yp);
    
    void dodim(double xp, double yp)
      { 
        double xa = xc+xp+0.0, ya = yc+yp+0.0, xb = xc+xp+0.0, yb = yc+yp+2.0;
        double agap = 0.5,  bgap = 1.5;
        double xa1 = xa-agap, ya1 = ya, xb1 = xb-bgap, yb1 = yb; 
        epswr_set_fill_color(epsf, 0.00,0.00,0.50);
        epswr_dot(epsf, xa, ya, 0.5, TRUE, FALSE);
        epswr_dot(epsf, xb, yb, 0.5, TRUE, FALSE);
        epswr_dot(epsf, xa1, ya1, 0.3, TRUE, FALSE);
        epswr_dot(epsf, xb1, yb1, 0.3, TRUE, FALSE);
        
        epswr_set_pen(epsf, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
        double dab, xr, yr, rot;
        epswr_dim_linear(epsf, xa, ya, xb, yb, &dab, agap, bgap, 2.0, FALSE, 1.0,5.0, 0.0,2.0, &xr, &yr, &rot);
        
        epswr_set_fill_color(epsf, 1.00,0.00,0.00);
        epswr_dot(epsf, xr, yr, 0.7, TRUE, FALSE);
      }

    epswr_segment(epsf,  1.0+xc,  5.0+yc,  1.0+xc, 13.0+yc);
    epswr_segment(epsf, 11.0+xc,  5.0+yc, 11.0+xc, 13.0+yc);
    epswr_segment(epsf,  2.0+xc, 14.0+yc, 10.0+xc, 14.0+yc);
    int32_t i;
    for (i = -7; i <= +7; i++)
      { double align = ((i % 2) == 0)*0.5 - ((i % 4) == 0)*0.25;
        double ticsz = 0.5 + ((i % 2) == 0)*0.5 + ((i % 4) == 0)*1.0; 
        epswr_tic(epsf, epswr_axis_VER,  1.0+xc, 9.0+0.5*i+yc, ticsz, align); 
        epswr_tic(epsf, epswr_axis_VER, 11.0+xc, 9.0+0.5*i+yc, ticsz, 1-align); 
        epswr_tic(epsf, epswr_axis_HOR, 6.0 + 0.5*i+xc, 14.0+yc, ticsz, 1-align); 
      }
    
    epswr_curve(epsf,
       4.0+xc,  4.0+yc,  
      11.0+xc, 14.0+yc, 
       1.0+xc, 14.0+yc, 
       8.0+xc,  4.0+yc
    );
    epswr_square(epsf,  4.0+xc,  4.0+yc, 1.0, FALSE, TRUE);  
    epswr_dot   (epsf, 11.0+xc, 14.0+yc, 0.5, FALSE, TRUE);
    epswr_dot   (epsf,  1.0+xc, 14.0+yc, 0.5, FALSE, TRUE); 
    epswr_square(epsf,  8.0+xc,  4.0+yc, 1.0, FALSE, TRUE);
    
    dodim(+10.0, +10.0);

  }
  
void DrawFigures(epswr_figure_t *epsf, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed)
  {
    /* Usable area [0 _ 12]×[0 _ 15] */
    DrawFrame(epsf, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);

    epswr_rectangle(epsf, 1.5+xc, 4.5+xc, 10.0+yc, 11.5+yc, fill, draw);
    epswr_circle(epsf, 3.0+xc, 8.0+yc, 1.75, fill, draw);
    epswr_lune(epsf, 9.0+xc, 7.0+yc, 1.5, 45.0, fill, draw);
    epswr_slice(epsf, 3.0+xc, 12.0+yc, 2.0, 60.0, 165.0, fill, draw);
    epswr_quadrilateral(epsf, 5.0+xc, 8.0+yc, 9.0+xc, 10.0+yc, 7.0+xc, 12.0+yc, 6.0+xc, 9.0+yc, fill, draw);
    epswr_triangle(epsf, 7.0+xc, 10.0+yc, 11.0+xc, 12.0+yc, 9.0+xc, 14.0+yc, fill, draw);

    double rmk = 1.5;
    epswr_dot(epsf, 3.0+xc, 8.0+yc, rmk, fill, draw);
    epswr_square(epsf, 5.0+xc, 6.0+yc, rmk, fill, draw);
    epswr_diamond(epsf, 6.5+xc, 6.0+yc, rmk, rmk, fill, draw);
    epswr_diamond(epsf, 10.0+xc, 9.0+yc, rmk, 2.0, fill, draw);
    
    epswr_cross(epsf, 5.5+xc, 13.0+yc, rmk, FALSE, draw);
    epswr_cross(epsf, 7.0+xc, 13.0+yc, rmk, TRUE, draw);
    epswr_asterisk(epsf, 5.5+xc, 11.0+yc, rmk, draw);
    
    /* Ordinary polygon {eo} filled: */
    DrawPolygon(epsf, 2.5+xc, 3.0+yc, fill, draw, eo, closed);
    
    /* Bézier polygon {eo} filled: */
    DrawBezierPolygon(epsf, 6.0+xc, 3.0+yc, fill, draw, eo, closed);
    
    /* Rounded polygon {eo} filled: */
    DrawRoundedPolygon(epsf, 9.5+xc, 3.0+yc, fill, draw, eo, closed);

    /* Rounded frame around area: */
    DrawRoundedFrame(epsf, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);
  }

void DrawFrame(epswr_figure_t *epsf, double xlo, double xhi, double ylo, double yhi)
  {
    epswr_comment(epsf, "enter DrawFrame");
    epswr_rectangle(epsf, xlo+0.25, xhi-0.25, ylo+0.25, yhi-0.25, FALSE, TRUE);
    epswr_comment(epsf, "exit DrawFrame");
  }

void DrawPolygon(epswr_figure_t *epsf, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed)
  {
    double r = 1.4;
    int32_t n = 7;
    double x[n], y[n];
    int32_t i;
    for (i = 0; i < n; i++)
      { double t = 2.0*M_PI*((double)i)/((double)n);
        x[i] = xc + r*cos(2*t); 
        y[i] = yc + 1.25*r*sin(2*t);
      }
    epswr_polygon(epsf, closed, x, y, n, fill, draw, eo);
  }

void DrawRoundedPolygon(epswr_figure_t *epsf, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed)
  {
    double ra = 1.4;
    double rb = 0.7;
    int32_t n = 10;
    double x[n], y[n];
    int32_t i;
    for (i = 0; i < n; i += 2)
      { double s = 2.0*M_PI*((double)i)/((double)n);
        x[i] = xc + ra*cos(s); 
        y[i] = yc + 1.25*ra*sin(s);
        double t = 2.0*M_PI*((double)i+1)/((double)n);
        x[i+1] = xc + rb*cos(t); 
        y[i+1] = yc + 1.25*rb*sin(t);
      }
    double round = 0.20;
    epswr_rounded_polygon(epsf, closed, x, y, n, round, fill, draw, eo);
  }

void DrawRoundedFrame(epswr_figure_t *epsf, double xlo, double xhi, double ylo, double yhi)
  {
    int32_t n = 10;
    double x[n], y[n];
    x[0] = xlo+0.50; y[0] = ylo+0.50;
    x[1] = xhi-0.50; y[1] = ylo+0.50;
    x[2] = xhi-0.50; y[2] = yhi-1.50;
    x[3] = xhi-1.50; y[3] = yhi-1.50;
    x[4] = xhi-1.50; y[4] = yhi-1.50;
    x[5] = xhi-1.50; y[5] = yhi-0.50;
    x[6] = xlo+1.50; y[6] = yhi-0.50;
    x[7] = xlo+1.50; y[7] = yhi-1.50;
    x[8] = xlo+1.50; y[8] = yhi-1.50;
    x[9] = xlo+0.50; y[9] = yhi-1.50;
    double round = 0.75;
    epswr_rounded_polygon(epsf, TRUE, x, y, n, round, FALSE, TRUE, FALSE);
  }

void DrawBezierPolygon(epswr_figure_t *epsf, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed)
  {
    /* Prototypical S-shaped arc: */
    double xs[4], ys[4];
    xs[0] = +1.200; ys[0] = 00.000;
    xs[1] = +0.400; ys[1] = -0.800;
    xs[2] = -0.400; ys[2] = +0.800;
    xs[3] = -1.200; ys[3] = 00.000;
    
    /* Prototypical convex arc: */
    double xv[4], yv[4];
    xv[0] = +1.200; yv[0] = 00.000;
    xv[1] = +0.400; yv[1] = +0.800;
    xv[2] = -0.400; yv[2] = +0.800;
    xv[3] = -1.200; yv[3] = 00.000;
    
    double mag = 1.500;   /* Magnification factor. */
    double shift = 0.600; /* Side shift from center. */
    
    int32_t n = 5;
    int32_t np = 4*n;
    double x[np], y[np];
    int32_t i, j;
    for (i = 0; i < n; i++)
      { /* Rotation angle: */
        double t = 2.0*M_PI*((double)i)/((double)n);
        double ct = cos(2*t), st = sin(2*t);
        for (j = 0; j < 4; j++)
          { /* Pick Bezier point from chosen arc: */
            double xb, yb;
            if (i != 2)
              { xb = xs[j]; yb = ys[j]; }
            else
              { xb = xv[j]; yb = yv[j]; }
            /* Displace in {y}: */
            yb += shift;
            /* Rotate and add center: */
            x[4*i + j] = xc + mag*(+ ct*xb - st*yb);
            y[4*i + j] = yc + mag*(+ st*xb + ct*yb);
          }
      }
    
    /* Draw the control polygon: */
    /* epswr_polygon(epsf, TRUE, x, y, np, FALSE, draw, FALSE); */
    
    /* Draw the Bézier polygon: */
    epswr_bezier_polygon(epsf, closed, x, y, n, fill, draw, eo);
  }
