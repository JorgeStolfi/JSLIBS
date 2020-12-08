#define PROG_NAME "testplot"
#define PROG_DESC "test of {pswr.h} (minus automatic layout)"
#define PROG_VERS "1.0"

/* Last edited on 2019-05-06 07:03:24 by jstolfi */

#define testplot_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

/* Created by J. Stolfi, UNICAMP sometime before 2003-09-30. */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <jsstring.h>
#include <pswr.h>
#include <pswr_iso.h>

#define OUT_PREFIX "out/"
  /* Prefix for output file names. */

int main (int argc, char **argv);
void DoEPSTests(void);
void DoPSTests(void);
void DrawThings(PSStream *ps);
void DrawTexts(PSStream *ps, double xc, double yc);
void DrawLines(PSStream *ps, double xc, double yc, bool_t arrowheads);
void DrawFigures(PSStream *ps, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed);
void DrawFrame(PSStream *ps, double xlo, double xhi, double ylo, double yhi);
void DrawPolygon(PSStream *ps, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed);
void DrawRoundedPolygon(PSStream *ps, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed);
void DrawBezierPolygon(PSStream *ps, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed);
void DrawRoundedFrame(PSStream *ps, double xlo, double xhi, double ylo, double yhi);

int main (int argc, char **argv)
  { DoEPSTests();
    DoPSTests();
    return 0;
  }
  
void DoEPSTests(void)
  {
    PSStream *ps = pswr_new_stream(OUT_PREFIX, NULL, TRUE, "fig", NULL, FALSE, 328.0, 246.0);
    pswr_new_canvas(ps, NULL);
    pswr_set_window
      ( ps, 
        -22.00, +22.00,  -16.50, +16.50, 
          4.00, 324.00,    3.00, 243.00
      );
    pswr_set_grid(ps, 44, 33);
    DrawThings(ps);
    pswr_close_stream(ps);
  }
      
void DoPSTests(void)
  { 
    double hSize, vSize; /* Paper dimensions in pt. */
    pswr_get_paper_dimensions("letter", TRUE, &hSize, &vSize);
    
    PSStream *ps = pswr_new_stream(OUT_PREFIX, NULL, FALSE, "doc", "letter", FALSE, 0, 0);
    pswr_new_canvas(ps, "one");
      pswr_set_window
        ( ps, 
          -22.00, +22.00,  -16.50, +16.50, 
          144.00, 464.00,  144.00, 384.00
        );
      pswr_set_grid(ps, 44, 33);
      DrawThings(ps);
      pswr_set_pen(ps, 0.000, 0.000, 1.000,  0.15,  0.0,0.0);
      pswr_add_caption(ps, "First Line Left Aligned\nSecond Line Left Aligned", 0.0);
      pswr_add_caption(ps, "Third Line Right Aligned", 1.0);
      pswr_add_caption(ps, "Fourth Line Centered", 0.5);

      pswr_set_window
        ( ps, 
          -22.00, +22.00,  -16.50, +16.50,
          144.00, 304.00,  432.00, 552.00
        );
      pswr_set_grid(ps, 44, 33);
      DrawThings(ps);
      pswr_add_caption(ps, "This Figure No Caption", 0.5);
    pswr_new_canvas(ps, "two");
      pswr_set_window
        ( ps, 
          -22.00, +22.00,  -16.50, +16.50, 
          144.00, 304.00,  432.00, 552.00
        );
      pswr_set_grid(ps, 44, 33);
      DrawThings(ps);
    pswr_close_stream(ps);
  }

void DrawThings(PSStream *ps)
  { pswr_comment(ps, "Thick solid red frame:");
    pswr_set_pen(ps, 1.000, 0.000, 0.000,  0.40,  0.0, 0.0);
    pswr_frame(ps);

    pswr_comment(ps, "Thin dashed light yellow gridlines:");
    pswr_set_pen(ps, 1.000, 1.000, 0.500,  0.10,  2.0, 1.0);
    pswr_grid_lines(ps);

    pswr_comment(ps, "Medium solid black coordinate lines:");
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
    pswr_coord_line(ps, HOR, 0.17);
    pswr_coord_line(ps, VER, 3.14);

    pswr_comment(ps, "Text in various positions:");
    DrawTexts(ps, -20.0, -15.5);

    pswr_comment(ps, "Medium solid black segments:");
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.20,  0.0, 0.0);
    DrawLines(ps,  -6.0, -15.5, FALSE);

    pswr_comment(ps, "Thicker blue segments with arrowheads:");
    pswr_set_pen(ps, 0.000, 0.000, 1.000,  0.40,  0.0, 0.0);
    DrawLines(ps,  +8.0, -15.5, TRUE);

    pswr_comment(ps, "Thin solid black figures, closed, yellow filled:");
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.10,  0.0, 0.0);
    pswr_set_fill_color(ps, 1.000, 1.000, 0.000);
    pswr_grid_cell(ps,  3, 2,  TRUE, TRUE);
    DrawFigures(ps, -20.0,  +0.5,  TRUE, TRUE, FALSE, TRUE);

    pswr_comment(ps, "Medium solid red figures, open, unfilled:");
    pswr_set_pen(ps, 1.000, 0.000, 0.000,  0.20,  0.0, 0.0);
    pswr_set_fill_color(ps, -1.00, -1.00, -1.00);
    pswr_grid_cell(ps,  5, 2,  TRUE, TRUE);
    DrawFigures(ps,  -6.0,  +0.5,  TRUE, TRUE, FALSE, FALSE);

    pswr_comment(ps, "Unstroked figures, closed, pink e-o-filled:");
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.50,  0.0, 0.0);
    pswr_set_fill_color(ps, 1.000, 0.800, 0.700);
    pswr_grid_cell(ps,  7, 2,  TRUE, FALSE);
    DrawFigures(ps,  +8.0,  +0.5,  TRUE, FALSE, TRUE, TRUE);
  }

void DrawTexts(PSStream *ps, double xc, double yc)
  {
    /* Usable area [0 _ 12]×[0 _ 15] */
    DrawFrame(ps, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);
    
    auto void do_lab
      ( char *fontname, double fontsize, 
        char *text, 
        double xd, double yd, 
        double rot, double xalign, double yalign, 
        double R, double G, double B
      );

    do_lab("Courier",      8.0, "rC8",    +1.00,   +6.00,    0.0, 0.0,0.0, 1.000,0.000,0.000);
    do_lab("Times-Roman", 12.0, "rTR12",  +6.00,   +5.50,  +90.0, 0.5,0.0, 1.000,0.000,0.000);
    do_lab("Helvetica",    6.0, "bH6",   +11.00,  +12.00,  -45.0, 1.0,0.5, 0.000,0.000,1.000);
    do_lab("CourierBold", 10.0, "gCB10",  +6.50,   +9.50, +180.0, 0.0,0.5, 0.000,0.800,0.000);
    
    return;
    
    void do_lab
      ( char *fontname, double fontsize, 
        char *text, 
        double xd, double yd, 
        double rot, double xalign, double yalign, 
        double R, double G, double B
      )
      { 
        pswr_set_label_font(ps, fontname, fontsize);
    
        double h = fontsize/8.0; /* Estimated intrline displacement. */
        double ddx = +h*sin(rot*M_PI/180);
        double ddy = -h*cos(rot*M_PI/180);
        
        pswr_set_pen(ps, R, G, B, 0.05,  0.0, 0.0);
        
        pswr_set_fill_color(ps, 0.00,0.00,0.00);
        pswr_dot(ps, xc+xd, yc+yd, 0.2,  TRUE, FALSE);

        pswr_set_fill_color(ps, 1-R, 1-G, 1-B);
        pswr_label(ps, txtcat("1-",text),  xc+xd, yc+yd,  rot, xalign, yalign);
        
        xd += ddx; yd += ddy;
        
        pswr_set_fill_color(ps, 0.00,0.00,0.00);
        pswr_dot(ps, xc+xd, yc+yd, 0.2,  TRUE, FALSE);
      
        pswr_set_fill_color(ps, 1-R, 1-G, 1-B);
        pswr_fill_draw_label(ps, txtcat("2-",text),  xc+xd, yc+yd,  rot, xalign, yalign, TRUE, FALSE);
        
        xd += ddx; yd += ddy;
        
        pswr_set_fill_color(ps, 0.00,0.00,0.00);
        pswr_dot(ps, xc+xd, yc+yd, 0.2,  TRUE, FALSE);
        
        pswr_set_fill_color(ps, 1-R, 1-G, 1-B);
        pswr_fill_draw_label(ps, txtcat("3-",text),  xc+xd, yc+yd,  rot, xalign, yalign, FALSE, TRUE);
        
        xd += ddx; yd += ddy;
        
        pswr_set_fill_color(ps, 0.00,0.00,0.00);
        pswr_dot(ps, xc+xd, yc+yd, 0.2,  TRUE, FALSE);
      
        pswr_set_fill_color(ps, 1-R, 1-G, 1-B);
        pswr_fill_draw_label(ps, txtcat("4-",text),  xc+xd, yc+yd,  rot, xalign, yalign, TRUE, TRUE);
      }
  }

void DrawLines(PSStream *ps, double xc, double yc, bool_t arrowheads)
  {
    /* Usable area [0 _ 12]×[0 _ 15] */
    DrawFrame(ps, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);
    
    auto void do_seg(double xa, double ya, double b, double yb);
    void do_seg(double xa, double ya, double xb, double yb)
      { pswr_segment(ps, xa+xc, ya+yc, xb+xc, yb+yc);
        pswr_dot(ps, xa+xc, ya+yc, 1.0, FALSE, TRUE);
        pswr_dot(ps, xb+xc, yb+yc, 1.0, FALSE, TRUE);
        if (arrowheads)
          { pswr_arrowhead(ps, xa+xc, ya+yc, xb+xc, yb+yc, 2.0, 3.0, 0.85, TRUE, TRUE); }
      }
    
    do_seg(1.0, 1.0, 11.0, 3.0);
    do_seg(1.0, 3.0, 11.0, 1.0);
    
    pswr_segment(ps,  1.0+xc,  5.0+yc,  1.0+xc, 13.0+yc);
    pswr_segment(ps, 11.0+xc,  5.0+yc, 11.0+xc, 13.0+yc);
    pswr_segment(ps,  2.0+xc, 14.0+yc, 10.0+xc, 14.0+yc);
    int i;
    for (i = -7; i <= +7; i++)
      { double align = ((i % 2) == 0)*0.5 - ((i % 4) == 0)*0.25;
        double ticsz = 0.5 + ((i % 2) == 0)*0.5 + ((i % 4) == 0)*1.0; 
        pswr_tic(ps, VER,  1.0+xc, 9.0+0.5*i+yc, ticsz, align); 
        pswr_tic(ps, VER, 11.0+xc, 9.0+0.5*i+yc, ticsz, 1-align); 
        pswr_tic(ps, HOR, 6.0 + 0.5*i+xc, 14.0+yc, ticsz, 1-align); 
      }
    
    pswr_curve(ps,
       4.0+xc,  4.0+yc,  
      11.0+xc, 14.0+yc, 
       1.0+xc, 14.0+yc, 
       8.0+xc,  4.0+yc
    );
    pswr_square(ps,  4.0+xc,  4.0+yc, 1.0, FALSE, TRUE);  
    pswr_dot   (ps, 11.0+xc, 14.0+yc, 0.5, FALSE, TRUE);
    pswr_dot   (ps,  1.0+xc, 14.0+yc, 0.5, FALSE, TRUE); 
    pswr_square(ps,  8.0+xc,  4.0+yc, 1.0, FALSE, TRUE);
  }
  
void DrawFigures(PSStream *ps, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed)
  {
    /* Usable area [0 _ 12]×[0 _ 15] */
    DrawFrame(ps, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);
    
    pswr_rectangle(ps, 1.5+xc, 4.5+xc, 10.0+yc, 11.5+yc, fill, draw);
    pswr_circle(ps, 3.0+xc, 8.0+yc, 1.75, fill, draw);
    pswr_lune(ps, 9.0+xc, 7.0+yc, 1.5, 45.0, fill, draw);
    pswr_slice(ps, 3.0+xc, 12.0+yc, 2.0, 60.0, 165.0, fill, draw);
    pswr_quadrilateral(ps, 5.0+xc, 8.0+yc, 9.0+xc, 10.0+yc, 7.0+xc, 12.0+yc, 6.0+xc, 9.0+yc, fill, draw);
    pswr_triangle(ps, 7.0+xc, 10.0+yc, 11.0+xc, 12.0+yc, 9.0+xc, 14.0+yc, fill, draw);
    
    double rmk = 1.5;
    pswr_dot(ps, 3.0+xc, 8.0+yc, rmk, fill, draw);
    pswr_square(ps, 5.0+xc, 6.0+yc, rmk, fill, draw);
    pswr_diamond(ps, 6.5+xc, 6.0+yc, rmk, rmk, fill, draw);
    pswr_diamond(ps, 10.0+xc, 9.0+yc, rmk, 2.0, fill, draw);
    
    pswr_cross(ps, 5.5+xc, 13.0+yc, rmk, FALSE, draw);
    pswr_cross(ps, 7.0+xc, 13.0+yc, rmk, TRUE, draw);
    pswr_asterisk(ps, 5.5+xc, 11.0+yc, rmk, draw);
    
    /* Ordinary polygon {eo} filled: */
    DrawPolygon(ps, 2.5+xc, 3.0+yc, fill, draw, eo, closed);
    
    /* Bézier polygon {eo} filled: */
    DrawBezierPolygon(ps, 6.0+xc, 3.0+yc, fill, draw, eo, closed);
    
    /* Rounded polygon {eo} filled: */
    DrawRoundedPolygon(ps, 9.5+xc, 3.0+yc, fill, draw, eo, closed);

    /* Rounded frame around area: */
    DrawRoundedFrame(ps, 0.00+xc, 12.00+xc, 0.00+yc, 15.00+yc);
  }

void DrawFrame(PSStream *ps, double xlo, double xhi, double ylo, double yhi)
  {
    pswr_rectangle(ps, xlo+0.25, xhi-0.25, ylo+0.25, yhi-0.25, FALSE, TRUE);
  }

void DrawPolygon(PSStream *ps, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed)
  {
    double r = 1.4;
    int n = 7;
    double x[n], y[n];
    int i;
    for (i = 0; i < n; i++)
      { double t = 2.0*M_PI*((double)i)/((double)n);
        x[i] = xc + r*cos(2*t); 
        y[i] = yc + 1.25*r*sin(2*t);
      }
    pswr_polygon(ps, closed, x, y, n, fill, draw, eo);
  }

void DrawRoundedPolygon(PSStream *ps, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed)
  {
    double ra = 1.4;
    double rb = 0.7;
    int n = 10;
    double x[n], y[n];
    int i;
    for (i = 0; i < n; i += 2)
      { double s = 2.0*M_PI*((double)i)/((double)n);
        x[i] = xc + ra*cos(s); 
        y[i] = yc + 1.25*ra*sin(s);
        double t = 2.0*M_PI*((double)i+1)/((double)n);
        x[i+1] = xc + rb*cos(t); 
        y[i+1] = yc + 1.25*rb*sin(t);
      }
    double round = 0.20;
    pswr_rounded_polygon(ps, closed, x, y, n, round, fill, draw, eo);
  }

void DrawRoundedFrame(PSStream *ps, double xlo, double xhi, double ylo, double yhi)
  {
    int n = 10;
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
    pswr_rounded_polygon(ps, TRUE, x, y, n, round, FALSE, TRUE, FALSE);
  }

void DrawBezierPolygon(PSStream *ps, double xc, double yc, bool_t fill, bool_t draw, bool_t eo, bool_t closed)
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
    
    int n = 5;
    int np = 4*n;
    double x[np], y[np];
    int i, j;
    for (i = 0; i < n; i++)
      { /* Rotation angle: */
        double t = 2.0*M_PI*((double)i)/((double)n);
        double ct = cos(2*t), st = sin(2*t);
        for (j = 0; j < 4; j++)
          { /* Pick Bezier point from chosen arc: */
            double xb, yb;
            if (i != 99)
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
    /* pswr_polygon(ps, TRUE, x, y, np, FALSE, draw, FALSE); */
    
    /* Draw the Bézier polygon: */
    pswr_bezier_polygon(ps, closed, x, y, n, fill, draw, eo);
  }
