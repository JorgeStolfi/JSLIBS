#define PROG_NAME "teststream"
#define PROG_DESC "test of {pswr.h} (with automatic layout)"
#define PROG_VERS "1.0"

/* Last edited on 2012-12-08 23:36:53 by stolfilocal */
/* Created by J. Stolfi, UNICAMP sometime before 2003-10-11. */

#define teststream_COPYRIGHT \
  "Copyright © 2003  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <pswr.h>
#include <bool.h>

#define OUT_PREFIX "out/"
  /* Prefix for output file names. */

int main (int argc, char **argv);
void DoEPSTests(bool_t numbered);
void DoPSTests(bool_t landscape);
void DrawPictures(PSStream *ps);
void DrawPicture(PSStream *ps, double rot);

int main (int argc, char **argv)
  { DoEPSTests(FALSE);
    DoEPSTests(TRUE);
    DoPSTests(FALSE);
    DoPSTests(TRUE);
    return 0;
  }
  
void DoEPSTests(bool_t numbered)
  {
    PSStream *ps = pswr_new_stream(OUT_PREFIX, NULL, TRUE, "fig", NULL, TRUE, 424.0, 298.0);
    int i;
    for (i = 0; i < 2; i++)
      { char *pageName = NULL;
        if (! numbered) { asprintf(&pageName, "pg%c", 'A' + i); }
        pswr_new_canvas(ps, pageName);
        ps->verbose = TRUE;
        pswr_set_canvas_layout(ps, 80.0, 80.0, FALSE, 2.0, 4.0, 1, 5, 3);
        DrawPictures(ps);
        free(pageName);
      }
    pswr_close_stream(ps);
  }
      
void DoPSTests(bool_t landscape)
  { 
    double hSize, vSize; /* Paper dimensions in pt. */
    char *paperSize = "letter";
    pswr_get_paper_dimensions(paperSize, landscape, &hSize, &vSize);
    char *docName = ( landscape ? "doc-L" : "doc-P");
    PSStream *ps = pswr_new_stream(OUT_PREFIX, NULL, FALSE, docName, paperSize, landscape, 0, 0);
    ps->verbose = TRUE;
    pswr_set_canvas_layout(ps, 80.0, 80.0, TRUE, 2.0, 4.0, 1, 3, 4);
    DrawPictures(ps);
    pswr_close_stream(ps);
  }

void DrawPictures(PSStream *ps)
  { 
    int NFrames = 15;
    int i; 
    for (i = 0; i < NFrames; i++)
      { double rot = ((double)i)/((double)NFrames);
        pswr_new_picture(ps, -1.0, 1.0, -1.0, 1.0);
        DrawPicture(ps, rot);
        if (i%2 == 0) { pswr_frame(ps); }
      }
  }

#define BUFSZ 100

void DrawPicture(PSStream *ps, double rot)
  {
    double t = rot*0.5*M_PI;
    double ct = cos(t), st = sin(t);
    char xrot[100];
    snprintf(xrot, BUFSZ, "%5.3f", rot);
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.15,  0.0,0.0);
    pswr_add_caption(ps, xrot, 0.5);
    
    double xscale = 1.05*sqrt(6.0);
    double yscale = 1.05*sqrt(18.0);
    
    pswr_set_pen(ps, 0.000, 0.000, 0.000,  0.10,  0.0,0.0);
    /* Enumerate the eight faces of the cube: */
    double v[3];
    int ax;
    for (ax = 0; ax < 3; ax++)
      { int bx = (ax + 1) % 3, cx = (ax + 2) % 3;
        int i;
        for (i = -1; i <= +1; i += 2)
          { int r = 0;
            int ar;
            double fn[3]; /* Face normal */
            double xp[4], yp[4];  /* Plot coordinates */
            /* Clear face normal: */
            for (ar = 0; ar < 3; ar++) { fn[ar] = 0.0; }
            /* Enumerate vertices of face {v[ax] == i}: */
            int j, k;
            for (j = -1; j <= +1; j += 2)
              { for (k = -1; k <= +1; k += 2)
                  { v[ax] = i;
                    v[bx] = j;
                    v[cx] = k*j;
                    /* rotate point {v} by angle {rot*2pi} around all axes: */
                    for (ar = 0; ar < 3; ar++)
                      { int br = (ar + 1) % 3, cr = (ar + 2) % 3;
                        double xt =  ct*v[br] + st*v[cr];
                        double yt = -st*v[br] + ct*v[cr];
                        v[br] = xt; v[cr] = yt;
                      }
                    /* Accumulate into face normal: */
                    for (ar = 0; ar < 3; ar++) { fn[ar] += v[ar]; }
                    /* Project and store in {xp[r],yp[r]}: */
                    xp[r] = (- v[0] + v[1])/xscale;
                    yp[r] = (- v[0] - v[1] + 2*v[2])/yscale;
                    r++;
                  }
              }
            /* Check visibility and plot: */
            if (fn[0] + fn[1] + fn[2] > 0.0)
              { /* Compute illumination: */
                double cosl = fn[2]/4.0;
                double shade = 0.5 + 0.5*(cosl <= 0.0 ? 0.0 : cosl);
                pswr_set_fill_color(ps, shade*1.000, shade*0.850, shade*0.500);
                pswr_polygon(ps, TRUE, xp, yp, 4, TRUE, TRUE, FALSE);
              }
          }
      }
  }
