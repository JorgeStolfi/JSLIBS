/* Last edited on 2011-05-15 22:49:41 by stolfi */

#define PROG_NAME "test_camera_plot"
#define PROG_DESC "tests the initial camera calibration guess algorithms"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include <r2.h>
#include <pswr.h>
#include <jsrandom.h>
#include <tf_camera_plot.h>

void write_test_plots(int Nx, int Ny, char *out_dir);

int main (int argc, char *argv[])
{
  fprintf(stderr, "Starting:\n");
  char *out_dir = "out";
  write_test_plots(512, 256, out_dir);
  fprintf(stderr, "Done.\n");
  return 0;
}

void write_test_plots(int Nx, int Ny, char *out_dir)
  {
    /* Camera aspect ratio: */
    double aspect = Nx/Ny;

    /* Postscript figure sizeand margin, in pt: */
    double hSize = 400.0*sqrt(aspect);
    double vSize = 400.0/sqrt(aspect);
    double hMrg = 4;
    double vMrg = 3;

    /* Open the output Encapsulated Postscript file: */
    char *ps_prefix = NULL;
    asprintf(&ps_prefix, "%s/", out_dir);
    PSStream *ps = pswr_new_stream(ps_prefix, NULL, TRUE, NULL, NULL, FALSE, hSize + 2*hMrg, vSize + 2*vMrg);

    pswr_new_canvas(ps, "test");

    /* Negate all Y coordinates to get the proper orientation for the Y axis: */
    pswr_set_window(ps, 0.0, Nx,  -Ny, 0.0, hMrg, hSize+hMrg, vMrg, vSize+vMrg);
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.10, 0,0); /* Black lines */
    /* Draw a frame around the image: */
    pswr_frame(ps);
    /* Draw the axes: */
    pswr_set_pen(ps, 0.0,0.0,1.0, 0.25, 0,0); /* Blue lines */
    pswr_axis(ps, HOR, -0.5*Ny, 0.0, +Nx);
    pswr_axis(ps, VER, +0.5*Nx, -Ny, 0.0);

    /* Plot random points in various styles: */
    /* frgb_t mColor = (frgb_t){{ 1.000, 1.000, 0.400 }}; */ /* Mark fill color */
    frgb_t mColor = (frgb_t){{ -1.00, -1.00, -1.00 }}; /* Mark fill color */
    pswr_set_fill_color(ps, mColor.c[0], mColor.c[1], mColor.c[2]);
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.20, 0,0); /* Black lines */

    pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0,0); /* Thick black lines */
    int nq = 10;
    int maxStyle = 5;
    r2_t q[nq];
    int m;
    for (m = 0; m <= maxStyle; m++) 
      { /* First point is at center to check alignment: */
        q[0] = (r2_t){{ 0.5*Nx, 0.5*Ny }};
        /* Second point is on X axis, position depends on mark: */
        double r = 0.1 + 0.3*((double)m)/((double)maxStyle);
        q[1] = (r2_t){{ (0.5 + r)*Nx, 0.5*Ny }};
        /* Other points are random: */
        int i;
        for (i = 2; i < nq; i++) { q[i] = (r2_t){{ 0.5*drandom()*Nx, drandom()*Ny }}; }
        tf_plot_marks(ps, nq, q, 1.0, m);
      }

    /* Finish the postscript file: */
    pswr_close_stream(ps);
  }
