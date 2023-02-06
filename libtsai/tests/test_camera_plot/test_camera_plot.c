/* Last edited on 2023-02-04 06:53:00 by stolfi */

#define PROG_NAME "test_camera_plot"
#define PROG_DESC "tests the initial camera calibration guess algorithms"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include <r2.h>
#include <epswr.h>
#include <jsrandom.h>
#include <tf_camera_plot.h>

void write_test_plots(int32_t Nx, int32_t Ny, char *out_dir);

int32_t main (int32_t argc, char *argv[])
{
  fprintf(stderr, "Starting:\n");
  char *out_dir = "out";
  write_test_plots(512, 256, out_dir);
  fprintf(stderr, "Done.\n");
  return 0;
}

void write_test_plots(int32_t Nx, int32_t Ny, char *out_dir)
  {
    /* Camera aspect ratio: */
    double aspect = Nx/Ny;

    /* Postscript figure sizeand margin, in pt: */
    double hSize = 400.0*sqrt(aspect);
    double vSize = 400.0/sqrt(aspect);
    double hMrg = 4;
    double vMrg = 3;

    /* Open the output Encapsulated Postscript file: */
    bool_t verbose = TRUE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( out_dir, NULL, "test", -1, NULL,
        hSize, vSize, hMrg, hMrg, vMrg, vMrg, 
        verbose
      );

    /* Negate all Y coordinates to get the proper orientation for the Y axis: */
    epswr_set_window(eps, hMrg, hSize+hMrg, vMrg, vSize+vMrg, FALSE, 0.0, Nx,  -Ny, 0.0);
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.10, 0,0); /* Black lines */
    /* Draw a frame around the image: */
    epswr_frame(eps);
    /* Draw the axes: */
    epswr_set_pen(eps, 0.0,0.0,1.0, 0.25, 0,0); /* Blue lines */
    epswr_axis(eps, epswr_axis_HOR, -0.5*Ny, 0.0, +Nx);
    epswr_axis(eps, epswr_axis_VER, +0.5*Nx, -Ny, 0.0);

    /* Plot random points in various styles: */
    /* frgb_t mColor = (frgb_t){{ 1.000, 1.000, 0.400 }}; */ /* Mark fill color */
    frgb_t mColor = (frgb_t){{ -1.00, -1.00, -1.00 }}; /* Mark fill color */
    epswr_set_fill_color(eps, mColor.c[0], mColor.c[1], mColor.c[2]);
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.20, 0,0); /* Black lines */

    epswr_set_pen(eps, 0.0,0.0,0.0, 0.25, 0,0); /* Thick black lines */
    int32_t nq = 10;
    int32_t maxStyle = 5;
    r2_t q[nq];
    int32_t m;
    for (m = 0; m <= maxStyle; m++) 
      { /* First point is at center to check alignment: */
        q[0] = (r2_t){{ 0.5*Nx, 0.5*Ny }};
        /* Second point is on X axis, position depends on mark: */
        double r = 0.1 + 0.3*((double)m)/((double)maxStyle);
        q[1] = (r2_t){{ (0.5 + r)*Nx, 0.5*Ny }};
        /* Other points are random: */
        int32_t i;
        for (i = 2; i < nq; i++) { q[i] = (r2_t){{ 0.5*drandom()*Nx, drandom()*Ny }}; }
        tf_plot_marks(eps, nq, q, 1.0, m);
      }

    /* Finish the postscript file: */
    epswr_end_figure(eps);
  }
