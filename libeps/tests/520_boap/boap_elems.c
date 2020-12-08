/* See {boap_elems.h} */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <jsfile.h>
#include <epswr.h>
#include <epswr_dim.h>

#include <boap_bar.h>

#include <boap_elems.h>

epswr_figure_t *boap_new_figure
  ( char *pname, 
    char *subname, 
    double xMin, double xMax, 
    double yMin, double yMax, 
    bool_t landscape
  )
  { 
    char *fname = NULL;
    if (subname != NULL)
      { asprintf(&fname, "out/porch_%s_%s.eps", pname, subname); }
    else
      { asprintf(&fname, "out/porch_%s.eps", pname); }
    FILE *wr = open_write(fname, TRUE);
    free(fname);
   
    /* Device coordinate extent (pt; aspect ratio 13:10 or 10:13): */
    double hSize, vSize;
    if (landscape)
      { hSize = 780.0; vSize = 600.0; }
    else
      { hSize = 600.0; vSize = 780.0; }
    double hvMarg = 6.0;

    bool_t verbose = FALSE;
    epswr_figure_t *epsf = epswr_new_figure(wr, hSize, vSize, hvMarg, hvMarg, hvMarg, hvMarg, verbose);
    epswr_set_fill_color(epsf, 1.000,0.800,0.400);
    epswr_set_label_font(epsf, "Courier-Bold", 10.0);
    epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0.0,0.0);

    /* Convert object bounding box to center+size: */
    double xCtr = (xMin + xMax)/2, xSize = xMax - xMin;
    double yCtr = (yMin + yMax)/2, ySize = yMax - yMin;

    /* Round client plot area size to 13:10 or 10:13 aspect ratio: */
    if (landscape)
      { double u = ceil(fmax(xSize/13, ySize/10));
        xSize = 13*u; ySize = 10*u;
      }
    else
      { double u = ceil(fmax(xSize/10, ySize/13));
        xSize = 10*u; ySize = 13*u;
      }

    /* Plot coordinates are {x=Y}, {y=Z}: */
    double plot_xMin = xCtr - xSize/2, plot_xMax = plot_xMin + xSize;
    double plot_yMin = yCtr - ySize/2, plot_yMax = plot_yMin + ySize;
    epswr_set_client_window(epsf, plot_xMin, plot_xMax, plot_yMin, plot_yMax);

    return epsf;
  }
