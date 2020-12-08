/* See epswr_def.h */
/* Last edited on 2020-10-27 18:51:47 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <bool.h>
#include <affirm.h>

#include <epswr_def.h>

void epswr_def_set_device_window
  ( epswr_def_figure_t *eps,
    double hMin, double hMax,
    double vMin, double vMax,
    const char *requestor
  )
  {
    eps->hMin = hMin;
    eps->hMax = hMax;
    eps->vMin = vMin;
    eps->vMax = vMax;
    
    if (eps->verbose)
      { epswr_def_report_window(eps, requestor, "Device", "set to", hMin, hMax, vMin, vMax);  }
  }

void epswr_def_undefine_client_window(epswr_def_figure_t *eps)
  { epswr_def_set_client_window(eps, NAN, NAN, NAN, NAN, NULL); }

void epswr_def_reset_client_window(epswr_def_figure_t *eps)
  { double xMin = 0.0, xMax = eps->hMax - eps->hMin;
    double yMin = 0.0, yMax = eps->vMax - eps->vMin;
    epswr_def_set_client_window(eps, xMin, xMax, yMin, yMax, NULL);
  }

void epswr_def_set_client_window
  ( epswr_def_figure_t *eps,
    double xMin, double xMax,
    double yMin, double yMax,
    const char *requestor
  )
  {
    eps->xMin = xMin;
    eps->xMax = xMax;
    eps->yMin = yMin;
    eps->yMax = yMax;
    
    if (eps->verbose)
      { epswr_def_report_window(eps, requestor, "Client", "set to", xMin, xMax, yMin, yMax);  }
  }    
    
void epswr_def_get_client_window
  ( epswr_def_figure_t *eps,
    double *xMinP, double *xMaxP,
    double *yMinP, double *yMaxP
  )
  {
    (*xMinP) = eps->xMin;
    (*xMaxP) = eps->xMax;
    (*yMinP) = eps->yMin;
    (*yMaxP) = eps->yMax;
  }

void epswr_def_report_window
  ( epswr_def_figure_t *eps,
    const char *requestor,
    const char *which,
    const char *action,
    double hxMin, double hxMax,
    double vyMin, double vyMax
  )
  { if (requestor != NULL) { fprintf(stderr, "{%s}: ", requestor); }
    fprintf(stderr, "%s window %s", which, action);
    fprintf(stderr, " [%6.1f _ %6.1f] × [%6.1f _ %6.1f]\n", hxMin, hxMax, vyMin, vyMax);
  }
