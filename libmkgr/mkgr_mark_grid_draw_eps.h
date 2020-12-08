#ifndef mkgr_mark_grid_draw_eps_H
#define mkgr_mark_grid_draw_eps_H

/* mkgr_mark_grid_draw_eps.h - functions to draw grid of marks as EPS files. */
/* Last edited on 2020-11-29 20:33:43 by jstolfi */

#define _GNU_SOURCE
#include <r2.h>
#include <epswr.h>

#include <mkgr_mark.h>
#include <mkgr_mark_grid.h>

void mkgr_mark_grid_draw_eps
  ( epswr_figure_t *eps, 
    mkgr_mark_grid_t *gr, 
    r2_t *org,
    double dpi
  );
  /* Draws the marks described in {gr} to the Encapsulated Postscript writer
    {eps}, with all centers shifted by {org} (if not {NULL}). All coordinates and
    dimensions are assumed to be millimeters. The mark centers are rounded
    to integer multiples of printer dots, as defined by the resolution {dpi}
    (dots per inch). */

#endif
