#ifndef mkgr_mark_grid_parms_H
#define mkgr_mark_grid__parmsH

/* mkgr_mark_grid_parms.h - paramters of marks in a complex grid. */
/* Last edited on 2020-11-30 13:58:42 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
    
typedef struct mkgr_mark_grid_parms_t 
  { double rcross;   /* Radius for cross marks (mm). */
    double rcircle;  /* Radius for circle marks (mm). */
    double rdot;     /* Radius for dot marks (mm).*/
    double rbgcirc;  /* Radius of background circle. */
    double lwd;      /* Line width (mm). */
    double spacing;  /* Spacing between marks (mm). */
  } mkgr_mark_grid_parms_t;
  /* A record with various grid parameters. */

void mkgr_mark_grid_parms_show(FILE *wr, char *kind, mkgr_mark_grid_parms_t *parms);
  /* Prints the grid parameters {*parms} to file {wr}, with {kind} in the title. */

void mkgr_mark_grid_parms_round(mkgr_mark_grid_parms_t *parms, double dots, double unit);
  /* Rounds the size parameters in {*parms} to multiples of prineter pixels,
    assuming the resolution {dots/unit} pixels per mm.
    The parameters {parms.lwd} and {parms.spacing} are rounded to an even number
    of pixels. */

#endif
