#ifndef mkgr_mark_grid_H
#define mkgr_mark_grid_H

/* mkgr_mark_grid.h - functions to create and fraw grids of marks. */
/* Last edited on 2020-11-29 21:23:39 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <frgb.h>

#include <mkgr_mark.h>

typedef struct mkgr_mark_grid_t 
  { int32_t NM;               /* Number of marks in grid. */
    r2_t pMin;            /* Low corner of bounding box (mm). */
    r2_t pMax;            /* High corner of bounding box (mm). */
    double rbgcirc;       /* Radius of background circle, or 0 if none. */
    mkgr_mark_vec_t mark; /* List of marks. */
  } mkgr_mark_grid_t;
  /* Geometry of the chart.  The bounding box tightly includes all marks. */

mkgr_mark_grid_t *mkgr_mark_grid_new(void);
  /* Allocates a grid data structure, initially with no marks. */

void mkgr_mark_grid_free(mkgr_mark_grid_t *gr);
  /* Frees the data structure, including the marks array and header. */

void mkgr_mark_grid_append_mark(mkgr_mark_grid_t *gr, mkgr_mark_t *mk);
  /* Appends the mark {mk} to the grid {gr}.  A no-op if the mark's 
    total radius (including outline, if any) is zero. */

typedef mkgr_mark_t mkgr_mark_grid_def_proc_t(int32_t ix, int32_t iy);
  /* Type of a procedure that creates the description of 
    a markin a mark grid, given its indices {ix,iy}. */

void mkgr_mark_grid_append_marks
  ( mkgr_mark_grid_t *gr,
    i2_t iMin, 
    i2_t iMax, 
    mkgr_mark_grid_def_proc_t *defmark
  );
  /* Appends to the grid {gr} a rectangular array of marks of marks with
    indices ranging from {iMin} to {iMax}.  The parameters
    of the mark with indices {ix,iy}  are obtained
    by calling {defmark(ix,iy,...)}. Expands the tables as needed and enlarges
    the bounding box {gr->pMin,gr->pMax}. */

void mkgr_mark_grid_write(FILE *wr, mkgr_mark_grid_t *gr);
  /* Writes the data of marks in grid {gr} to file {wr}. */

#endif
