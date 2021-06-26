/* test_voxb_mark.h --- marks box extent for visualization. */
/* Last edited on 2021-06-22 13:45:52 by jstolfi */

#ifndef test_voxb_mark_H
#define test_voxb_mark_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <ppv_array.h>

void test_voxb_mark_corners(ppv_array_t *A, r3_t *ctr, r3_t *rad);
  /* Splats eight balls at the corners of the box with
    center {ctr} and half-side {rad}, for visualiztion. */

void test_voxb_mark_edges(ppv_array_t *A, r3_t *ctr, r3_t *rad, int ax);
  /* Splats four rods paralel to coordinate axis {ax} along the edges of the box with
    center {ctr} and half-side {rad}, for visualiztion.
    The edges */

#endif
