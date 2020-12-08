/* test_voxm_mark.h --- marks box extent for visualization. */
/* Last edited on 2016-04-03 14:25:16 by stolfilocal */

#ifndef test_voxm_mark_H
#define test_voxm_mark_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <ppv_array.h>

void test_voxm_mark_corners(ppv_array_t *a, r3_t *ctr, r3_t *rad, double fuzzR);
  /* Splats eight balls at the corners of the box with
    center {ctr} and half-side {rad}, for visualiztion. */

void test_voxm_mark_edges(ppv_array_t *a, r3_t *ctr, r3_t *rad, int ax, double fuzzR);
  /* Splats four rods paralel to coordinate axis {ax} along the edges of the box with
    center {ctr} and half-side {rad}, for visualiztion.
    The edges */

#endif
