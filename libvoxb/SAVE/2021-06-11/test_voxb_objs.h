/* test_voxb_objs.h --- tests the primitive shapes in {voxb_obj.h}. */
/* Last edited on 2021-06-11 02:14:17 by jstolfi */

#ifndef test_voxb_objs_H
#define test_voxb_objs_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <ppv_array.h>

void test_voxb_objs(ppv_array_t *a, r3_t *ctr, r3_t *rad);
  /* Splats some primitive objects onto {a}, spanning the box whith 
    half-size {rad} and center {ctr} */

void test_voxb_objs_ball(ppv_array_t *a, r3_t *ctr, double rad);
  /* Splats a ball near the +X side of the array. */

void test_voxb_objs_donut(ppv_array_t *a, r3_t *ctr, double rad);
  /* Splats a donut at the center, in natural orientation. */

void test_voxb_objs_rod(ppv_array_t *a, r3_t *ctr, double rad);
  /* Splats a rod through the donut hole. */
 
void test_voxb_objs_tube(ppv_array_t *a, r3_t *ctr, double rad);
  /* Splats a tube adjacent to the donut. */

void test_voxb_objs_cube_hole(ppv_array_t *a, r3_t *ctr, double rad);
  /* Splats a tilted cube above the center with a drilled cubical hole. */
 
void test_voxb_objs_box(ppv_array_t *a, r3_t *ctr, double rad);
  /* Splats a box below center. */
 
void test_voxb_objs_rounded_box(ppv_array_t *a, r3_t *ctr, double rad);
  /* Splats a box below center. */

void test_voxb_objs_cup(ppv_array_t *a, r3_t *ctr, double rad);
  /* Splats a cylindrical cup near the corner of an array. */

#endif
