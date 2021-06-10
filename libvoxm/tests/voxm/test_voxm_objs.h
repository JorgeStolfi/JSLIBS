/* test_voxm_obj.h --- tests the primitive shapes in {voxm_obj.h}. */
/* Last edited on 2021-06-09 16:39:34 by jstolfi */

#ifndef test_voxm_obj_H
#define test_voxm_obj_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <ppv_array.h>

void test_voxm_objs(ppv_array_t *a, r3_t *ctr, r3_t *rad, double fuzzR);
  /* Splats some primitive objects onto {a}, spanning the box whith 
    half-size {rad} and center {ctr} */

void test_voxm_objs_ball(ppv_array_t *a, r3_t *ctr, double rad, double fuzzR);
  /* Splats a ball near the +X side of the array. */

void test_voxm_objs_donut(ppv_array_t *a, r3_t *ctr, double rad, double fuzzR);
  /* Splats a donut at the center, in natural orientation. */

void test_voxm_objs_rod(ppv_array_t *a, r3_t *ctr, double rad, double fuzzR);
  /* Splats a rod through the donut hole. */
 
void test_voxm_objs_tube(ppv_array_t *a, r3_t *ctr, double rad, double fuzzR);
  /* Splats a tube adjacent to the donut. */

void test_voxm_objs_cube_hole(ppv_array_t *a, r3_t *ctr, double rad, double fuzzR);
  /* Splats a tilted cube above the center with a drilled cubical hole. */
 
void test_voxm_objs_box(ppv_array_t *a, r3_t *ctr, double rad, double fuzzR);
  /* Splats a box below center. */
 
void test_voxm_objs_rounded_box(ppv_array_t *a, r3_t *ctr, double rad, double fuzzR);
  /* Splats a box below center. */

void test_voxm_objs_cup(ppv_array_t *a, r3_t *ctr, double rad, double fuzzR);
  /* Splats a cylindrical cup near the corner of an array. */

#endif
