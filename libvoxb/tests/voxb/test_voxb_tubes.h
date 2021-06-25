/* test_voxb_tubes.h --- tests the tube-splatting functions in {voxb_splat_tube.h}. */
/* Last edited on 2021-06-22 13:46:26 by jstolfi */

#ifndef test_voxb_tubes_H
#define test_voxb_tubes_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <ppv_array.h>

#include <r3_path.h>

void test_voxb_tubes(ppv_array_t *A, r3_t *ctr, r3_t *rad);
  /* Splats some tubes in {A}, spanning the box whith half-size {rad} 
    and center {ctr}. */

void test_voxb_tubes_helix(ppv_array_t *A, r3_t *ctr, r3_t *rad);
  /* Tests helical tubes. */

void test_voxb_tubes_segment(ppv_array_t *A, r3_t *ctr, r3_t *rad);
  /* Tests tubes defined by two {r3_path_state_t}s. */

void test_voxb_tubes_bezier(ppv_array_t *A, r3_t *ctr, r3_t *rad);
  /* Tests tubes defined by bezier control points. */

/* INTERNAL TOOLS */

void test_voxb_rescale_r3(r3_t *p, double scale, r3_t *shift);
  /* Affinely maps the point {p} from the cube {[-1 _ +1]^3} to {shift + [-scale _ +scale]^3}
    plus the given {shift} in the {X} coordinate. */

void test_voxb_rescale_r3_path_state(r3_path_state_t *P, double scale, r3_t *shift);
  /* Affinely maps the state {P} from the cube {[-1 _ +1]^3} to {shift + [-scale _ +scale]^3}
    plus the given {shift} in the {X} coordinate. */

void test_voxb_rescale_r3_motion_state(r3_motion_state_t *P, double scale, r3_t *shift);
  /* Affinely maps the state {P} from the cube {[-1 _ +1]^3} to {shift + [-scale _ +scale]^3}
    plus the given {shift} in the {X} coordinate. */

#endif
