/* test_voxm_tubes.h --- tests the tube-splatting functions in {voxm_splat_tube.h}. */
/* Last edited on 2021-06-12 10:43:03 by jstolfi */

#ifndef test_voxm_tubes_H
#define test_voxm_tubes_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <r3_motion.h>
#include <ppv_array.h>

#include <r3_path.h>

void test_voxm_tubes(ppv_array_desc_t *A, r3_t *ctr, r3_t *rad, double fuzzR);
  /* Splats some tubes in {A}, spanning the box whith half-size {rad} 
    and center {ctr}. */

void test_voxm_tubes_helix(ppv_array_desc_t *A, r3_t *ctr, r3_t *rad, double fuzzR);
  /* Tests helical tubes. */

void test_voxm_tubes_segment(ppv_array_desc_t *A, r3_t *ctr, r3_t *rad, double fuzzR);
  /* Tests tubes defined by two {r3_path_state_t}s. */

void test_voxm_tubes_bezier(ppv_array_desc_t *A, r3_t *ctr, r3_t *rad, double fuzzR);
  /* Tests tubes defined by bezier control points. */

/* INTERNAL TOOLS */

void test_voxm_rescale_r3(r3_t *p, double scale, r3_t *shift);
  /* Affinely maps the point {p} from the cube {[-1 _ +1]^3} to {shift + [-scale _ +scale]^3}
    plus the given {shift} in the {X} coordinate. */

void test_voxm_rescale_r3_path_state(r3_path_state_t *P, double scale, r3_t *shift);
  /* Affinely maps the state {P} from the cube {[-1 _ +1]^3} to {shift + [-scale _ +scale]^3}
    plus the given {shift} in the {X} coordinate. */

void test_voxm_rescale_r3_motion_state(r3_motion_state_t *P, double scale, r3_t *shift);
  /* Affinely maps the state {P} from the cube {[-1 _ +1]^3} to {shift + [-scale _ +scale]^3}
    plus the given {shift} in the {X} coordinate. */

#endif
