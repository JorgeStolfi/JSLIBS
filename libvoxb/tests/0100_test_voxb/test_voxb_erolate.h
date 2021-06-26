/* test_voxb_erolate.h --- tests the primitive shapes in {voxb_obj.h}. */
/* Last edited on 2021-06-22 13:45:38 by jstolfi */

#ifndef test_voxb_erolate_H
#define test_voxb_erolate_H

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <ppv_array.h>
#include <voxb_erolate.h>

void test_voxb_erolate(ppv_array_t *A, r3_t *ctr, r3_t *rad, double smr);
  /* Splats into {A} a pierced cube with center {ctr} and radius {rad}.
    Then calls {voxb_erolate_with_ball} with erosion/dilation radius {smr}. */

#endif
