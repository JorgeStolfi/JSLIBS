/* Test tools for {hr3_pmap.h} test program  */
/* Last edited on 2024-09-17 17:21:43 by stolfi */

#ifndef hr3_pmap_test_tools_H
#define hr3_pmap_test_tools_H

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <flt.h>
#include <jsrandom.h>
#include <affirm.h>

#include <rn.h>
#include <r3x3.h>
#include <r3.h>
#include <r2x2.h>
#include <r2.h>
#include <hr3.h>

#include <hr3_pmap.h>

#define NH 4
  /* Number of homogeneous coords in a point and coeffs in a plane. */

#define NG 6
  /* Number of Grassmann coordinates in a line. */

#define NC 3
  /* Number of Cartesian coordinates in a point. */

void hr3_test_print_pmap(char *name, hr3_pmap_t *M);
  /* Prints {M} to stderr. */

void hr3_test_throw_pmap(hr3_pmap_t *M);
  /* Fills {M} with a random projective map. */

void hr3_test_throw_aff_map(hr3_pmap_t *M);
  /* Fills {M} with a random affine map. */

void hr3_test_do_check_r4_map
  ( char *name, 
    r4_t *p, 
    bool_t twirl, 
    bool_t row,
    r4x4_t *M, 
    r4_t *q,
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If {twirl} is false, computes a point {pM} by mapping {p} through
    {M}. If {twirl} is true, considers all points that differ from {p} by
    the signs of any or all coordinates, maps each through {M}, and lets
    {pM} be the one that is closest to {q}, apart from a positive
    scaling factor.
    
    Either way, the image of a point {u} by {M} is {u*M} if {row} true,
    or {M*u} if {row} is false.
    
    In any case, returns silently if the final {pM} is equal to {q}
    (modulo positive homogeneous scaling and a rounding error
    tolerance), and aborts with error {msg} if not. */ 

#define hr3_test_check_r4_map(name, p, twirl, row, M, q, msg) \
  hr3_test_do_check_r4_map((name), (p), (twirl), (row), (M), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr3_test_do_check_pmap_point
  ( char *name, 
    hr3_point_t *p, 
    bool_t twirl,
    hr3_pmap_t *M, 
    bool_t inv,
    hr3_point_t *q, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If {twirl} is FALSE, maps the point {p} by the given projective map
    and compares it with {q}. If {twirl} is true, tries reversing the
    signs of the three coordinates of {p} in all combinations, and
    proceeds as above. In either case, returns silently if some attempt
    produced a match (modulo rounding errors); aborts with error {msg}
    if all attempts failed. The projective map used is {M} if {inv} is
    false, or its inverse if {inv} is true. */

#define hr3_test_check_pmap_point(name, p, twirl, M, inv, q, msg) \
  hr3_test_do_check_pmap_point((name), (p), (twirl), (M), (inv), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr3_test_do_check_pmap_plane
  ( char *name, 
    hr3_plane_t *A, 
    bool_t twirl,
    hr3_pmap_t *M,
    bool_t inv,
    hr3_plane_t *B, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If {twirl} is FALSE, maps the line {A} by the specified projective
    map and compares it with {B}. If {twirl} is true, tries reversing the
    signs of the three coordinates of {A} in all combinations, and
    proceeds as above. In either case, returns silently if some attempt
    produced a match (modulo rounding errors); aborts with error {msg}
    if all attempts failed. The projective map used is {M} if {inv} is
    false, or its inverse if {inv} is true. */

#define hr3_test_check_pmap_plane(name, A, twirl, M, inv, B, msg) \
  hr3_test_do_check_pmap_plane((name), (A), (twirl), (M), (inv), (B), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr3_test_do_check_pmap_r3_point
  ( char *name, 
    r3_t *p, 
    hr3_pmap_t *M, 
    bool_t inv,
    r3_t *q,
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* Same as {hr3_test_do_check_pmap_point}, but {p} and {q} are
    Cartesian points that are converted to/from homogeneous 
    in order to apply the projective map. */

#define hr3_test_check_pmap_r3_point(name, p, M, inv, q, msg) \
  hr3_test_do_check_pmap_r3_point((name), (p), (M), (inv), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

#endif
