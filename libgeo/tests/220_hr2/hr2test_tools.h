/* hr2test_tools.h --- test tools for program for hr2_test.c  */
/* Last edited on 2024-08-29 22:15:11 by stolfi */

#ifndef hr2test_tools_H
#define hr2test_tools_H

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

#include <hr2.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

void h2tt_print_pmap(char *name, hr2_pmap_t *M);
  /* Prints {M} to stderr. */

void h2tt_throw_pmap(hr2_pmap_t *M);
  /* Fills {M} with a random projective map. */

void h2tt_throw_aff_map(hr2_pmap_t *M);
  /* Fills {M} with a random affine map. */

void h2tt_do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ, prints them, prints {msg}, and stops. */

#define check_eq(x,y,msg) \
  h2tt_do_check_eq((x), (y), (msg), __FILE__, __LINE__, __FUNCTION__)
  
void h2tt_do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ by more than {eps}, prints them, prints {msg}, and stops. */

#define h2tt_check_eps(x, y, eps, msg) \
  h2tt_do_check_eps((x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_num_eps
  ( char *name,
    double x,
    double y,
    double eps,
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If {x} and {y} differ by more than {eps}, prints {name}, {x}, {y}, and {msg}, and stops. */

#define h2tt_check_num_eps(name, x, y, eps, msg) \
  h2tt_do_check_num_eps((name), (x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_r2_eps
  ( char *name, 
    r2_t *b_cmp, 
    r2_t *b_exp,
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If points {b_cmp} and {b_exp} differ by more than {eps} prints
    {name}, {b_cmp}, {b_eps}, and {msg}, and stops.
    Otherwise does nothing. */

#define h2tt_check_r2_eps(name, b_cmp, b_exp, eps, msg) \
  h2tt_do_check_r2_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_r3_norm_eps
  ( char *name, 
    r3_t *b_cmp, 
    r3_t *b_exp,
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If the vectors {b_cmp} and {b_exp} differ by
    more than {eps} apart from a positive scaling factor, prints {name},
    {b_cmp}, {b_eps}, and {msg}, and stops. Otherwise does nothing. */ 

#define h2tt_check_r3_norm_eps(name, b_cmp, b_exp, eps, msg) \
  h2tt_do_check_r3_norm_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_hr2_point_eps
  ( char *name, 
    hr2_point_t *b_cmp, 
    hr2_point_t *b_exp, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If points {b_cmp} and {b_exp} differ by more than {eps} 
    (apart from homogeneous scaling), prints {name},
    {b_cmp}, {b_exp}, {msg}, and stops. Otherwise does nothing. */

#define h2tt_check_hr2_point_eps(name, b_cmp, b_exp, eps, msg) \
  h2tt_do_check_hr2_point_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_pmap
  ( char *name,
    hr2_pmap_t *M, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* Checks the validity of the map {M}, namely if {M.dir} and {M.inv} are 
    close inverses of each other, within tolerance {eps}. */

#define h2tt_check_pmap(name, M, eps, msg) \
  h2tt_do_check_pmap((name), (M), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_map_r3
  ( char *name, 
    r3_t *p, 
    bool_t anti, 
    bool_t twirl, 
    bool_t row,
    r3x3_t *M, 
    r3_t *q,
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If {anti} and {twirl} are false, computes a point {pM} by mapping
    {p} through {M}.
    
    If {anti} is true, tries reversing the signs of the three
    coordinates of {p} simultaneously, and proceeds as above.
    
    If {twirl} is true, tries reversing the signs of the three
    coordinates of {p} in all combinations, and proceeds as above.
    
    The parameters {anti} and {twirl} are mutually exclusive. If either
    is true, takes {pM} to be the mapped point that is closest to {q},
    apart from a positive scaling factor.
    
    Either way, the image of a point {u} by {M} is {u*M} if {row} true,
    or {M*u} if {row} is false.
    
    In any case, returns silently if the chosen {pM} is equal to {q}
    (modulo positive homogeneous scaling and a rounding error
    tolerance), and aborts with error {msg} if not. */ 

#define h2tt_check_map_r3(name, p, anti, twirl, row, M, q, msg) \
  h2tt_do_check_map_r3((name), (p), (anti), (twirl), (row), (M), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_pmap_point
  ( char *name, 
    hr2_point_t *p, 
    bool_t anti,
    bool_t twirl,
    hr2_pmap_t *M, 
    bool_t inv,
    hr2_point_t *q, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If {anti} and {twirl} are FALSE, maps the point {p} by the given projective map
    and compares it with {q}. 
    
    If {anti} is true, tries reversing the signs of the three
    coordinates of {p} simultaneously, and proceeds as above.
    
    If {twirl} is true, tries reversing the
    signs of the three coordinates of {p} in all combinations, and
    proceeds as above. 
    
    The parameters {anti} and {twirl} are mutually exclusive. If either
    is true, takes {pM} to be the mapped point that is closest to {q}
    in the sperical model sense.
     
    In any case, returns silently if the chosen {pM} matches {q}
    (modulo rounding errors); aborts with error {msg} if all attempts
    failed. The projective map used is {M} if {inv} is false, or its
    inverse if {inv} is true. */

#define h2tt_check_pmap_point(name, p, anti, twirl, M, inv, q, msg) \
  h2tt_do_check_pmap_point((name), (p), (anti), (twirl), (M), (inv), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_pmap_line
  ( char *name, 
    hr2_line_t *A, 
    bool_t anti,
    bool_t twirl,
    hr2_pmap_t *M,
    bool_t inv,
    hr2_line_t *B, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If {twirl} is FALSE, maps the line {A} by the specified projective
    map and compares it with {B}. 
    
    If {anti} is true, tries reversing the signs of the three
    coefficients of {A} simultaneously, and proceeds as above.
    
    If {twirl} is true, tries reversing the signs of the three
    coefficients of {A} in all combinations, and proceeds as above.
    
    In either case, returns silently if some attempt produced a match
    (modulo rounding errors); aborts with error {msg} if all attempts
    failed. The projective map used is {M} if {inv} is false, or its
    inverse if {inv} is true. */

#define h2tt_check_pmap_line(name, A, anti, twirl, M, inv, B, msg) \
  h2tt_do_check_pmap_line((name), (A), (anti), (twirl), (M), (inv), (B), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_do_check_pmap_r2_point
  ( char *name, 
    r2_t *p, 
    hr2_pmap_t *M, 
    bool_t inv,
    r2_t *q,
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* Same as {h2tt_do_check_pmap_point}, but {p} and {q} are
    Cartesian points that are converted to homogeneous in order to
    apply the projective map, and back to Cartesian in order to
    compare them. */

#define h2tt_check_pmap_r2_point(name, p, M, inv, q, msg) \
  h2tt_do_check_pmap_r2_point((name), (p), (M), (inv), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

void h2tt_show_hr2_point(char *tag, hr2_point_t *p);
  /* If {p} is not {NULL}, prints the homogeneous coordinates of point {p} on a single line,
    prefixed by {tag}. If finite, also prints the Cartesian
    coordinates. If {p} is {NULL}, does nothing. */

void h2tt_show_hr2_line(char *tag, hr2_line_t *A);
  /* If {A} is not {NULL}, prints the homogeneous coefficients line {A} on a single line,
    prefixed by {tag}.  If {A} is {NULL}, does nothing. */

#endif
