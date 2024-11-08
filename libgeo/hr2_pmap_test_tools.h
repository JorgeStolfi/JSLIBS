/* Test tools for {hr2_pmap.h} test program  */
/* Last edited on 2024-11-07 23:43:08 by stolfi */

#ifndef hr2_pmap_test_tools_H
#define hr2_pmap_test_tools_H

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <jsrandom.h>
#include <affirm.h>

#include <rn.h>
#include <r3x3.h>
#include <r3.h>
#include <r2x2.h>
#include <r2.h>

#include <hr2.h>
#include <hr2_pmap.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

void hr2_test_print_pmap(char *name, hr2_pmap_t *M);
  /* Prints {M} to stderr. */

void hr2_test_throw_pmap(hr2_pmap_t *M);
  /* Fills {M} with a random projective map. */

void hr2_test_throw_aff_map(hr2_pmap_t *M);
  /* Fills {M} with a random affine map. */

void hr2_test_choose_r2_point_pairs
  ( hr2_pmap_type_t type,
    sign_t sgn,
    bool_t tight,
    bool_t ident,
    bool_t verbose,
    int32_t *np_P,
    r2_t **p1_P,
    r2_t **p2_P,
    double **w_P,
    hr2_pmap_t *M_P
  );
  /* Chooses a suitable number {np} of point pairs {p1[i],p2[i]} and
    weights {w[i]}, with {i} in {0..np-1}, for testing the quadratic
    optimization with maps of the given {type} and {sgn}.
    The {sgn} must be {-1} or {+1}. 
    
    If {tight}, {np} will be the minimum number {nr} of points needed to
    define such a map, and that map will be exact.
    
    The procedure will choose two maps {M1,M2} of the given {type} and
    {sgn}, then generate each pair of points {p1[i],p2[i]} by mapping
    the same random point {p0[i]} through these two maps, and adding to
    both small amounts of random noise. Thus, for all {i}, {p1[i]} and
    {p2[i]} will be approximately related by the nominal map {M=M1^{-1}
    M2}.
    
    If {ident} is true, the nominal map {M} will be the identity if
    {sgn} is {+1}, or the {XY} swap map if {sgn} is {-1}. Otherwise {M}
    will be a random map of the given {type} and {sgn}.
    
    Returns {np,p1,p2,w} in {*np_P,*p1_P,*p2_P,*w_P}.  Sometimes,
    the weight table {w} will be {NULL}.
    
    Also returns in {*M_P} the nominal projective map {M=M1^{-1} M2}. */

void hr2_test_do_check_pmap
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

#define hr2_test_check_pmap(name, M, eps, msg) \
  hr2_test_do_check_pmap((name), (M), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_do_check_map_r3
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

#define hr2_test_check_map_r3(name, p, anti, twirl, row, M, q, msg) \
  hr2_test_do_check_map_r3((name), (p), (anti), (twirl), (row), (M), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_do_check_pmap_point
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

#define hr2_test_check_pmap_point(name, p, anti, twirl, M, inv, q, msg) \
  hr2_test_do_check_pmap_point((name), (p), (anti), (twirl), (M), (inv), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_do_check_pmap_line
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

#define hr2_test_check_pmap_line(name, A, anti, twirl, M, inv, B, msg) \
  hr2_test_do_check_pmap_line((name), (A), (anti), (twirl), (M), (inv), (B), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_do_check_pmap_r2_point
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
  /* Same as {hr2_test_do_check_pmap_point}, but {p} and {q} are
    Cartesian points that are converted to homogeneous in order to
    apply the projective map, and back to Cartesian in order to
    compare them. */

#define hr2_test_check_pmap_r2_point(name, p, M, inv, q, msg) \
  hr2_test_do_check_pmap_r2_point((name), (p), (M), (inv), (q), (msg), __FILE__, __LINE__, __FUNCTION__)

typedef bool_t hr2_test_pmap_check_proc_t(hr2_pmap_t *M);
  /* Tipe of a procedure that checks whether a map {M} is OK. */

void hr2_test_perturbed_pmap(hr2_pmap_t *M, hr2_test_pmap_check_proc_t *ok, r3x3_t *P, char *fname);
  /* Makes perturbations of size {P[i,j]} in turn to every element {[i,j]} of {M.dir}
    and {M.inv}, and requires that {ok} returns {FALSE} in all cases.
    If {ok} returns {TRUE} for any case, prints "{fname} failed". */

#endif
