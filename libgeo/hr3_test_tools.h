/* Test tools for {hr3.h} test programs */
/* Last edited on 2024-09-17 17:34:31 by stolfi */

#ifndef hr3_test_tools_H
#define hr3_test_tools_H

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

#define NH 4
  /* Number of homogeneous coords in a point and coeffs in a plane. */

#define NG 6
  /* Number of Grassmann coordinates in a line. */

#define NC 3
  /* Number of Cartesian coordinates in a point. */

void hr3_test_do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ, prints them, prints {msg}, and stops. */

#define check_eq(x,y,msg) \
  hr3_test_do_check_eq((x), (y), (msg), __FILE__, __LINE__, __FUNCTION__)
  
void hr3_test_do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ by more than {eps}, prints them, prints {msg}, and stops. */

#define hr3_test_check_eps(x, y, eps, msg) \
  hr3_test_do_check_eps((x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr3_test_do_check_num_eps
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

#define hr3_test_check_num_eps(name, x, y, eps, msg) \
  hr3_test_do_check_num_eps((name), (x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr3_test_do_check_r3_eps
  ( char *name, 
    r3_t *b_cmp, 
    r3_t *b_exp,
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If points {b_cmp} and {b_exp} differ by more than {eps} prints
    {name}, {b_cmp}, {b_eps}, and {msg}, and stops.
    Otherwise does nothing. */

#define hr3_test_check_r3_eps(name, b_cmp, b_exp, eps, msg) \
  hr3_test_do_check_r3_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr3_test_do_check_r4_norm_eps
  ( char *name, 
    r4_t *b_cmp, 
    r4_t *b_exp,
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If the vectors {b_cmp} and {b_exp} differ by
    more than {eps} apart from a positive scaling factor, prints {name},
    {b_cmp}, {b_eps}, and {msg}, and stops. Otherwise does nothing. */ 

#define hr3_test_check_r4_norm_eps(name, b_cmp, b_exp, eps, msg) \
  hr3_test_do_check_r4_norm_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr3_test_do_check_hr3_point_eps
  ( char *name, 
    hr3_point_t *b_cmp, 
    hr3_point_t *b_exp, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  );
  /* If points {b_cmp} and {b_exp} differ by more than {eps} 
    (apart from homogeneous scaling), prints {name},
    {b_cmp}, {b_exp}, {msg}, and stops. Otherwise does nothing. */

#define hr3_test_check_hr3_point_eps(name, b_cmp, b_exp, eps, msg) \
  hr3_test_do_check_hr3_point_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr3_test_show_hr3_point(char *tag, hr3_point_t *p);
  /* If {p} is not {NULL}, prints the homogeneous coordinates of point {p} on a single line,
    prefixed by {tag}. If finite, also prints the Cartesian
    coordinates. If {p} is {NULL}, does nothing. */

void hr3_test_show_hr3_line(char *tag, hr3_line_t *L);
  /* If {A} is not {NULL}, prints the Grassmann coefficients line {L} on a single line,
    prefixed by {tag}.  If {L} is {NULL}, does nothing. */

void hr3_test_show_hr3_plane(char *tag, hr3_plane_t *A);
  /* If {A} is not {NULL}, prints the homogeneous coefficients of plane {A} on a single line,
    prefixed by {tag}.  If {A} is {NULL}, does nothing. */

#endif
