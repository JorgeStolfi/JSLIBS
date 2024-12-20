/* Test tools for {hr2.h} test programs  */
/* Last edited on 2024-12-05 10:27:17 by stolfi */

#ifndef hr2_test_tools_H
#define hr2_test_tools_H

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

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

void hr2_test_tools_do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ, prints them, prints {msg}, and stops. */

#define hr2_test_tools_check_eq(x,y,msg) \
  hr2_test_tools_do_check_eq((x), (y), (msg), __FILE__, __LINE__, __FUNCTION__)
  
void hr2_test_tools_do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func);
  /* If {x} and {y} differ by more than {eps}, prints them, prints {msg}, and stops. */

#define hr2_test_tools_check_eps(x, y, eps, msg) \
  hr2_test_tools_do_check_eps((x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_tools_do_check_num_eps
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

#define hr2_test_tools_check_num_eps(name, x, y, eps, msg) \
  hr2_test_tools_do_check_num_eps((name), (x), (y), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_tools_do_check_r2_eps
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

#define hr2_test_tools_check_r2_eps(name, b_cmp, b_exp, eps, msg) \
  hr2_test_tools_do_check_r2_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_tools_do_check_r3_norm_eps
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

#define hr2_test_tools_check_r3_norm_eps(name, b_cmp, b_exp, eps, msg) \
  hr2_test_tools_do_check_r3_norm_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_tools_do_check_hr2_point_eps
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

#define hr2_test_tools_check_hr2_point_eps(name, b_cmp, b_exp, eps, msg) \
  hr2_test_tools_do_check_hr2_point_eps((name), (b_cmp), (b_exp), (eps), (msg), __FILE__, __LINE__, __FUNCTION__)

void hr2_test_tools_show_hr2_point(char *tag, hr2_point_t *p);
  /* If {p} is not {NULL}, prints the homogeneous coordinates of point {p} on a single line,
    prefixed by {tag}. If finite, also prints the Cartesian
    coordinates. If {p} is {NULL}, does nothing. */

void hr2_test_tools_show_hr2_line(char *tag, hr2_line_t *A);
  /* If {A} is not {NULL}, prints the homogeneous coefficients line {A} on a single line,
    prefixed by {tag}.  If {A} is {NULL}, does nothing. */

#endif
