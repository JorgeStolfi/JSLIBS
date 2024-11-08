/* See {hr2_test_tools.h}. */
/* Last edited on 2024-11-04 06:36:53 by stolfi */

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

#include <hr2_test_tools.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

void hr2_test_do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func)
  { if (x != y)
      { fprintf(stderr, " ** %+20.16e %+20.16e differ\n", x, y);
        programerror(msg, file, lnum, func);
      }
  }

void hr2_test_do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %+20.16e %+20.16e", x, y);
        fprintf(stderr, " diff = %+20.16e  max = %20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

void hr2_test_do_check_hr2_point_eps
  ( char *name,
    hr2_point_t *q_cmp, 
    hr2_point_t *q_exp, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  {
    double diff = hr2_pt_pt_diff(q_cmp, q_exp);
    if (diff > eps)
      { fprintf(stderr, " ** ");
        if (name != NULL) { fprintf(stderr, "(%s) ", name); }
        fprintf(stderr, "\n");

        hr2_test_show_hr2_point("was", q_cmp);
        hr2_test_show_hr2_point("should be", q_exp);
        fprintf(stderr, " diff = %20.16e  max = %20.16e\n", diff, eps);

        programerror(msg, file, lnum, func);
      }
  }
        
void hr2_test_show_hr2_point(char *tag, hr2_point_t *p)
  { 
    if (p != NULL)
      { /* Normalize the {hr2_point_t} to unit {\RR^3} norm: */
        hr2_point_t xx;
        r3_dir(&(p->c), &(xx.c));
        fprintf(stderr, "  %10s", tag); 
        rn_gen_print(stderr, NH, xx.c.c, "%+11.8f", "[ ", " ", " ]");
        if (xx.c.c[0] != 0.0) 
          { /* Print the Cartesian coordinates: */
            r2_t p = r2_from_hr2(&xx); r2_gen_print(stderr, &p, "%+11.8f", " = ( ", " ", " )");
          }
        fprintf(stderr, "\n");
      }
  }    

void hr2_test_show_hr2_line(char *tag, hr2_line_t *A)
  { 
    if (A != NULL)
      { /* Normalize the {hr2_point_t} to unit {\RR^3} norm: */
        hr2_line_t AA;
        (void)r3_dir(&(A->f), &(AA.f));
        fprintf(stderr, "  %10s", tag); 
        rn_gen_print(stderr, NH, AA.f.c, "%+11.8f", "< ", " ", " >");
        fprintf(stderr, "\n");
      }
  }    

void hr2_test_do_check_num_eps
  ( char *name,
    double x,
    double y,
    double eps,
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %s: %+20.16e %+20.16e", name, x, y);
        fprintf(stderr, " diff = %+20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

void hr2_test_do_check_r2_eps
  ( char *name, 
    r2_t *b_cmp, 
    r2_t *b_exp, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  {
    double diff = r2_dist(b_exp, b_cmp);
    if (diff > eps)
      { fprintf(stderr, " ** %s: ", name);
        rn_gen_print(stderr, NC, b_cmp->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " should be "); 
        rn_gen_print(stderr, NC, b_exp->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " sdiff = %20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

void hr2_test_do_check_r3_norm_eps
  ( char *name, 
    r3_t *b_cmp, 
    r3_t *b_exp,
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  {
    r3_t bn_cmp; (void)r3_dir(b_cmp, &bn_cmp);
    r3_t bn_exp; (void)r3_dir(b_exp, &bn_exp);
    
    double diff = r3_dist(&bn_exp, &bn_cmp);
    if (diff > eps)
      { fprintf(stderr, " ** %s: ", name);
        rn_gen_print(stderr, NH, b_cmp->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " should be "); 
        rn_gen_print(stderr, NH, b_exp->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " sdiff = %20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

