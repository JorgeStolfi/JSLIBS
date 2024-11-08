/* See {hr3_test_tools.h}. */
/* Last edited on 2024-11-04 06:37:10 by stolfi */

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

#include <hr3.h>
#include <hr3_test_tools.h>

void hr3_test_do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func)
  { if (x != y)
      { fprintf(stderr, " ** %+20.16e %+20.16e differ\n", x, y);
        programerror(msg, file, lnum, func);
      }
  }

void hr3_test_do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %+20.16e %+20.16e", x, y);
        fprintf(stderr, " diff = %+20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

void hr3_test_do_check_hr3_point_eps
  ( char *name,
    hr3_point_t *q_cmp, 
    hr3_point_t *q_exp, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  {
    double diff = hr3_pt_pt_diff(q_cmp, q_exp);
    if (diff > eps)
      { fprintf(stderr, " ** ");
        if (name != NULL) { fprintf(stderr, "(%s) ", name); }
        fprintf(stderr, "\n");

        hr3_test_show_hr3_point("was", q_cmp);
        hr3_test_show_hr3_point("should be", q_exp);
        fprintf(stderr, " diff = %20.16e  max = %+20.16e\n", diff, eps);

        programerror(msg, file, lnum, func);
      }
  }
        
void hr3_test_show_hr3_point(char *tag, hr3_point_t *p)
  { 
    if (p != NULL)
      { /* Normalize the {hr3_point_t} to unit {\RR^4} norm: */
        hr3_point_t xx;
        r4_dir(&(p->c), &(xx.c));
        fprintf(stderr, "  %10s", tag); 
        rn_gen_print(stderr, NH, xx.c.c, "%+11.8f", "[ ", " ", " ]");
        if (xx.c.c[0] != 0.0) 
          { /* Print the Cartesian coordinates: */
            r3_t p = r3_from_hr3(&xx); r3_gen_print(stderr, &p, "%+11.8f", " = ( ", " ", " )");
          }
        fprintf(stderr, "\n");
      }
  }    

void hr3_test_show_hr3_line(char *tag, hr3_line_t *A)
  { 
    if (A != NULL)
      { /* Normalize the Grassmann coordinates to unit {\RR^4} norm: */
        hr3_line_t AA;
        (void)r6_dir(&(A->k), &(AA.k));
        fprintf(stderr, "  %10s", tag); 
        rn_gen_print(stderr, NG, AA.k.c, "%+11.8f", "{ ", " ", " }");
        fprintf(stderr, "\n");
      }
  }

void hr3_test_show_hr3_plane(char *tag, hr3_plane_t *A)
  { 
    if (A != NULL)
      { /* Normalize the homogeneous coefficients to unit {\RR^4} norm: */
        hr3_plane_t AA;
        (void)r4_dir(&(A->f), &(AA.f));
        fprintf(stderr, "  %10s", tag); 
        rn_gen_print(stderr, NG, AA.f.c, "%+11.8f", "< ", " ", " >");
        fprintf(stderr, "\n");
      }
  } 

void hr3_test_do_check_num_eps
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

void hr3_test_do_check_r3_eps
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
    double diff = r3_dist(b_exp, b_cmp);
    if (diff > eps)
      { fprintf(stderr, " ** %s: ", name);
        rn_gen_print(stderr, NC, b_cmp->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " should be "); 
        rn_gen_print(stderr, NC, b_exp->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " sdiff = %20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

void hr3_test_do_check_r4_norm_eps
  ( char *name, 
    r4_t *b_cmp, 
    r4_t *b_exp,
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  {
    r4_t bn_cmp; (void)r4_dir(b_cmp, &bn_cmp);
    r4_t bn_exp; (void)r4_dir(b_exp, &bn_exp);
    
    double diff = r4_dist(&bn_exp, &bn_cmp);
    if (diff > eps)
      { fprintf(stderr, " ** %s: ", name);
        rn_gen_print(stderr, NH, b_cmp->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " should be "); 
        rn_gen_print(stderr, NH, b_exp->c, "%+11.8f", "[ ", " ", " ]");
        fprintf(stderr, " sdiff = %20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

