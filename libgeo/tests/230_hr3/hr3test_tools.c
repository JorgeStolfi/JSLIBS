/* See {hr3test_tools.h}. */
/* Last edited on 2023-10-09 21:28:56 by stolfi */

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
#include <hr3test_tools.h>

void h3tt_do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func)
  { if (x != y)
      { fprintf(stderr, " ** %+20.16e %+20.16e differ\n", x, y);
        programerror(msg, file, lnum, func);
      }
  }

void h3tt_do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %+20.16e %+20.16e", x, y);
        fprintf(stderr, " diff = %+20.16e  max = %+20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

void h3tt_do_check_hr3_point_eps
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

        h3tt_show_hr3_point("was", q_cmp);
        h3tt_show_hr3_point("should be", q_exp);
        fprintf(stderr, " diff = %20.16e  max = %+20.16e\n", diff, eps);

        programerror(msg, file, lnum, func);
      }
  }
        
void h3tt_show_hr3_point(char *tag, hr3_point_t *p)
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

void h3tt_show_hr3_line(char *tag, hr3_line_t *A)
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

void h3tt_show_hr3_plane(char *tag, hr3_plane_t *A)
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

void h3tt_throw_pmap(hr3_pmap_t *M)
  {
    r4_t a;
    for (int32_t i = 0; i < NH; i++)
      { r4_throw_cube(&a);
        for (int32_t j = 0; j < NH; j++) { M->dir.c[i][j] = a.c[j]; }
      }
    r4x4_inv(&(M->dir), &(M->inv));
  }

void h3tt_throw_aff_map(hr3_pmap_t *M)
  {
    for (int32_t i = 0; i < NH; i++)
      { r3_t p;
        r3_throw_cube(&p);
        M->dir.c[i][0] = (i == 0 ? 1.0 : 0.0);
        M->dir.c[i][1] = p.c[0];
        M->dir.c[i][2] = p.c[1];
        M->dir.c[i][3] = p.c[2];
      }
    r4x4_inv(&(M->dir), &(M->inv));
  }
  
void h3tt_print_pmap(char *name, hr3_pmap_t *M)
  { fprintf(stderr, "%s = ", name);
    hr3_pmap_gen_print(stderr, M, "%+10.6f", NULL,  NULL,NULL,NULL, NULL,NULL,NULL, "\n");
  }

void h3tt_do_check_r4_map
  ( char *name, 
    r4_t *p,
    bool_t flip, 
    bool_t row, 
    r4x4_t *M, 
    r4_t *q,
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  {
    r4_t qn; (void)r4_dir(q, &qn);
    r4_t pM;
    if (flip)
      { /* {b_exp} stands for any sign-flipped version of itself: */
        double d2min = +INF;
        for (int32_t sz = -1; sz <= +1; sz += 2)
          for (int32_t sy = -1; sy <= +1; sy += 2)
            for (int32_t sx = -1; sx <= +1; sx += 2)
              for (int32_t sw = -1; sw <= +1; sw += 2)
                { r4_t u = (r4_t){{sw*p->c[0], sx*p->c[1], sy*p->c[2], sz*p->c[3]}};
                  r4_t uM;
                  if (row) 
                    { r4x4_map_row(&u, M, &pM); }
                  else
                    { r4x4_map_col(M, &u, &pM); }
                  (void)r4_dir(&uM, &uM);
                  double d2 = r4_dist_sqr(&uM, &qn);
                  if (d2 < d2min) { d2min = d2; pM = uM; }
                }
      }
    else
      { if (row) 
          { r4x4_map_row(p, M, &pM); }
        else
          { r4x4_map_col(M, p, &pM); }
        (void)r4_dir(&pM, &pM);
      }
    h3tt_do_check_r4_norm_eps(name, &pM, &qn, 0.00000001, msg, file, lnum, func);
  }  

void h3tt_do_check_pmap_point
  ( char *name, 
    hr3_point_t *p, 
    bool_t flip,
    hr3_pmap_t *M, 
    bool_t inv,
    hr3_point_t *q, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  { 
    r4x4_t *N = (inv ? &(M->inv) : &(M->dir));
    h3tt_do_check_r4_map(name, &(p->c), flip, TRUE, N, &(q->c), msg, file, lnum, func);
  }

void h3tt_do_check_pmap_plane
  ( char *name, 
    hr3_plane_t *A, 
    bool_t flip,
    hr3_pmap_t *M, 
    bool_t inv, 
    hr3_plane_t *B, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  {
    r4x4_t *N = (inv ? &(M->dir) : &(M->inv));
    h3tt_do_check_r4_map(name, &(A->f), flip, FALSE, N, &(B->f), msg, file, lnum, func);
  }
  
void h3tt_do_check_pmap_r3_point
  ( char *name, 
    r3_t *p, 
    hr3_pmap_t *M, 
    bool_t inv,
    r3_t *q,
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  { 
    r3_t pM;
    if (inv)
      { pM = hr3_pmap_inv_r3_point(p, M); }
    else
      { pM = hr3_pmap_r3_point(p, M); }
    h3tt_do_check_r3_eps(name, &pM, q, 0.00000001, msg, file, lnum, func);
  }

void h3tt_do_check_num_eps
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

void h3tt_do_check_r3_eps
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

void h3tt_do_check_r4_norm_eps
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

