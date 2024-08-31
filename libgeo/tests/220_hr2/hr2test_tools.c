/* See {hr2test_tools.h}. */
/* Last edited on 2024-08-30 04:25:30 by stolfi */

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
#include <hr2test_tools.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

void h2tt_do_check_eq(double x, double y, char *msg, char *file, int32_t lnum, const char *func)
  { if (x != y)
      { fprintf(stderr, " ** %+20.16e %+20.16e differ\n", x, y);
        programerror(msg, file, lnum, func);
      }
  }

void h2tt_do_check_eps(double x, double y, double eps, char *msg, char *file, int32_t lnum, const char *func)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %+20.16e %+20.16e", x, y);
        fprintf(stderr, " diff = %+20.16e  max = %20.16e\n", diff, eps);
        programerror(msg, file, lnum, func);
      }
  }

void h2tt_do_check_hr2_point_eps
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

        h2tt_show_hr2_point("was", q_cmp);
        h2tt_show_hr2_point("should be", q_exp);
        fprintf(stderr, " diff = %20.16e  max = %20.16e\n", diff, eps);

        programerror(msg, file, lnum, func);
      }
  }
        
void h2tt_show_hr2_point(char *tag, hr2_point_t *p)
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

void h2tt_show_hr2_line(char *tag, hr2_line_t *A)
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

void h2tt_throw_pmap(hr2_pmap_t *M)
  {
    r3_t a;
    for (int32_t i = 0; i < NH; i++)
      { r3_throw_cube(&a);
        for (int32_t j = 0; j < NH; j++) { M->dir.c[i][j] = a.c[j]; }
      }
    r3x3_inv(&(M->dir), &(M->inv));
  }

void h2tt_throw_aff_map(hr2_pmap_t *M)
  {
    for (int32_t i = 0; i < NH; i++)
      { r2_t p;
        r2_throw_cube(&p);
        M->dir.c[i][0] = (i == 0 ? 1.0 : 0.0);
        M->dir.c[i][1] = p.c[0];
        M->dir.c[i][2] = p.c[1];
      }
    r3x3_inv(&(M->dir), &(M->inv));
  }
  
void h2tt_print_pmap(char *name, hr2_pmap_t *M)
  { fprintf(stderr, "%s = ", name);
    hr2_pmap_gen_print(stderr, M, "%+10.6f", NULL,  NULL,NULL,NULL, NULL,NULL,NULL, "\n");
  }

void h2tt_do_check_pmap
  ( char *name,
    hr2_pmap_t *M, 
    double eps, 
    char *msg, 
    char *file, 
    int32_t lnum,
    const char *func
  )
  { r3x3_t P; r3x3_mul(&(M->dir), &(M->inv), &P);
    double mv = r3x3_L_inf_normalize(&P);
    if (mv == 0)
      { fprintf(stderr, " ** %s: {dir} x {inv} is the zero matrix\n", name);
        programerror(msg, file, lnum, func);
      }
    for (int32_t i = 0; i < NH; i++)
      for (int32_t j = 0; j < NH; j++)
        { double vexp = (i == j ? 1.0 : 0.0);
          double err = P.c[i][j] - vexp;
          if (fabs(err) > eps)
            { fprintf(stderr, " ** %s: {dir} x {inv} is not the identity\n", name);
              fprintf(stderr, " element [%d][%d] of product is %24.16e", i, j, P.c[i][j]);
              fprintf(stderr, " error = %+20.16e (max %20.16e)\n", err, eps);
              programerror(msg, file, lnum, func);
            }
        }
  } 

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
  )
  {
    r3_t qn; (void)r3_dir(q, &qn);
    r3_t pM;
    if (anti)
      { double d2min = +INF;
        for (int32_t s = -1; s <= +1; s += 2)
          { r3_t u = (r3_t){{s*p->c[0], s*p->c[1], s*p->c[2]}};
            r3_t uM;
            if (row) 
              { r3x3_map_row(&u, M, &uM); }
            else
              { r3x3_map_col(M, &u, &uM); }
            (void)r3_dir(&uM, &uM);
            double d2 = r3_dist_sqr(&uM, &qn);
            if (d2 < d2min) { d2min = d2; pM = uM; }
          }
      }
    else if (twirl)
      { double d2min = +INF;
        for (int32_t sy = -1; sy <= +1; sy += 2)
          for (int32_t sx = -1; sx <= +1; sx += 2)
            for (int32_t sw = -1; sw <= +1; sw += 2)
              { r3_t u = (r3_t){{sw*p->c[0], sx*p->c[1], sy*p->c[2]}};
                r3_t uM;
                if (row) 
                  { r3x3_map_row(&u, M, &uM); }
                else
                  { r3x3_map_col(M, &u, &uM); }
                (void)r3_dir(&uM, &uM);
                double d2 = r3_dist_sqr(&uM, &qn);
                if (d2 < d2min) { d2min = d2; pM = uM; }
              }
      }
    else
      { if (row) 
          { r3x3_map_row(p, M, &pM); }
        else
          { r3x3_map_col(M, p, &pM); }
        (void)r3_dir(&pM, &pM);
      }
    h2tt_do_check_r3_norm_eps(name, &pM, &qn, 0.00000001, msg, file, lnum, func);
  }  

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
  )
  { 
    r3x3_t *N = (inv ? &(M->inv) : &(M->dir));
    h2tt_do_check_map_r3(name, &(p->c), anti, twirl, TRUE, N, &(q->c), msg, file, lnum, func);
  }

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
  )
  {
    r3x3_t *N = (inv ? &(M->dir) : &(M->inv));
    h2tt_do_check_map_r3(name, &(A->f), anti, twirl, FALSE, N, &(B->f), msg, file, lnum, func);
  }
  
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
  )
  { 
    hr2_point_t hp = hr2_from_r2(p);
    hr2_point_t hq = hr2_from_r2(q);
    h2tt_do_check_pmap_point(name, &hp, TRUE, FALSE, M, inv, &hq, msg, file, lnum, func);
  }

void h2tt_do_check_num_eps
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

void h2tt_do_check_r2_eps
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

void h2tt_do_check_r3_norm_eps
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

