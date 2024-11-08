/* See {hr2test_tools.h}. */
/* Last edited on 2024-11-04 06:36:44 by stolfi */

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
#include <hr2_test_tools.h>

#include <hr2_pmap_test_tools.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

void hr2_test_throw_pmap(hr2_pmap_t *M)
  {
    r3_t a;
    for (int32_t i = 0; i < NH; i++)
      { r3_throw_cube(&a);
        for (int32_t j = 0; j < NH; j++) { M->dir.c[i][j] = a.c[j]; }
      }
    r3x3_inv(&(M->dir), &(M->inv));
  }

void hr2_test_throw_aff_map(hr2_pmap_t *M)
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
  
void hr2_test_print_pmap(char *name, hr2_pmap_t *M)
  { fprintf(stderr, "%s = \n", name);
    hr2_pmap_gen_print(stderr, M, "%+20.15f", NULL,  NULL,NULL,NULL, NULL,NULL,NULL, "\n");
  }

void hr2_test_do_check_pmap
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
    hr2_test_do_check_r3_norm_eps(name, &pM, &qn, 0.00000001, msg, file, lnum, func);
  }  

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
  )
  { 
    r3x3_t *N = (inv ? &(M->inv) : &(M->dir));
    hr2_test_do_check_map_r3(name, &(p->c), anti, twirl, TRUE, N, &(q->c), msg, file, lnum, func);
  }

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
  )
  {
    r3x3_t *N = (inv ? &(M->dir) : &(M->inv));
    hr2_test_do_check_map_r3(name, &(A->f), anti, twirl, FALSE, N, &(B->f), msg, file, lnum, func);
  }
  
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
  )
  { 
    hr2_point_t hp = hr2_from_r2(p);
    hr2_point_t hq = hr2_from_r2(q);
    hr2_test_do_check_pmap_point(name, &hp, TRUE, FALSE, M, inv, &hq, msg, file, lnum, func);
  }

void hr2_test_perturbed_pmap(hr2_pmap_t *M, hr2_test_pmap_check_proc_t *ok, r3x3_t *P, char *fname)
  { for (int32_t dir = 0; dir <= 1; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        double mag = 0.1 + 0.2*drandom();
        r3x3_scale(mag, A, A);
        for (int32_t i = 0; i < 3; i++)
          for (int32_t j = 0; j < 3; j++)
            { double Asave = A->c[i][j];
              double eps = P->c[i][j];
              if (eps != 0)
                { A->c[i][j] += eps*mag;
                  if (ok(M))
                    { char *dtag = (dir == 0 ? "dir" : "inv");
                      char *elm = NULL; asprintf(&elm, "{%s} is true even perturbing {M.%s[%d][%d]} by %24.16e", fname, dtag, i, j, eps*mag);
                      hr2_test_print_pmap(elm, M);
                      char *msg = NULL; asprintf(&msg, "{%s} failed (0)", fname);
                      demand(FALSE, msg);
                    }
                  A->c[i][j] = Asave;
                }
            }
      }
  }
