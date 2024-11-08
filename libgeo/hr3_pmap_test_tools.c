/* See {hr3_pmap_test_tools.h}. */
/* Last edited on 2024-11-04 06:37:01 by stolfi */

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
#include <hr3_pmap_test_tools.h>

void hr3_test_throw_pmap(hr3_pmap_t *M)
  {
    r4_t a;
    for (int32_t i = 0; i < NH; i++)
      { r4_throw_cube(&a);
        for (int32_t j = 0; j < NH; j++) { M->dir.c[i][j] = a.c[j]; }
      }
    r4x4_inv(&(M->dir), &(M->inv));
  }

void hr3_test_throw_aff_map(hr3_pmap_t *M)
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
  
void hr3_test_print_pmap(char *name, hr3_pmap_t *M)
  { fprintf(stderr, "%s = ", name);
    hr3_pmap_gen_print(stderr, M, "%+10.6f", NULL,  NULL,NULL,NULL, NULL,NULL,NULL, "\n");
  }

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
  )
  {
    r4_t qn; (void)r4_dir(q, &qn);
    r4_t pM;
    if (twirl)
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
    hr3_test_do_check_r4_norm_eps(name, &pM, &qn, 0.00000001, msg, file, lnum, func);
  }  

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
  )
  { 
    r4x4_t *N = (inv ? &(M->inv) : &(M->dir));
    hr3_test_do_check_r4_map(name, &(p->c), twirl, TRUE, N, &(q->c), msg, file, lnum, func);
  }

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
  )
  {
    r4x4_t *N = (inv ? &(M->dir) : &(M->inv));
    hr3_test_do_check_r4_map(name, &(A->f), twirl, FALSE, N, &(B->f), msg, file, lnum, func);
  }
  
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
  )
  { 
    r3_t pM;
    if (inv)
      { pM = hr3_pmap_inv_r3_point(p, M); }
    else
      { pM = hr3_pmap_r3_point(p, M); }
    hr3_test_do_check_r3_eps(name, &pM, q, 0.00000001, msg, file, lnum, func);
  }


