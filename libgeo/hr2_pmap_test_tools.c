/* See {hr2_pmap_test_tools.h}. */
/* Last edited on 2024-11-21 16:47:23 by stolfi */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <jsrandom.h>
#include <jsprintf.h>
#include <affirm.h>

#include <rn.h>
#include <r3x3.h>
#include <r3.h>
#include <r2x2.h>
#include <r2.h>
#include <hr2.h>
#include <hr2_pmap.h>
#include <hr2_test_tools.h>
#include <rmxn_throw.h>

#include <rmxn_test_tools.h>

#include <hr2_pmap_test_tools.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */
  
#define Pr fprintf
#define Er stderr


void hr2_pmap_test_tools_throw_aff_map(hr2_pmap_t *M)
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
  
hr2_pmap_t hr2_pmap_test_tools_throw_almost_singular(double detMax)
  { 
    bool_t debug = FALSE;
    if (debug) { Pr(Er, "    > --- %s ---\n", __FUNCTION__); }
    hr2_pmap_t M;
    rmxn_throw_almost_singular_pair(NH, &(M.dir.c[0][0]), &(M.inv.c[0][0]), detMax);
    if (debug) { Pr(Er, "    < --- %s ---\n", __FUNCTION__); }
    return M;
  }

hr2_pmap_t hr2_pmap_test_tools_throw_non_singular(double detMin)
  { 
    bool_t debug = FALSE;
    if (debug) { Pr(Er, "    > --- %s ---\n", __FUNCTION__); }
    hr2_pmap_t M;
    rmxn_throw_non_singular_pair(NH, &(M.dir.c[0][0]), &(M.inv.c[0][0]), detMin);
    if (debug) { Pr(Er, "    < --- %s ---\n", __FUNCTION__); }
    return M;
  }
  
void hr2_pmap_test_tools_print(uint32_t indent, char *name, hr2_pmap_t *M)
  { char *olp = jsprintf("%*s", indent, "");
    if (name != NULL) { Pr(Er, "%*s%s = \n", indent, "", name); }
    hr2_pmap_gen_print(Er, M, "%+24.16e", NULL,  olp, "  ", "\n", "[ ", " ", " ]", "\n");
    free(olp);
  }
  
void hr2_pmap_test_tools_print_matrix(uint32_t indent, char *name, r3x3_t *A)
  { char *ilp = jsprintf("%*s [ ", indent, "");
    if (name != NULL) { Pr(Er, "%*s%s = \n", indent, "", name); }
    r3x3_gen_print(Er, A, "%+24.16e", "",  "", "", ilp, " ", " ]\n");
    free(ilp);
  }

void hr2_pmap_test_tools_do_check
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
      { Pr(Er, " ** %s: {dir} x {inv} is the zero matrix\n", name);
        programerror(msg, file, lnum, func);
      }
    for (int32_t i = 0; i < NH; i++)
      for (int32_t j = 0; j < NH; j++)
        { double vexp = (i == j ? 1.0 : 0.0);
          double err = P.c[i][j] - vexp;
          if (fabs(err) > eps)
            { Pr(Er, " ** %s: {dir} x {inv} is not the identity\n", name);
              Pr(Er, " element [%d][%d] of product is %24.16e", i, j, P.c[i][j]);
              Pr(Er, " error = %+20.16e (max %20.16e)\n", err, eps);
              programerror(msg, file, lnum, func);
            }
        }
  } 

void hr2_pmap_test_tools_do_check_map_r3
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
    char *xsgn = "";
    char *xroc = (row ? "row" : "col");
    if (anti)
      { xsgn = "anti ";
        double d2min = +INF;
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
      { xsgn = "twirl ";
        double d2min = +INF;
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
    char *msg1 = jsprintf("%s (%s%s)", msg, xsgn, xroc);
    hr2_test_tools_do_check_r3_norm_eps(name, &pM, &qn, 0.00000001, msg1, file, lnum, func);
    free(msg1);
  }  

void hr2_pmap_test_tools_do_check_map_point
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
    char *xiod = (inv ? "inv" : "dir");
    char *msg1 = jsprintf("%s (%s)", msg, xiod);
    hr2_pmap_test_tools_do_check_map_r3(name, &(p->c), anti, twirl, TRUE, N, &(q->c), msg1, file, lnum, func);
    free(msg1);
  }

void hr2_pmap_test_tools_do_check_map_line
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
    hr2_pmap_test_tools_do_check_map_r3(name, &(A->f), anti, twirl, FALSE, N, &(B->f), msg, file, lnum, func);
  }
  
void hr2_pmap_test_tools_do_check_map_r2_point
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
    hr2_pmap_test_tools_do_check_map_point(name, &hp, TRUE, FALSE, M, inv, &hq, msg, file, lnum, func);
  }

void hr2_pmap_test_tools_check_perturbed
  ( hr2_pmap_t *M,
    hr2_pmap_test_tools_check_proc_t *OK,
    r3x3_t *P,
    bool_t OK_exp,
    char *OK_name,
    bool_t verbose
  )
  { if (verbose) { Pr(Er, "        checking unperturbed map"); }
    hr2_pmap_test_tools_check_matrix(M, "unperturbed map", OK, OK_name, OK_exp);
    if (verbose) { hr2_pmap_test_tools_print_matrix(8, "perturbation matrix {P}:", P); }
    double Pmax = r3x3_L_inf_norm(P);
    if (verbose) { Pr(Er, "        max perturbation = %24.16e\n", Pmax); }
    uint32_t NT = 5; /* Perturbation trials. */
    for (int32_t trial = 0; trial < NT; trial++)
      { for (int32_t dir = 0; dir <= 1; dir++)
          { hr2_pmap_t Q = (*M);
            r3x3_t *A = (dir == 0 ? &(Q.dir) : &(Q.inv));
            char *dtag = (dir == 0 ? "dir" : "inv");
            double Aenorm = r3x3_norm(A)/3; /* RMS elem value. */
            char *Q_name = jsprintf("map {Q} that is {M} with {M.%s} perturbed by {%23.16e * random * P}", dtag, Aenorm);
            if (verbose) { Pr(Er, "        checking %s\n", Q_name); }
            for (int32_t i = 0; i < NH; i++)
              { for (int32_t j = 0; j < NH; j++)
                 { A->c[i][j] += Aenorm*dabrandom(-1,+1)*P->c[i][j]; }
              }
            r3x3_t *B = (dir == 0 ? &(Q.inv) : &(Q.dir));
            r3x3_inv(A, B);
            if (verbose) { hr2_pmap_test_tools_print(8, NULL, &Q); }
            hr2_pmap_test_tools_check_matrix(&Q, Q_name, OK, OK_name, OK_exp);
            free(Q_name);
          }
      }
  }
          
void  hr2_pmap_test_tools_check_matrix
  ( hr2_pmap_t *Q,
    char *Q_name,
    hr2_pmap_test_tools_check_proc_t *OK,
    char *OK_name,
    bool_t OK_exp
  )     
  {
    bool_t OK_cmp = OK(Q);
    if (OK_cmp != OK_exp)
      { char *xok_cmp = (OK_cmp ? "true" : "false");
        char *xok_exp = (OK_exp ? "true" : "false");
        char *msg1 = jsprintf("{%s} is %s (should be %s) for %s", OK_name, xok_cmp, xok_exp, Q_name);
        hr2_pmap_test_tools_print(8, msg1, Q);
        char *msg2 = jsprintf("{%s} failed (0)", OK_name);
        demand(FALSE, msg2);
      }
  }
