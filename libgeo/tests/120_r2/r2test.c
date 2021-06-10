/* r2test --- test program for r2.h, r2x2.h  */
/* Last edited on 2021-06-09 20:37:51 by jstolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>
#include <flt.h>
#include <r2.h>
#include <r2_extra.h>
#include <r2x2.h>
#include <rn_test_tools.h>


#define N 2
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_r2(int32_t verbose);
void test_r2maps(int32_t verbose);
void test_r2x2(int32_t verbose);
void throw_matrix(r2x2_t *m);
void throw_diag_matrix(r2x2_t *m);
void throw_symmetric_matrix(r2x2_t *m);

int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_r2(i < 3);
    for (i = 0; i < 100; i++) test_r2maps(i < 3);
    for (i = 0; i < 100; i++) test_r2x2(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_r2(int32_t verbose)
  {
    r2_t a, b, c, d, e, para, perp;
    double r, s, t;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    if (verbose)
      { fprintf(stderr,
          "sizeof(r2_t) = %lud  %d*sizeof(double) = %lud\n",
          sizeof(r2_t), N, N*sizeof(double)
        );
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_zero ---\n"); }
    r2_zero(&a);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i],0.0, NO, NO, "r2_zero error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_all ---\n"); }
    r2_all(3.14, &a);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i],3.14, NO, NO, "r2_all error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_axis ---\n"); }
    for (k = 0; k < N; k++)
      { r2_axis(k, &a);
        for (i = 0; i < N; i++)
          { rn_check_eq(a.c[i],(i == k ? 1.0 : 0.0), NO, NO, "r2_axis error"); }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_throw_cube ---\n"); }
    r2_throw_cube(&a);
    for (i = 0; i < N; i++)
      { affirm(a.c[i] != a.c[(i+1)%N], "r2_throw probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "r2_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "r2_throw error(2)"); 
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_throw_dir ---\n"); }
    /* Should check uniformity... */
    r2_throw_dir(&a);
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r2_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rn_check_eps(1,rr,0.000000001 * rr, NO, NO, "r2_throw_dir error (2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_throw_ball ---\n"); }
    /* Should check uniformity... */
    r2_throw_ball(&a);
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r2_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.000000001*rr, "r2_throw_ball error (2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_add ---\n"); }
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_add(&a, &b, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "r2_add error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_sub ---\n"); }
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_sub(&a, &b, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "r2_sub error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_neg ---\n"); }
    r2_throw_cube(&a);
    r2_neg(&a, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],- a.c[i], NO, NO, "r2_neg error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_scale ---\n"); }
    s = drandom();
    r2_throw_cube(&a);
    r2_scale(s, &a, &d);
    for (i = 0; i < N; i++)
      { double zi = s*a.c[i];
        rn_check_eq(d.c[i],zi, NO, NO, "r2_scale error(1)");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_mix ---\n"); }
    s = drandom();
    t = drandom();
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_mix(s, &a, t, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = s * a.c[i] + t * b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r2_mix error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_mix_in ---\n"); }
    s = drandom();
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    d = b;
    r2_mix_in(s, &a, &d);
    for (i = 0; i < N; i++)
      { double ddi = b.c[i] + s * a.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r2_mix_in error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_weigh ---\n"); }
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r2_weigh(&a, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = a.c[i] * b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r2_weigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_norm, r2_norm_sqr, r2_L_inf_norm ---\n"); }
    r2_throw_cube(&a);
    r = r2_norm(&a);
    s = r2_norm_sqr(&a);
    t = r2_L_inf_norm(&a);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
        if (ai > tt) { tt = ai; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "r2_norm error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "r2_norm_sqr error");
    rn_check_eq(t,tt, NO, NO, "r2_L_inf_norm error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_dist, r2_dist_sqr, r2_L_inf_dist ---\n"); }
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r = r2_dist(&a, &b);
    s = r2_dist_sqr(&a, &b);
    t = r2_L_inf_dist(&a, &b);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "r2_dist error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "r2_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "r2_L_inf_dist error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_dir, r2_L_inf_dir ---\n"); }
    r2_throw_cube(&a);
    r = r2_dir(&a, &b);
    s = r2_L_inf_dir(&a, &d);
    ss = r2_norm(&a);
    tt = r2_L_inf_norm(&a);
    for (i = 0; i < N; i++)
      { rn_check_eps(b.c[i],a.c[i]/ss,0.000000001 * ss, NO, NO, "r2_dir error");
        rn_check_eps(d.c[i],a.c[i]/tt,0.000000001 * tt, NO, NO, "r2_L_inf_dir error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_dot, r2_cos, r2_sin, r2_angle ---\n"); }
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r = r2_dot(&a, &b);
    double S = r2_sin(&a, &b);
    double C = r2_cos(&a, &b);
    double A = r2_angle(&a, &b);
    mag = sqrt(r2_dot(&a,&a)*r2_dot(&b,&b));
    rr = 0.0;
    for (i = 0; i < N; i++) { rr += a.c[i]*b.c[i]; }
    double CC = rr/(r2_norm(&a)*r2_norm(&b));
    rn_check_eps(r,rr,0.000000001 * mag, NO, NO, "r2_dot error(1)");
    rn_check_eps(C,CC,0.000000001, NO, NO, "r2_cos error(1)");
    d = a;
    r2_mix_in(-rr/r2_norm_sqr(&b), &b, &d);
    double SS = r2_norm(&d)/r2_norm(&a);
    rn_check_eps(S,SS,0.000000001, NO, NO, "r2_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.000000001, NO, NO, "r2_angle error(1)");
    for (i = 0; i < N; i++)
      { r2_axis(i, &a);
        for (j = 0; j < N; j++)
          { r2_axis(j, &b);
            r = r2_dot(&a, &b);
            s = r2_sin(&a, &b);
            t = r2_cos(&a, &b);
            rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, NO, NO, "r2_dot error(2)");
            rn_check_eq(t,rr, NO, NO, "r2_dot error(2)");
            rn_check_eq(s,1.0 - rr, NO, NO, "r2_dot error(2)");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        int32_t p;
        r2_axis(i0, &a);
        r2_cross(&a, &d);
        r2_axis(i1, &e);
        for (p = 0; p < N; p++)
          { double ep = sign*e.c[p];
            rn_check_eq(d.c[p],ep, NO, NO, "r2_cross error(x)");
          }
      }
    /* Test on random vectors: */
    r2_throw_cube(&a);
    r2_cross(&a, &d);
    mag = r2_norm(&a);
    r = r2_dot(&a, &d);
    rn_check_eps(r,0.0,0.00000001 * mag, NO, NO, "r2_cross error(1)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        double sign = ((i % 2) == 0 ? 1.0 : -1.0);
        r2_axis(i0, &a);
        r2_axis(i1, &b);
        r = r2_det(&a, &b);
        rn_check_eq(r,sign, NO, NO, "r2_det error(2)");
      }
    /* Test on random vectors: */
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r = r2_det(&a, &b);
    r2_cross(&a, &e);
    rr = r2_dot(&e, &b);
    mag = r2_norm(&a)*r2_norm(&b);
    rn_check_eps(r,rr,0.00000001 * mag, NO, NO, "r2_det error(1)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_decomp ---\n"); }
    r2_throw_cube(&a);
    r2_throw_cube(&b);
    r = r2_decomp(&a, &b, &para, &perp);
    rr = r2_dot(&a, &b)/r2_norm_sqr(&b);  
    rn_check_eps(r,rr,0.000000001 * (fabs(r) + fabs(rr)), NO, NO, "r2_decomp error(1)");
    r2_add(&para, &perp, &c);
    s = r2_dist(&a, &c);
    affirm (s <= 0.000000001 * r2_norm(&a), "r2_decomp error(2)");
    s = r2_dot(&perp, &b);
    rn_check_eps(s,0.0,0.000000001 * r2_norm(&b), NO, NO, "r2_decomp error(3)");
    t = r2_dot(&para, &perp);
    rn_check_eps(t,0.0,0.000000001 * r2_norm(&a), NO, NO, "r2_decomp error(4)");
    
    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_print ---\n"); }
    if (verbose)
      { r2_throw_cube (&a);
        fprintf(stderr, "a = ");
        r2_print(stderr, &a);
        fputc('\n', stderr);
      }

    if (verbose)
      { 
        fprintf(stderr, "!! r2_is_finite NOT TESTED\n");
        fprintf(stderr, "!! r2_eq NOT TESTED\n");
        fprintf(stderr, "!! r2_circumcenter NOT TESTED\n");
        fprintf(stderr, "!! r2_orient NOT TESTED\n");
        fprintf(stderr, "!! r2_incircle NOT TESTED\n");
        fprintf(stderr, "!! r2_throw_normal NOT TESTED\n");
      }
  }

void test_r2maps(int32_t verbose)
  {
    r2_t a, b, c;
    double r, s;
    int32_t i, j, k, ii, jj;
    bool_t debug = verbose;

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_map_projective ---\n"); }
    r3x3_t M = (r3x3_t){{{ 1.0, 2.0, 3.0 }, { 1.0, 1.0, 0.0}, { 0.0, 0.0, 1.0 }}};
    
    auto void do_project(r2_t *p, r2x2_t *J);
    void do_project(r2_t *p, r2x2_t *J)
      { r2_map_projective(p, &M, J); }
    
    for (ii = 0; ii <= 2; ii++)
      { for (jj = 0; jj <= 2; jj++)
          { double X = ii/2.0;
            double Y = jj/2.0;
            a = (r2_t){{ X, Y }};
            r2x2_t J;    
            b = a;
            r2x2_ident(&J);
            r2_map_projective(&b, &M, &J);
            r2_map_check_jacobian(&a, &do_project, "r2_map_projective", 1.0e-6, debug);
            c = (r2_t){{ (X + 2.0)/(X + 1.0), (Y + 3.0)/(X + 1.0) }};
            for (k = 0; k < N; k++)
              { rn_check_eps(b.c[k],c.c[k], 1.0e-8, NO, NO, "r2_map_projective error"); }
          }
      }
            
    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_map_twirl ---\n"); }
    int32_t NX = 640;
    int32_t NY = 480;
    r2_t ctr = (r2_t){{ 0.5*NX, 0.5*NY }};
    double rad = 0.25*fmin(NX,NY);
    double ang = 0.5*M_PI*(2*drandom()-1);
    
    auto void do_twirl(r2_t *p, r2x2_t *J);
    void do_twirl(r2_t *p, r2x2_t *J)
      { r2_map_twirl(p, &ctr, rad, ang, J); }
    
    for (ii = 0; ii <= 2; ii++)
      { for (jj = 0; jj <= 2; jj++)
          { double X = (0.10 + 0.80*ii/2.0)*NX;
            double Y = (0.10 + 0.80*jj/2.0)*NY;
            a = (r2_t){{ X, Y }};
            r2x2_t J;    
            b = a;
            r2x2_ident(&J);
            r2_map_twirl(&b, &ctr, rad, ang, &J);
            r2_map_check_jacobian(&a, &do_twirl, "r2_map_twirl", 1.0e-6, debug);
          }
      }
            
    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2_map_expand,r2_map_contract ---\n"); }
    double xlo = 2.0, xhi = 5.0;
    double ylo = 1.0, yhi = 3.0;
    
    auto void do_expand(r2_t *p, r2x2_t *J);
    void do_expand(r2_t *p, r2x2_t *J)
      { r2_map_expand(p, xlo, xhi, ylo, yhi, J); }
    
    auto void do_contract(r2_t *p, r2x2_t *J);
    void do_contract(r2_t *p, r2x2_t *J)
      { r2_map_contract(p, xlo, xhi, ylo, yhi, J); }
    
    for (ii = 0; ii <= 2; ii++)
      { for (jj = 0; jj <= 2; jj++)
          { double X = xlo + (0.05 + 0.90*jj/2)*(xhi - xlo);
            double Y = ylo + (0.05 + 0.90*ii/2)*(yhi - ylo);
            
            a = (r2_t){{ X, Y }};
            
            r2_map_check_jacobian(&a, &do_expand, "r2_map_expand", 1.0e-6, debug);
            
            r2x2_t JE;    
            b = a;
            r2x2_ident(&JE);
            r2_map_expand(&b, xlo, xhi, ylo, yhi, &JE);
            
            r2_map_check_jacobian(&b, &do_contract, "r2_map_contract", 1.0e-6, debug);
            
            c = b;
            r2x2_t JC = JE;    
            r2_map_contract(&c, xlo, xhi, ylo, yhi, &JC);

            r = r2_dist(&a, &c);  
            rn_check_eps(r,0.0,0.0000001, NO, NO, "r2_map_expand/r2_map_contract error(1)");

            r2x2_t K;
            r2x2_ident(&K);
            s = 0;
            for(i = 0; i < 2; i++)
              { for (j = 0; j < 2; j++)
                  { s += fabs(K.c[i][j] - JC.c[i][j]); }
              }
            rn_check_eps(s,0.0,0.000000001 * r2x2_norm(&JE), NO, NO, "r2_map_expand/r2_map_contract error(2)");
          }
      }

    if (verbose)
      { 
        fprintf(stderr, "!! r2_map_radial NOT TESTED\n");
      }
  }

void test_r2x2(int32_t verbose)
  {
    r2x2_t A, B, C;
    r2_t a, b, c, bb, cc;
    double r, s, t;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- Size and allocation ---\n"); }
    if (verbose)
      { fprintf(stderr,
          "sizeof(r2x2_t) = %lud  %d*%d*sizeof(double) = %lud\n",
          sizeof(r2x2_t), N, N, N*N*sizeof(double)
        );
        fprintf(stderr, "&B = %016lx\n", (long unsigned)&B);
        fprintf(stderr, "&A-&B = %lud\n", ((long unsigned)(&A))-((long unsigned)(&B)));
        fprintf(stderr, "&B-&C = %lud\n", ((long unsigned)(&B))-((long unsigned)(&C)));
        fprintf(stderr, "&(B.c) = %016lx\n", (long unsigned)&(B.c));
        fprintf(stderr, "B.c = %016lx\n", (long unsigned)(B.c));
        fprintf(stderr, "&(B.c[0]) = %016lx\n", (long unsigned)&(B.c[0]));
        fprintf(stderr, "B.c[0] = %016lx\n", (long unsigned)(B.c[0]));
        fprintf(stderr, "&(B.c[0][0]) = %016lx\n", (long unsigned)&(B.c[0][0]));
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- Indexing and addressing ---\n"); }
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
        { double *Aij = &(A.c[i][j]); 
          affirm(Aij = ((double *)&A)+(N*i)+j, "r2x2_t indexing error");
        }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_zero, r2x2_ident ---\n"); }
    r2x2_zero(&A);
    r2x2_ident(&B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(A.c[i][j],0.0, NO, NO, "r2x2_zero error"); 
            rn_check_eq(B.c[i][j],(i == j ? 1.0 : 0.0), NO, NO, "r2x2_ident error");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_get_row, r2x2_set_row, r2x2_get_col, r2x2_set_col ---\n"); }
    throw_matrix(&A);
    int32_t dir; /* 0 for row, 1 for col. */
    for (dir = 0; dir < 2; dir++)
      { for (i = 0; i < N; i++)
          { /* Check {r2x2_get_row,r2x2_get_col}: */
            r2_throw_cube(&a);
            if (dir == 0) { r2x2_get_row(&A, i, &a); } else { r2x2_get_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r2x2_get_row/r2x2_get_col error");
              }
            /* Check {r2x2_set_row,r2x2_set_col}: */
            r2_throw_cube(&a);
            if (dir == 0) { r2x2_set_row(&A, i, &a); } else { r2x2_set_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r2x2_set_row/r2x2_set_col error");
              }
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_map_row, r2x2_map_col ---\n"); }
    throw_matrix(&A);
    r2_throw_cube(&a);
    r2x2_map_row(&a, &A, &b);
    r2x2_map_col(&A, &a, &c);
    r2_zero(&bb);
    r2_zero(&cc);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { bb.c[j] += a.c[i] * A.c[i][j];
            cc.c[i] += A.c[i][j] * a.c[j];
          }
      }
    r = r2_dist(&b, &bb);
    affirm(r < 0.000000001 * r2_norm(&bb), "r2_map_row error");
    s = r2_dist(&c, &cc);
    affirm(s < 0.000000001 * r2_norm(&cc), "r2_map_col error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_scale ---\n"); }
    throw_matrix(&A);
    r = drandom();
    r2x2_scale(r, &A, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sel = r * A.c[i][j];
            rn_check_eps(C.c[i][j],sel,0.000000001 * fabs(sel), NO, NO,
              "r2x2_scale error"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_mul ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    r2x2_mul(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
            rn_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r2x2_mul error"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_mul_tr ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    r2x2_mul_tr(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[j][k]; }
            rn_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r2x2_mul error"
            );
          }
      }

     /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_transp ---\n"); }
    throw_matrix(&A);
    r2x2_transp(&A, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r2x2_transp error (1)"); }
      }
    /* In-place transpose: */
    B = A;
    r2x2_transp(&B, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r2x2_transp error (2)"); }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_det ---\n"); }
    throw_matrix(&A);
    for (i = 0; i < N; i++)
      { int32_t k = (i + 1) % N;
        for (j = 0; j < N; j++)
          { /* Check for linearity */
            r = drandom();
            A.c[i][j] = r;
            rr = r2x2_det(&A);

            s = drandom();
            A.c[i][j] = s;
            ss = r2x2_det(&A);

            t = drandom();
            A.c[i][j] = r*(1-t) + s*t;
            tt = r2x2_det(&A);
            mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_check_eps(tt,rr*(1.0 - t) + ss*t,000000001 * mag, NO, NO,
              "r2x2_det error(1)"
            );
          }

        /* Row swap test: */
        r = r2x2_det(&A);
        for (j = 0; j < N; j++)
          { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
        rr = r2x2_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "r2x2_det error(2)");

        /* Col swap test: */
        r = r2x2_det(&A);
        for (j = 0; j < N; j++)
          { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
        rr = r2x2_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "r2x2_det error(3)");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_inv ---\n"); }
    throw_matrix(&A);
    r2x2_inv(&A, &B);
    r2x2_mul(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double val = (i == j ? 1.0 : 0.0);
            affirm(fabs(C.c[i][j] - val) < 000000001, "r2x2_inv error");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_sqrt ---\n"); }
    do 
      { throw_matrix(&A); }
    while ((r2x2_det(&A) < 0) || (A.c[0][0] + A.c[1][1] + 2*sqrt(r2x2_det(&A)) <= 1.0e-10));
    r2x2_sqrt(&A, &B);
    r2x2_mul(&B, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { affirm(fabs(C.c[i][j] - A.c[i][j]) < 000000001, "r2x2_sqrt error"); }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_norm,r2x2_norm_sqr,r2x2_mod_norm,r2x2_mod_norm_sqr ---\n"); }
    throw_matrix(&A);
    s = r2x2_norm_sqr(&A);
    r = r2x2_norm(&A);
    t = r2x2_mod_norm_sqr(&A);
    ss = 0; tt = 0;
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            double Dij = (i == j ? Aij - 1 : Aij);
            ss += Aij*Aij;
            tt += Dij*Dij;
          }
      }
    affirm(ss >= 0, "r2x2_norm_sqr error");
    affirm(fabs(ss - s) < 000000001, "r2x2_norm_sqr error");
    rr = sqrt(ss);
    affirm(fabs(rr - r) < 000000001, "r2x2_norm error");
    affirm(tt >= 0, "r2x2_mod_norm_sqr error");
    affirm(fabs(tt - t) < 000000001, "r2x2_mod_norm_sqr error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_sym_eigen ---\n"); }
    if (drandom() < 0.2) 
      { /* Test with a diagonal matrix: */
        throw_diag_matrix(&A);
      }
    else
      { /* Test with a general symmetric matrix: */
        throw_symmetric_matrix(&A);
      }
    r2x2_sym_eigen(&A, &a, &B);
    // if (verbose)
    //   { fprintf(stderr, "testing r2x2_sym_eigen:\n");
    //     fprintf(stderr, "A = ");
    //     r2x2_print(stderr, &A);
    //     fputc('\n', stderr);
    //     fprintf(stderr, "eigenvalues = ");
    //     r2_print(stderr, &a);
    //     fputc('\n', stderr);
    //     fprintf(stderr, "eigenvectors = ");
    //     r2x2_print(stderr, &B);
    //     fputc('\n', stderr);
    //   }
    /* Check order of eigenvalues: */
    for (i = 1; i < N; i++)
      { affirm(a.c[i-1] >= a.c[i], "r2x2_sym_eigen error: order"); }
    /* Check whether {B} is orthonormal: */
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) 
              { sum += B.c[k][i]*B.c[k][j]; }
            double val = (i == j ? 1.0 : 0.0);
            rn_check_eps(val,sum,0.000000001 * fabs(sum), NO, NO, 
              "r2x2_sym_eigen error: not orthormal"
            );
          }
      }
    /* Check whether {B} is right-handed: */
    rr = r2x2_det(&B);
    affirm(fabs(rr - 1.0) < 000000001, "r2x2_sym_eigen error: not right-handed");
    /* Check whether {A = B'*DIAG(e)*B}: */
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) 
              { sum += B.c[k][i]*a.c[k]*B.c[k][j]; }
            rn_check_eps(A.c[i][j],sum,0.000000001 * fabs(sum), NO, NO, 
              "r2x2_sym_eigen error: decomp"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r2x2_print ---\n"); }
    if (verbose)
      { throw_matrix (&A);
        fprintf(stderr, "A = ");
        r2x2_print(stderr, &A);
        fputc('\n', stderr);
      }

    if (verbose)
      { 
        fprintf(stderr, "!! r2x2_add NOT TESTED\n");
        fprintf(stderr, "!! r2x2_sub NOT TESTED\n");
        fprintf(stderr, "!! r2x2_neg NOT TESTED\n");
        fprintf(stderr, "!! r2x2_mix NOT TESTED\n");
        fprintf(stderr, "!! r2x2_adj NOT TESTED\n");
        fprintf(stderr, "!! r2x2_is_unif_scaling NOT TESTED\n");
        fprintf(stderr, "!! r2x2_rot90 NOT TESTED\n");
        fprintf(stderr, "!! r2x2_rot_and_scale NOT TESTED\n");
        fprintf(stderr, "!! r2x2_moments NOT TESTED\n");
        fprintf(stderr, "!! r2x2_from_point_pairs NOT TESTED\n");
      }
  }  

void throw_matrix(r2x2_t *m)
  {
    int32_t i, j;
    r2_t a;
    for (i = 0; i < N; i++)
      { r2_throw_cube(&a);
        for (j = 0; j < N; j++) { m->c[i][j] = a.c[j]; }
      }
  }

void throw_diag_matrix(r2x2_t *m)
  {
    int32_t i, j;
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { m->c[i][j] = (i == j ? 2*drandom() - 1.0 : 0.0); }
      }
  }

void throw_symmetric_matrix(r2x2_t *m)
  {
    int32_t i, j;
    for (i = 0; i < N; i++)
      { /* Note: {j} runs to {i} not {N-1}! */
        for (j = 0; j <= i; j++)
          { m->c[i][j] = 2*drandom() - 1.0;
            if (j != i) { m->c[j][i] = m->c[i][j]; }
          }
      }
  }
