/* r3test --- test program for r3.h, r3x3.h  */
/* Last edited on 2023-10-01 19:25:48 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>
#include <flt.h>

#include <r3.h>
#include <r3_extra.h>
#include <r3x3.h>
#include <rn_test_tools.h>

#define N 3
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_r3(int32_t verbose);
void test_r3x3(int32_t verbose);
void throw_matrix(r3x3_t *m);
void print_matrix(char *name, r3x3_t *A);

void check_num_eps(char *name, double x, double y, double eps, char *msg);
  /* If {x} and {y} differ by more than {eps}, prints {name}, {x}, {y}, and {msg}, and stops. */

void check_regular_polyhedron(char *func, double R, double L, int32_t n, r3_t p[], int32_t deg);
  /* Check that {p[0..n-1]} are the vertices of a regular 
    polyhedron with radius {R}, side {L}, and vertex degree {deg}. */

void test_r3x3_diff_sqr(bool_t verbose);

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_r3(i < 3);
    for (i = 0; i < 100; i++) test_r3x3(i < 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_r3(int32_t verbose)
  { r3_t a, b, c, d, e, para, perp;
    double r, s, t;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    if (verbose)
      { fprintf(stderr,
          "sizeof(r3_t) = %lu  %d*sizeof(double) = %lu\n",
          sizeof(r3_t), N, N*sizeof(double)
        );
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_zero ---\n"); }
    r3_zero(&a);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i],0.0, NO, NO, "r3_zero error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_all ---\n"); }
    r3_all(3.14, &a);
    for (i = 0; i < N; i++)
      { rn_check_eq(a.c[i],3.14, NO, NO, "r3_all error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_axis ---\n"); }
    for (k = 0; k < N; k++)
      { r3_axis(k, &a);
        for (i = 0; i < N; i++)
          { rn_check_eq(a.c[i],(i == k ? 1.0 : 0.0), NO, NO, "r3_axis error"); }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_throw_cube ---\n"); }
    r3_throw_cube(&a);
    for (i = 0; i < N; i++)
      { affirm(a.c[i] != a.c[(i+1)%N], "r3_throw probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        double vv = a.c[i]*256.0;
        affirm(vv != floor(vv), "r3_throw error(3)"); 
        affirm((a.c[i] > -1.0) && (a.c[i] < 1.0), "r3_throw error(2)"); 
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_throw_dir ---\n"); }
    /* Should check uniformity... */
    r3_throw_dir(&a);
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r3_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    rn_check_eps(1,rr,0.000000001 * rr, NO, NO, "r3_throw_dir error (2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_throw_ball ---\n"); }
    /* Should check uniformity... */
    r3_throw_ball(&a);
    /* Check variation: */
    for (i = 0; i < N; i++) { affirm(a.c[i] != a.c[(i+1)%N], "r3_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (i = 0; i < N; i++) { double ai = a.c[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.000000001*rr, "r3_throw_ball error (2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_add ---\n"); }
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_add(&a, &b, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] + b.c[i], NO, NO, "r3_add error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_sub ---\n"); }
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_sub(&a, &b, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],a.c[i] - b.c[i], NO, NO, "r3_sub error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_neg ---\n"); }
    r3_throw_cube(&a);
    r3_neg(&a, &d);
    for (i = 0; i < N; i++)
      { rn_check_eq(d.c[i],- a.c[i], NO, NO, "r3_neg error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_scale ---\n"); }
    s = drandom();
    r3_throw_cube(&a);
    r3_scale(s, &a, &d);
    for (i = 0; i < N; i++)
      { double zi = s*a.c[i];
        rn_check_eq(d.c[i],zi, NO, NO, "r3_scale error(1)");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_mix ---\n"); }
    s = drandom();
    t = drandom();
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_mix(s, &a, t, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = s * a.c[i] + t * b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r3_mix error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_mix_in ---\n"); }
    s = drandom();
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    d = b;
    r3_mix_in(s, &a, &d);
    for (i = 0; i < N; i++)
      { double ddi = b.c[i] + s * a.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r3_mix_in error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_weigh ---\n"); }
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_weigh(&a, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = a.c[i] * b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r3_weigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_unweigh ---\n"); }
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_unweigh(&a, &b, &d);
    for (i = 0; i < N; i++)
      { double ddi = a.c[i] / b.c[i];
        rn_check_eq(d.c[i],ddi, NO, NO, "r3_unweigh error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_rot_axis ---\n"); }
    { r3_throw_cube(&a);
      int32_t i = int32_abrandom(0, N-1);
      int32_t j = int32_abrandom(0, N-2); if (j >= i) { j++; }
      double ang = 2.1*M_PI*drandom();
      r3_rot_axis(&a, i, j, ang, &d);
      rn_test_rot_axis(N, a.c, i, j, ang, d.c, "r3_rot_axis error");
    }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_norm, r3_norm_sqr, r3_L_inf_norm ---\n"); }
    r3_throw_cube(&a);
    r = r3_norm(&a);
    s = r3_norm_sqr(&a);
    t = r3_L_inf_norm(&a);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double ai = fabs(a.c[i]);
        ss += ai*ai; 
        if (ai > tt) { tt = ai; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "r3_norm error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "r3_norm_sqr error");
    rn_check_eq(t,tt, NO, NO, "r3_L_inf_norm error");

    if (verbose) { fprintf(stderr, "--- r3_dist, r3_dist_sqr, r3_L_inf_dist ---\n"); }
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r = r3_dist(&a, &b);
    s = r3_dist_sqr(&a, &b);
    t = r3_L_inf_dist(&a, &b);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < N; i++)
      { double di = fabs(a.c[i] - b.c[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000000001 * rr, NO, NO, "r3_dist error");
    rn_check_eps(s,ss,0.000000001 * ss, NO, NO, "r3_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "r3_L_inf_dist error");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_dir, r3_L_inf_dir ---\n"); }
    r3_throw_cube(&a);
    r = r3_dir(&a, &b);
    s = r3_L_inf_dir(&a, &d);
    ss = r3_norm(&a);
    tt = r3_L_inf_norm(&a);
    for (i = 0; i < N; i++)
      { rn_check_eps(b.c[i],a.c[i]/ss,0.000000001 * ss, NO, NO, "r3_dir error");
        rn_check_eps(d.c[i],a.c[i]/tt,0.000000001 * tt, NO, NO, "r3_L_inf_dir error");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_dot, r3_cos, r3_sin, r3_angle ---\n"); }
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r = r3_dot(&a, &b);
    double S = r3_sin(&a, &b);
    double C = r3_cos(&a, &b);
    double A = r3_angle(&a, &b);
    mag = sqrt(r3_dot(&a,&a)*r3_dot(&b,&b));
    rr = 0.0;
    for (i = 0; i < N; i++) { rr += a.c[i]*b.c[i]; }
    double CC = rr/(r3_norm(&a)*r3_norm(&b));
    rn_check_eps(r,rr,0.000000001 * mag, NO, NO, "r3_dot error(1)");
    rn_check_eps(C,CC,0.000000001, NO, NO, "r3_cos error(1)");
    d = a;
    r3_mix_in(-rr/r3_norm_sqr(&b), &b, &d);
    double SS = r3_norm(&d)/r3_norm(&a);
    rn_check_eps(S,SS,0.000000001, NO, NO, "r3_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.000000001, NO, NO, "r3_angle error(1)");
    for (i = 0; i < N; i++)
      { r3_axis(i, &a);
        for (j = 0; j < N; j++)
          { r3_axis(j, &b);
            r = r3_dot(&a, &b);
            s = r3_sin(&a, &b);
            t = r3_cos(&a, &b);
            rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, NO, NO, "r3_dot error(2)");
            rn_check_eq(t,rr, NO, NO, "r3_dot error(2)");
            rn_check_eq(s,1.0 - rr, NO, NO, "r3_dot error(2)");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        int32_t p;
        r3_axis(i0, &a);
        r3_axis(i1, &b);
        r3_cross(&a, &b, &d);
        r3_axis(i2, &e);
        for (p = 0; p < N; p++)
          { double ep = e.c[p];
            rn_check_eq(d.c[p],ep, NO, NO, "r3_cross error(x)");
          }
      }
    /* Test on random vectors: */
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_cross(&a, &b, &d);
    mag = r3_norm(&a)*r3_norm(&b);
    r = r3_dot(&a, &d);
    rn_check_eps(r,0.0,0.00000001 * mag*r3_norm(&a), NO, NO, "r3_cross error(1)");
    r = r3_dot(&b, &d);
    rn_check_eps(r,0.0,0.00000001 * mag*r3_norm(&b), NO, NO, "r3_cross error(2)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_pick_ortho ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        r3_axis(i0, &a);
        double ma = r3_L_inf_norm(&a);
        double mo = r3_pick_ortho(&a, &d);
        rn_check_eq(ma, mo, NO, NO, "r3_pick_ortho error(0)");
        r = r3_dot(&a, &d);
        rn_check_eps(r, 0.0, 0.00000001*ma, NO, NO, "r3_pick_ortho error(1)");
      }
    /* Test on random vectors: */
    { r3_throw_cube(&a);
      double ma = r3_L_inf_norm(&a);
      double mo = r3_pick_ortho(&a,&d);
      rn_check_eq(ma, mo, NO, NO, "r3_pick_ortho error(2)");
      r = r3_dot(&a, &d);
      rn_check_eps(r, 0.0, 0.00000001*ma, NO, NO, "r3_pick_ortho error(3)");
    }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        r3_axis(i0, &a);
        r3_axis(i1, &b);
        r3_axis(i2, &c);
        r = r3_det(&a, &b, &c);
        rn_check_eq(r,1.0, NO, NO, "r3_det error(2)");
      }
    /* Test on random vectors: */
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r3_throw_cube(&c);
    r = r3_det(&a, &b, &c);
    r3_cross(&a, &b, &e);
    rr = r3_dot(&e, &c);
    mag = r3_norm(&a)*r3_norm(&b)*r3_norm(&c);
    rn_check_eps(r,rr,0.00000001 * mag, NO, NO, "r3_det error(1)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_decomp ---\n"); } 
    r3_throw_cube(&a);
    r3_throw_cube(&b);
    r = r3_decomp(&a, &b, &para, &perp);
    rr = r3_dot(&a, &b)/r3_norm_sqr(&b);  
    rn_check_eps(r,rr,0.000000001 * (fabs(r) + fabs(rr)), NO, NO, "r3_decomp error(1)");
    r3_add(&para, &perp, &c);
    s = r3_dist(&a, &c);
    affirm (s <= 0.000000001 * r3_norm(&a), "r3_decomp error(2)");
    s = r3_dot(&perp, &b);
    rn_check_eps(s,0.0,0.000000001 * r3_norm(&b), NO, NO, "r3_decomp error(3)");
    t = r3_dot(&para, &perp);
    rn_check_eps(t,0.0,0.000000001 * r3_norm(&a), NO, NO, "r3_decomp error(4)");

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3_print ---\n"); }
    if (verbose)
      { r3_throw_cube (&a);
        fprintf(stderr, "a = ");
        r3_print(stderr, &a);
        fputc('\n', stderr);
      }
      
    /* Checking the regular polyhedra: */
    { r3_t pv[20];
      double R = drandom();

      /* ---------------------------------------------------------------------- */
      if (verbose) { fprintf(stderr, "--- r3_tetrahedron_vertices ---\n"); }
      r3_tetrahedron_vertices(R, 4, pv);
      double tetra_L = R*sqrt(8.0/3.0);
      check_regular_polyhedron("r3_tetrahedron_vertices", R, tetra_L, 4, pv, 3);

      /* ---------------------------------------------------------------------- */
      if (verbose) { fprintf(stderr, "--- r3_octahedron_vertices ---\n"); }
      r3_octahedron_vertices(R, 6, pv);
      double octa_L = R*sqrt(2.0);
      check_regular_polyhedron("r3_octahedron_vertices", R, octa_L, 6, pv, 4);

      /* ---------------------------------------------------------------------- */
      if (verbose) { fprintf(stderr, "--- r3_hexahedron_vertices ---\n"); }
      r3_hexahedron_vertices(R, 8, pv);
      double hexa_L = 2*R/sqrt(3.0);
      check_regular_polyhedron("r3_hexahedron_vertices", R, hexa_L, 8, pv, 3);

      /* ---------------------------------------------------------------------- */
      if (verbose) { fprintf(stderr, "--- r3_icosahedron_vertices ---\n"); }
      r3_icosahedron_vertices(R, 12, pv);
      double icosa_L = R*sqrt(2 - 2/sqrt(5));
      check_regular_polyhedron("r3_icosahedron_vertices", R, icosa_L, 12, pv, 5);

      /* ---------------------------------------------------------------------- */
      if (verbose) { fprintf(stderr, "--- r3_dodecahedron_vertices ---\n"); }
      r3_dodecahedron_vertices(R, 20, pv);
      double dodeca_L = R*(sqrt(5)-1)/sqrt(3);
      check_regular_polyhedron("r3_dodecahedron_vertices", R, dodeca_L, 20, pv, 3);
   }

      /* ---------------------------------------------------------------------- */
    /* Checking the cylindrical mesh: */
    if (verbose) { fprintf(stderr, "!! warning: r3_cylindrical_grid not tested\n"); }
  }

void test_r3x3(int32_t verbose)
  {
    r3x3_t A, B, C;
    r3_t a, b, c, bb, cc;
    double r, s, t;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- Size and allocation ---\n"); }
    if (verbose)
      { fprintf(stderr,
          "sizeof(r3x3_t) = %lu  %d*%d*sizeof(double) = %lu\n",
          sizeof(r3x3_t), N, N, N*N*sizeof(double)
        );
        fprintf(stderr, "&B = %016lx\n", (long unsigned)&B);
        fprintf(stderr, "&A-&B = %lu\n", ((long unsigned)(&A))-((long unsigned)(&B)));
        fprintf(stderr, "&B-&C = %lu\n", ((long unsigned)(&B))-((long unsigned)(&C)));
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
          affirm(Aij == ((double *)&A)+(N*i)+j, "r3x3_t indexing error");
        }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_zero, r3x3_ident ---\n"); }
    r3x3_zero(&A);
    r3x3_ident(&B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(A.c[i][j],0.0, NO, NO, "r3x3_zero error"); 
            rn_check_eq(B.c[i][j],(i == j ? 1.0 : 0.0), NO, NO, "r3x3_ident error");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_get_row, r3x3_set_row, r3x3_get_col, r3x3_set_col ---\n"); }
    throw_matrix(&A);
    int32_t dir; /* 0 for row, 1 for col. */
    for (dir = 0; dir < 2; dir++)
      { for (i = 0; i < N; i++)
          { /* Check {r3x3_get_row,r3x3_get_col}: */
            r3_throw_cube(&a);
            if (dir == 0) { r3x3_get_row(&A, i, &a); } else { r3x3_get_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r3x3_get_row/r3x3_get_col error");
              }
            /* Check {r3x3_set_row,r3x3_set_col}: */
            r3_throw_cube(&a);
            if (dir == 0) { r3x3_set_row(&A, i, &a); } else { r3x3_set_col(&A, i, &a); }
            for (j = 0; j < N; j++)
              { double vj = (dir == 0 ? A.c[i][j] : A.c[j][i]);
                affirm(vj = a.c[j], "r3x3_set_row/r3x3_set_col error");
              }
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_map_row, r3x3_map_col ---\n"); }
    throw_matrix(&A);
    r3_throw_cube(&a);
    r3x3_map_row(&a, &A, &b);
    r3x3_map_col(&A, &a, &c);
    r3_zero(&bb);
    r3_zero(&cc);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { bb.c[j] += a.c[i] * A.c[i][j];
            cc.c[i] += A.c[i][j] * a.c[j];
          }
      }
    r = r3_dist(&b, &bb);
    affirm(r < 0.000000001 * r3_norm(&bb), "r3_map_row error");
    s = r3_dist(&c, &cc);
    affirm(s < 0.000000001 * r3_norm(&cc), "r3_map_col error");

    if (verbose) { fprintf(stderr, "--- r3x3_scale ---\n"); }
    throw_matrix(&A);
    r = drandom();
    r3x3_scale(r, &A, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sel = r * A.c[i][j];
            rn_check_eps(C.c[i][j],sel,0.000000001 * fabs(sel), NO, NO,
              "r3x3_scale error"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_mul ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    r3x3_mul(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
            rn_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r3x3_mul error"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_mul_tr ---\n"); }
    throw_matrix(&A);
    throw_matrix(&B);
    r3x3_mul_tr(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double sum = 0.0;
            for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[j][k]; }
            rn_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
              "r3x3_mul error"
            );
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_transp ---\n"); }
    throw_matrix(&A);
    r3x3_transp(&A, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r3x3_transp error (1)"); }
      }
    /* In-place transpose: */
    B = A;
    r3x3_transp(&B, &B);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { rn_check_eq(B.c[i][j],A.c[j][i], NO, NO, "r3x3_transp error (2)"); }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_det ---\n"); }
    throw_matrix(&A);
    for (i = 0; i < N; i++)
      { int32_t k = (i + 1) % N;
        for (j = 0; j < N; j++)
          { /* Check for linearity */
            r = drandom();
            A.c[i][j] = r;
            rr = r3x3_det(&A);

            s = drandom();
            A.c[i][j] = s;
            ss = r3x3_det(&A);

            t = drandom();
            A.c[i][j] = r*(1-t) + s*t;
            tt = r3x3_det(&A);
            mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_check_eps(tt,(rr*(1 - t) + ss*t),000000001 * mag, NO, NO,
              "r3x3_det error(1)"
            );
          }

        /* Row swap test: */
        r = r3x3_det(&A);
        for (j = 0; j < N; j++)
          { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
        rr = r3x3_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "r3x3_det error(2)");

        /* Col swap test: */
        r = r3x3_det(&A);
        for (j = 0; j < N; j++)
          { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
        rr = r3x3_det(&A);
        mag = fabs(r) + fabs(rr);
        rn_check_eps(r,-rr,000000001 * mag, NO, NO, "r3x3_det error(3)");
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_inv ---\n"); }
    throw_matrix(&A);
    r3x3_inv(&A, &B);
    r3x3_mul(&A, &B, &C);
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double val = (i == j ? 1.0 : 0.0);
            affirm(fabs(C.c[i][j] - val) < 000000001, "r3x3_inv error");
          }
      }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_norm,r3x3_norm_sqr,r3x3_mod_norm ---\n"); }
    throw_matrix(&A);
    s = r3x3_norm_sqr(&A);
    r = r3x3_norm(&A);
    t = r3x3_mod_norm_sqr(&A);
    ss = 0; tt = 0;
    for (i = 0; i < N; i++)
      { for (j = 0; j < N; j++)
          { double Aij = A.c[i][j];
            double Dij = (i == j ? Aij - 1 : Aij);
            ss += Aij*Aij;
            tt += Dij*Dij;
          }
      }
    affirm(ss >= 0, "r3x3_norm_sqr error");
    affirm(fabs(ss - s) < 000000001, "r3x3_norm_sqr error");
    rr = sqrt(ss);
    affirm(fabs(rr - r) < 000000001, "r3x3_norm error");
    affirm(tt >= 0, "r3x3_mod_norm_sqr error");
    affirm(fabs(tt - t) < 000000001, "r3x3_mod_norm_sqr error");

    /* ---------------------------------------------------------------------- */
    /* ---------------------------------------------------------------------- */
    test_r3x3_diff_sqr(verbose);

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_u_v_rotation ---\n"); }
    r3_throw_dir(&a);
    r3_throw_dir(&b);
    r3x3_u_v_rotation(&a, &b, &A);
    r3x3_map_row(&a, &A, &c);
    for (i = 0; i < N; i++)
      { affirm(fabs(b.c[i] - c.c[i]) < 000000001, "r3x3_u_v_rotation error"); }

    /* ---------------------------------------------------------------------- */
    if (verbose) { fprintf(stderr, "--- r3x3_print ---\n"); }
    if (verbose)
      { throw_matrix (&A);
        fprintf(stderr, "A = ");
        r3x3_print(stderr, &A);
        fputc('\n', stderr);
      }

    if (verbose)
      { 
        fprintf(stderr, "!! r3x3_add NOT TESTED\n");
        fprintf(stderr, "!! r3x3_sub NOT TESTED\n");
        fprintf(stderr, "!! r3x3_neg NOT TESTED\n");
        fprintf(stderr, "!! r3x3_mix NOT TESTED\n");
        fprintf(stderr, "!! r3x3_adj NOT TESTED\n");
        fprintf(stderr, "!! r3x3_is_unif_scaling NOT TESTED\n");
        fprintf(stderr, "!! r3_bezier.h NOT TESTED\n");
        fprintf(stderr, "!! r3_path.h NOT TESTED\n");
        fprintf(stderr, "!! r3_motion.h NOT TESTED\n");
      }
  }  

void test_r3x3_diff_sqr(bool_t verbose)
  {
    bool_t debug = FALSE;
    if (verbose) { fprintf(stderr, "--- r3x3_diff_sqr ---\n"); }
    
    r3x3_t A, B, R;
    throw_matrix(&A);
    throw_matrix(&B);
    throw_matrix(&R);
    int32_t iz = int32_abrandom(0, 2);
    int32_t jz = int32_abrandom(0, 2);
    R.c[iz][jz] = 0;
    if (debug)
      { print_matrix("A", &A);
        print_matrix("B", &B);
        print_matrix("R", &R);
      }
    double dabs2, drel2;
    r3x3_diff_sqr(&A, &B, &R, &dabs2, &drel2);
    double cabs2 = 0, crel2 = 0;
    for (int32_t j = 0; j < N;  j++)
      { for (int32_t i = 0; i < N;  i++)
          { double rij = R.c[i][j];
            if (rij != 0.0)
              { double d = A.c[i][j] - B.c[i][j];
                cabs2 += d*d;
                crel2 += (d/rij)*(d/rij);
              }
          }
      }
    check_num_eps("drel2", drel2, crel2, 0.0000001, "r3x3_diff_sqr failed");
    check_num_eps("dabs2", dabs2, cabs2, 0.0000001, "r3x3_diff_sqr failed");
  }

void throw_matrix(r3x3_t *m)
  { int32_t i, j;
    r3_t a;
    for (i = 0; i < N; i++)
      { r3_throw_cube(&a);
        for (j = 0; j < N; j++) { m->c[i][j] = a.c[j]; }
      }
  }
  
void print_matrix(char *name, r3x3_t *A)
  { fprintf(stderr, "%s = ", name);
    r3x3_gen_print(stderr, A, "%+10.6f", NULL,NULL,NULL, NULL,NULL,NULL);
    fprintf(stderr, "\n");
  }

void check_regular_polyhedron(char *func, double R, double L, int32_t n, r3_t p[], int32_t deg)
  {
    /* Not a complete test... */
    R = fabs(R); /* We can't distinguish {R} from {-R}. */
    int32_t i, j;
    /* Find the smallest nonzero vertex-vertex distance {dmin}: */
    double dmin = +INF;
    for (i = 0; i < n; i++)
      { for (j = 0; j < i; j++)
          { dmin = fmin(dmin, r3_dist(&(p[i]), &(p[j]))); }
      }
    rn_check_eps(dmin, L, 0.000001*R, NULL, NULL, "polyhedron has wrong side");
    /* Check each vertex: */
    for (i = 0; i < n; i++)
      { /* Check distance from origin: */
        double Ri = r3_norm(&(p[i]));
        rn_check_eps(Ri, R, 0.000001*R, &i, NULL, "vertex has wrong radius");
        /* Check number of nearest neighbors: */
        int32_t degi = 0;
        for (j = 0; j < n; j++)
          { double dij = r3_dist(&(p[i]), &(p[j]));
            if (fabs(dij - L) < 0.0001*R) { degi++; }
          }
        rn_check_eq(degi, deg, &i, NULL, "vertex has wrong degree");
      }
  }

void check_num_eps(char *name, double x, double y, double eps, char *msg)
  { double diff = fabs(x - y);
    if (diff > eps)
      { fprintf(stderr, " ** %s: %+20.16e %+20.16e", name, x, y);
        fprintf(stderr, " diff = %+20.16e  max = %+20.16e - %s\n", diff, eps, msg);
        exit(1);
      }
  }
