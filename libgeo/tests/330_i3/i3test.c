/* i3test --- test program for i3.h, i3x3.h  */
/* Last edited on 2021-06-09 19:52:47 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>
#include <in_test_tools.h>

#include <i3.h>
/* #include <i3x3.h> */

#define N 3
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_i3(int32_t verbose);
/* void test_i3x3(int32_t verbose); */
/* void throw_matrix(i3x3_t *m); */

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    srand(1993);
    srandom(1993);
    for (i = 0; i < 100; i++) test_i3(i < 3);
    /* for (i = 0; i < 100; i++) test_i3x3(i < 3); */
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_i3(int32_t verbose)
  { i3_t a, b, c, d, e;
    int32_t t32;
    int32_t tt32;
    int64_t r64, s64;
    int64_t rr64, ss64;
    int32_t trad = 4615;
    int32_t i, j, k;

    if (verbose)
      { fprintf(stderr,
          "sizeof(i3_t) = %lud  %d*sizeof(int32_t) = %lud\n",
          sizeof(i3_t), N, N*sizeof(int32_t)
        );
      }

    if (verbose) { fprintf(stderr, "--- i3_zero ---\n"); }
    i3_zero(&a);
    for (i = 0; i < N; i++)
      { in_check_eq(a.c[i], 0, NO, NO, "i3_zero error"); }

    if (verbose) { fprintf(stderr, "--- i3_all ---\n"); }
    i3_all(4634, &a);
    for (i = 0; i < N; i++)
      { in_check_eq(a.c[i], 4634, NO, NO, "i3_all error"); }

    if (verbose) { fprintf(stderr, "--- i3_axis ---\n"); }
    for (k = 0; k < N; k++)
      { i3_axis(k, &a);
        for (i = 0; i < N; i++)
          { in_check_eq(a.c[i], (i == k), NO, NO, "i3_axis error"); }
      }

    if (verbose) { fprintf(stderr, "--- i3_throw_cube ---\n"); }
    i3_throw_cube(trad, &a);
    for (i = 0; i < N; i++)
      { /* There is a small chance that this test will fail by chance: */
        affirm(a.c[i] != a.c[(i+1)%N], "i3_throw probable error(1)"); 
        affirm((a.c[i] >= -trad) && (a.c[i] <= +trad), "i3_throw error(2)"); 
      }

    if (verbose) { fprintf(stderr, "--- i3_add ---\n"); }
    i3_throw_cube(trad, &a);
    i3_throw_cube(trad, &b);
    i3_add(&a, &b, &d);
    for (i = 0; i < N; i++)
      { in_check_eq(d.c[i], a.c[i] + b.c[i], NO, NO, "i3_add error"); }

    if (verbose) { fprintf(stderr, "--- i3_sub ---\n"); }
    i3_throw_cube(trad, &a);
    i3_throw_cube(trad, &b);
    i3_sub(&a, &b, &d);
    for (i = 0; i < N; i++)
      { in_check_eq(d.c[i], a.c[i] - b.c[i], NO, NO, "i3_sub error"); }

    if (verbose) { fprintf(stderr, "--- i3_neg ---\n"); }
    i3_throw_cube(trad, &a);
    i3_neg(&a, &d);
    for (i = 0; i < N; i++)
      { in_check_eq(d.c[i], -a.c[i], NO, NO, "i3_neg error"); }

    if (verbose) { fprintf(stderr, "--- i3_norm_sqr, i3_L_inf_norm ---\n"); }
    i3_throw_cube(trad, &a);
    s64 = i3_norm_sqr(&a);
    t32 = i3_L_inf_norm(&a);
    ss64 = 0;
    tt32 = 0;
    for (i = 0; i < N; i++)
      { int32_t ai = abs(a.c[i]);
        ss64 += ((int64_t)ai)*((int64_t)ai); 
        if (ai > tt32) { tt32 = ai; }
      }
    in_check_eq(s64, ss64, NO, NO, "i3_norm_sqr error");
    in_check_eq(t32, tt32, NO, NO, "i3_L_inf_norm error");

    if (verbose) { fprintf(stderr, "--- i3_dist_sqr, i3_L_inf_dist ---\n"); }
    i3_throw_cube(trad, &a);
    i3_throw_cube(trad, &b);
    s64 = i3_dist_sqr(&a, &b);
    t32 = i3_L_inf_dist(&a, &b);
    ss64 = 0;
    tt32 = 0;
    for (i = 0; i < N; i++)
      { int32_t di = abs(a.c[i] - b.c[i]);
        ss64 += di*(int64_t)di; 
        if (di > tt32) { tt32 = di; }
      }
    in_check_eq(s64, ss64, NO, NO, "i3_dist_sqr error");
    in_check_eq(t32, tt32, NO, NO, "i3_L_inf_dist error");

    if (verbose) { fprintf(stderr, "--- i3_dot ---\n"); }
    i3_throw_cube(trad, &a);
    i3_throw_cube(trad, &b);
    r64 = i3_dot(&a, &b);
    rr64 = 0;
    for (i = 0; i < N; i++) { rr64 += a.c[i]*(int64_t)b.c[i]; }
    in_check_eq(r64, rr64, NO, NO, "i3_dot error(1)");
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { i3_axis(i, &a);
        for (j = 0; j < N; j++)
          { i3_axis(j, &b);
            r64 = i3_dot(&a, &b);
            rr64 = (i == j);
            in_check_eq(r64, rr64, NO, NO, "i3_dot error(2)");
          }
      }

    if (verbose) { fprintf(stderr, "--- i3_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        int32_t p;
        i3_axis(i0, &a);
        i3_axis(i1, &b);
        i3_cross(&a, &b, &d);
        i3_axis(i2, &e);
        for (p = 0; p < N; p++)
          { int64_t ep = e.c[p];
            in_check_eq(d.c[p], ep, NO, NO, "i3_cross error(x)");
          }
      }
    /* Test on random vectors: */
    i3_throw_cube(trad, &a);
    i3_throw_cube(trad, &b);
    i3_cross(&a, &b, &d);
    r64 = i3_dot(&a, &d);
    in_check_eq(r64, 0, NO, NO, "i3_cross error(1)");

    if (verbose) { fprintf(stderr, "--- i3_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t i2 = (i + 2) % N;
        i3_axis(i0, &a);
        i3_axis(i1, &b);
        i3_axis(i2, &c);
        r64 = i3_det(&a, &b, &c);
        in_check_eq(r64, 1, NO, NO, "i3_det error(2)");
      }
    /* Test on random vectors: */
    i3_throw_cube(trad, &a);
    i3_throw_cube(trad, &b);
    i3_throw_cube(trad, &c);
    r64 = i3_det(&a, &b, &c);
    i3_cross(&a, &b, &e);
    rr64 = i3_dot(&e, &c);
    in_check_eq(r64, rr64, NO, NO, "i3_det error(1)");

    if (verbose) { fprintf(stderr, "--- i3_print ---\n"); }
    if (verbose)
      { i3_throw_cube(trad, &a);
        fprintf(stderr, "a = ");
        i3_print(stderr, &a);
        fputc('\n', stderr);
      }

    if (verbose)
      { 
        fprintf(stderr, "!! i3_eq NOT TESTED\n");
      }
  }

// void test_i3x3(int32_t verbose)
//   {
//     i3x3_t A, B, C;
//     i3_t a, b, c, bb, cc;
//     double r, s, t;
//     double rr, ss, tt;
//     double mag;
//     int32_t i, j, k;
// 
//     if (verbose) { fprintf(stderr, "--- Size and allocation ---\n"); }
//     if (verbose)
//       { fprintf(stderr,
//           "sizeof(i3x3_t) = %d  %d*%d*sizeof(double) = %d\n",
//           sizeof(i3x3_t), N, N, N*N*sizeof(double)
//         );
//         fprintf(stderr, "&B = %08x\n", (unsigned)&B);
//         fprintf(stderr, "&A-&B = %d\n", ((unsigned)(&A))-((unsigned)(&B)));
//         fprintf(stderr, "&B-&C = %d\n", ((unsigned)(&B))-((unsigned)(&C)));
//         fprintf(stderr, "&(B.c) = %08x\n", (unsigned)&(B.c));
//         fprintf(stderr, "B.c = %08x\n", (unsigned)(B.c));
//         fprintf(stderr, "&(B.c[0]) = %08x\n", (unsigned)&(B.c[0]));
//         fprintf(stderr, "B.c[0] = %08x\n", (unsigned)(B.c[0]));
//         fprintf(stderr, "&(B.c[0][0]) = %08x\n", (unsigned)&(B.c[0][0]));
//       }
// 
//     if (verbose) { fprintf(stderr, "--- Indexing and addressing ---\n"); }
//     for (i = 0; i < N; i++)
//       for (j = 0; j < N; j++)
//         { double *Aij = &(A.c[i][j]); 
//           affirm(Aij = ((double *)&A)+(N*i)+j, "i3x3_t indexing error");
//         }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_zero, i3x3_ident ---\n"); }
//     i3x3_zero(&A);
//     i3x3_ident(&B);
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { in_check_eq(A.c[i][j],0, NO, NO, "i3x3_zero error"); 
//             in_check_eq(B.c[i][j],(i == j ? 1.0 : 0), NO, NO, "i3x3_ident error");
//           }
//       }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_map_row, i3x3_map_col ---\n"); }
//     throw_matrix(&A);
//     i3_throw_cube(trad, &a);
//     i3x3_map_row(&a, &A, &b);
//     i3x3_map_col(&A, &a, &c);
//     i3_zero(&bb);
//     i3_zero(&cc);
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { bb.c[j] += a.c[i] * A.c[i][j];
//             cc.c[i] += A.c[i][j] * a.c[j];
//           }
//       }
//     r = i3_dist(&b, &bb);
//     affirm(r < 0.000000001 * i3_norm(&bb), "i3_map_row error");
//     s = i3_dist(&c, &cc);
//     affirm(s < 0.000000001 * i3_norm(&cc), "i3_map_col error");
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_scale ---\n"); }
//     throw_matrix(&A);
//     r = drandom();
//     i3x3_scale(r, &A, &C);
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { double sel = r * A.c[i][j];
//             	in_check_eps(C.c[i][j],sel,0.000000001 * fabs(sel), NO, NO,
//               "i3x3_scale error"
//             );
//           }
//       }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_mul ---\n"); }
//     throw_matrix(&A);
//     throw_matrix(&B);
//     i3x3_mul(&A, &B, &C);
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { double sum = 0.0;
//             for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
//             	in_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
//               "i3x3_mul error"
//             );
//           }
//       }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_mul_tr ---\n"); }
//     throw_matrix(&A);
//     throw_matrix(&B);
//     i3x3_mul_tr(&A, &B, &C);
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { double sum = 0.0;
//             for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[j][k]; }
//             	in_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
//               "i3x3_mul error"
//             );
//           }
//       }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_transp ---\n"); }
//     throw_matrix(&A);
//     i3x3_transp(&A, &B);
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { in_check_eq(B.c[i][j],A.c[j][i], NO, NO, "i3x3_transp error (1)"); }
//       }
//     /* In-place transpose: */
//     B = A;
//     i3x3_transp(&B, &B);
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { in_check_eq(B.c[i][j],A.c[j][i], NO, NO, "i3x3_transp error (2)"); }
//       }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_det ---\n"); }
//     throw_matrix(&A);
//     for (i = 0; i < N; i++)
//       { int32_t k = (i + 1) % N;
//         for (j = 0; j < N; j++)
//           { /* Check for linearity */
//             r = drandom();
//             A.c[i][j] = r;
//             rr = i3x3_det(&A);
// 
//             s = drandom();
//             A.c[i][j] = s;
//             ss = i3x3_det(&A);
// 
//             t = drandom();
//             A.c[i][j] = r*(1-t) + s*t;
//             tt = i3x3_det(&A);
//             mag = fabs(rr) + fabs(ss) + fabs(tt);
//             	in_check_eps(tt,(rr*(1 - t) + ss*t),000000001 * mag, NO, NO,
//               "i3x3_det error(1)"
//             );
//           }
// 
//         /* Row swap test: */
//         r = i3x3_det(&A);
//         for (j = 0; j < N; j++)
//           { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
//         rr = i3x3_det(&A);
//         mag = fabs(r) + fabs(rr);
//         	in_check_eps(r,-rr,000000001 * mag, NO, NO, "i3x3_det error(2)");
// 
//         /* Col swap test: */
//         r = i3x3_det(&A);
//         for (j = 0; j < N; j++)
//           { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
//         rr = i3x3_det(&A);
//         mag = fabs(r) + fabs(rr);
//         	in_check_eps(r,-rr,000000001 * mag, NO, NO, "i3x3_det error(3)");
//       }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_inv ---\n"); }
//     throw_matrix(&A);
//     i3x3_inv(&A, &B);
//     i3x3_mul(&A, &B, &C);
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { double val = (i == j ? 1.0 : 0.0);
//             affirm(fabs(C.c[i][j] - val) < 000000001, "i3x3_inv error");
//           }
//       }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_norm,i3x3_norm_sqr,i3x3_mod_norm ---\n"); }
//     throw_matrix(&A);
//     s = i3x3_norm_sqr(&A);
//     r = i3x3_norm(&A);
//     t = i3x3_mod_norm_sqr(&A);
//     ss = 0; tt = 0;
//     for (i = 0; i < N; i++)
//       { for (j = 0; j < N; j++)
//           { double Aij = A.c[i][j];
//             double Dij = (i == j ? Aij - 1 : Aij);
//             ss += Aij*Aij;
//             tt += Dij*Dij;
//           }
//       }
//     affirm(ss >= 0, "i3x3_norm_sqr error");
//     affirm(fabs(ss - s) < 000000001, "i3x3_norm_sqr error");
//     rr = sqrt(ss);
//     affirm(fabs(rr - r) < 000000001, "i3x3_norm error");
//     affirm(tt >= 0, "i3x3_mod_norm_sqr error");
//     affirm(fabs(tt - t) < 000000001, "i3x3_mod_norm_sqr error");
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_u_v_rotation ---\n"); }
//     i3_throw_dir(&a);
//     i3_throw_dir(&b);
//     i3x3_u_v_rotation(&a, &b, &A);
//     i3x3_map_row(&a, &A, &c);
//     for (i = 0; i < N; i++)
//       { affirm(fabs(b.c[i] - c.c[i]) < 000000001, "i3x3_u_v_rotation error"); }
// 
//     if (verbose) { fprintf(stderr, "--- i3x3_print ---\n"); }
//     if (verbose)
//       { throw_matrix (&A);
//         fprintf(stderr, "A = ");
//         i3x3_print(stderr, &A);
//         fputc('\n', stderr);
//       }
//   }  
// 
// void throw_matrix(i3x3_t *m)
//   { int32_t i, j;
//     i3_t a;
//     for (i = 0; i < N; i++)
//       { i3_throw_cube(trad, &a);
//         for (j = 0; j < N; j++) { m->c[i][j] = a.c[j]; }
//       }
//   }

