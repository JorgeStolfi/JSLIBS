/* test_i2 --- test program for i2.h, i2x2.h  */
/* Last edited on 2024-11-08 12:43:35 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include <affirm.h>
#include <jsrandom.h>
#include <in_test_tools.h>

#include <i2.h>
/* #include <i2x2.h> */

#define N 2
#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_i2(int32_t verbose);
/* void test_i2x2(int32_t verbose); */
/* void throw_matrix(i2x2_t *m); */
/* void throw_diag_matrix(i2x2_t *m); */
/* void throw_symmetric_matrix(i2x2_t *m); */

int32_t main (int32_t argc, char **argv)
  {
    int32_t i;
    srand(1993);
    srandom(1993);

    for (i = 0; i < 100; i++) test_i2(i < 3);
    /* for (i = 0; i < 100; i++) test_i2x2(i < 3); */
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_i2(int32_t verbose)
  {
    i2_t a, b, d, e;
    int32_t t32;
    int32_t tt32;
    int64_t r64, s64;
    int64_t rr64, ss64;
    int32_t trad = 4634;
    int32_t i, j, k;

    if (verbose)
      { fprintf(stderr,
          "sizeof(i2_t) = %lu  %d*sizeof(int32_t) = %lu\n",
          sizeof(i2_t), N, N*sizeof(int32_t)
        );
      }

    if (verbose) { fprintf(stderr, "--- i2_zero ---\n"); }
    i2_zero(&a);
    for (i = 0; i < N; i++) { in_check_eq(a.c[i], 0, NO, NO, "i2_zero error"); }

    if (verbose) { fprintf(stderr, "--- i2_all ---\n"); }
    i2_all(4615, &a);
    for (i = 0; i < N; i++) { in_check_eq(a.c[i], 4615, NO, NO, "i2_all error"); }

    if (verbose) { fprintf(stderr, "--- i2_axis ---\n"); }
    for (k = 0; k < N; k++)
      { i2_axis(k, &a);
        for (i = 0; i < N; i++) { in_check_eq(a.c[i], (i == k), NO, NO, "i2_axis error"); }
      }

    if (verbose) { fprintf(stderr, "--- i2_throw_cube ---\n"); }
    i2_throw_cube(trad, &a);
    for (i = 0; i < N; i++)
      { /* There is a small chance that this test will fail by chance: */
        affirm(a.c[i] != a.c[(i+1)%N], "i2_throw probable error(1)"); 
        affirm((a.c[i] >= -trad) && (a.c[i] <= +trad), "i2_throw error(2)"); 
      }
    
    if (verbose) { fprintf(stderr, "--- i2_add ---\n"); }
    i2_throw_cube(trad, &a);
    i2_throw_cube(trad, &b);
    i2_add(&a, &b, &d);
    for (i = 0; i < N; i++)
      { in_check_eq(d.c[i] ,a.c[i] + b.c[i], NO, NO, "i2_add error"); }

    if (verbose) { fprintf(stderr, "--- i2_sub ---\n"); }
    i2_throw_cube(trad, &a);
    i2_throw_cube(trad, &b);
    i2_sub(&a, &b, &d);
    for (i = 0; i < N; i++)
      { in_check_eq(d.c[i], a.c[i] - b.c[i], NO, NO, "i2_sub error"); }

    if (verbose) { fprintf(stderr, "--- i2_neg ---\n"); }
    i2_throw_cube(trad, &a);
    i2_neg(&a, &d);
    for (i = 0; i < N; i++)
      { in_check_eq(d.c[i], -a.c[i], NO, NO, "i2_neg error"); }

    if (verbose) { fprintf(stderr, "--- i2_norm_sqr, i2_L_inf_norm ---\n"); }
    i2_throw_cube(trad, &a);
    s64 = i2_norm_sqr(&a);
    t32 = i2_L_inf_norm(&a);
    ss64 = 0;
    tt32 = 0;
    for (i = 0; i < N; i++)
      { int32_t ai = abs(a.c[i]);
        ss64 += ai*(int64_t)ai; 
        if (ai > tt32) { tt32 = ai; }
      }
    in_check_eq(s64, ss64, NO, NO, "i2_norm_sqr error");
    in_check_eq(t32, tt32, NO, NO, "i2_L_inf_norm error");

    if (verbose) { fprintf(stderr, "--- i2_dist_sqr, i2_L_inf_dist ---\n"); }
    i2_throw_cube(trad, &a);
    i2_throw_cube(trad, &b);
    s64 = i2_dist_sqr(&a, &b);
    t32 = i2_L_inf_dist(&a, &b);
    ss64 = 0;
    tt32 = 0;
    for (i = 0; i < N; i++)
      { int32_t di = abs(a.c[i] - b.c[i]);
        ss64 += di*(int64_t)di; 
        if (di > tt32) { tt32 = di; }
      }
    in_check_eq(s64, ss64, NO, NO, "i2_dist_sqr error");
    in_check_eq(t32, tt32, NO, NO, "i2_L_inf_dist error");

    if (verbose) { fprintf(stderr, "--- i2_dot ---\n"); }
    i2_throw_cube(trad, &a);
    i2_throw_cube(trad, &b);
    r64 = i2_dot(&a, &b);
    rr64 = 0;
    for (i = 0; i < N; i++) { rr64 += a.c[i]*(int64_t)b.c[i]; }
    in_check_eq(r64, rr64, NO, NO, "i2_dot error(1)");
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { i2_axis(i, &a);
        for (j = 0; j < N; j++)
          { i2_axis(j, &b);
            r64 = i2_dot(&a, &b);
            rr64 = (i == j);
            in_check_eq(r64,rr64, NO, NO, "i2_dot error(2)");
          }
      }

    if (verbose) { fprintf(stderr, "--- i2_cross ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t sign = ((i % 2) == 0 ? +1 : -1);
        int32_t p;
        i2_axis(i0, &a);
        i2_cross(&a, &d);
        i2_axis(i1, &e);
        for (p = 0; p < N; p++)
          { int64_t ep = sign*e.c[p];
            in_check_eq(d.c[p], ep, NO, NO, "i2_cross error(x)");
          }
      }
    /* Test on random vectors: */
    i2_throw_cube(trad, &a);
    i2_cross(&a, &d);
    r64 = i2_dot(&a, &d);
    in_check_eq(r64, 0, NO, NO, "i2_cross error(1)");

    if (verbose) { fprintf(stderr, "--- i2_det ---\n"); }
    /* Test on basis vectors: */
    for (i = 0; i < N; i++)
      { int32_t i0 = (i + 0) % N;
        int32_t i1 = (i + 1) % N;
        int32_t sign = ((i % 2) == 0 ? +1 : -1);
        i2_axis(i0, &a);
        i2_axis(i1, &b);
        r64 = i2_det(&a, &b);
        in_check_eq(r64, sign, NO, NO, "i2_det error(2)");
      }
    /* Test on random vectors: */
    i2_throw_cube(trad, &a);
    i2_throw_cube(trad, &b);
    r64 = i2_det(&a, &b);
    i2_cross(&a, &e);
    rr64 = i2_dot(&e, &b);
    in_check_eq(r64,rr64, NO, NO, "i2_det error(1)");

    if (verbose) { fprintf(stderr, "--- i2_print ---\n"); }
    if (verbose)
      { i2_throw_cube(trad, &a);
        fprintf(stderr, "a = ");
        i2_print(stderr, &a);
        fputc('\n', stderr);
      }

    if (verbose)
      { 
        fprintf(stderr, "!! i2_eq NOT TESTED\n");
      }
  }

//  void test_i2x2(int32_t verbose)
//    {
//      i2x2_t A, B, C;
//      i2_t a, b, c, bb, cc;
//      double r, s, t;
//      double rr, ss, tt;
//      double mag;
//      int32_t i, j, k;
//  
//      if (verbose) { fprintf(stderr, "--- Size and allocation ---\n"); }
//      if (verbose)
//        { fprintf(stderr,
//            "sizeof(i2x2_t) = %d  %d*%d*sizeof(double) = %d\n",
//            sizeof(i2x2_t), N, N, N*N*sizeof(double)
//          );
//          fprintf(stderr, "&B = %08x\n", (unsigned)&B);
//          fprintf(stderr, "&A-&B = %d\n", ((unsigned)(&A))-((unsigned)(&B)));
//          fprintf(stderr, "&B-&C = %d\n", ((unsigned)(&B))-((unsigned)(&C)));
//          fprintf(stderr, "&(B.c) = %08x\n", (unsigned)&(B.c));
//          fprintf(stderr, "B.c = %08x\n", (unsigned)(B.c));
//          fprintf(stderr, "&(B.c[0]) = %08x\n", (unsigned)&(B.c[0]));
//          fprintf(stderr, "B.c[0] = %08x\n", (unsigned)(B.c[0]));
//          fprintf(stderr, "&(B.c[0][0]) = %08x\n", (unsigned)&(B.c[0][0]));
//        }
//  
//      if (verbose) { fprintf(stderr, "--- Indexing and addressing ---\n"); }
//      for (i = 0; i < N; i++)
//        for (j = 0; j < N; j++)
//          { double *Aij = &(A.c[i][j]); 
//            affirm(Aij == ((double *)&A)+(N*i)+j, "i2x2_t indexing error");
//          }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_zero, i2x2_ident ---\n"); }
//      i2x2_zero(&A);
//      i2x2_ident(&B);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { in_check_eq(A.c[i][j],0.0, NO, NO, "i2x2_zero error"); 
//              in_check_eq(B.c[i][j],(i == j ? 1.0 : 0.0), NO, NO, "i2x2_ident error");
//            }
//        }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_map_row, i2x2_map_col ---\n"); }
//      throw_matrix(&A);
//      i2_throw_cube(trad, &a);
//      i2x2_map_row(&a, &A, &b);
//      i2x2_map_col(&A, &a, &c);
//      i2_zero(&bb);
//      i2_zero(&cc);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { bb.c[j] += a.c[i] * A.c[i][j];
//              cc.c[i] += A.c[i][j] * a.c[j];
//            }
//        }
//      r = i2_dist(&b, &bb);
//      affirm(r < 0.000000001 * i2_norm(&bb), "i2_map_row error");
//      s = i2_dist(&c, &cc);
//      affirm(s < 0.000000001 * i2_norm(&cc), "i2_map_col error");
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_scale ---\n"); }
//      throw_matrix(&A);
//      r = drandom();
//      i2x2_scale(r, &A, &C);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { double sel = r * A.c[i][j];
//              in_check_eps(C.c[i][j],sel,0.000000001 * fabs(sel), NO, NO,
//                "i2x2_scale error"
//              );
//            }
//        }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_mul ---\n"); }
//      throw_matrix(&A);
//      throw_matrix(&B);
//      i2x2_mul(&A, &B, &C);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { double sum = 0.0;
//              for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[k][j]; }
//              in_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
//                "i2x2_mul error"
//              );
//            }
//        }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_mul_tr ---\n"); }
//      throw_matrix(&A);
//      throw_matrix(&B);
//      i2x2_mul_tr(&A, &B, &C);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { double sum = 0.0;
//              for (k = 0; k < N; k++) { sum += A.c[i][k]*B.c[j][k]; }
//              in_check_eps(C.c[i][j],sum,0.000000001 * fabs(sum), NO, NO,
//                "i2x2_mul error"
//              );
//            }
//        }
//  
//       if (verbose) { fprintf(stderr, "--- i2x2_transp ---\n"); }
//      throw_matrix(&A);
//      i2x2_transp(&A, &B);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { in_check_eq(B.c[i][j],A.c[j][i], NO, NO, "i2x2_transp error (1)"); }
//        }
//      /* In-place transpose: */
//      B = A;
//      i2x2_transp(&B, &B);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { in_check_eq(B.c[i][j],A.c[j][i], NO, NO, "i2x2_transp error (2)"); }
//        }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_det ---\n"); }
//      throw_matrix(&A);
//      for (i = 0; i < N; i++)
//        { int32_t k = (i + 1) % N;
//          for (j = 0; j < N; j++)
//            { /* Check for linearity */
//              r = drandom();
//              A.c[i][j] = r;
//              rr = i2x2_det(&A);
//  
//              s = drandom();
//              A.c[i][j] = s;
//              ss = i2x2_det(&A);
//  
//              t = drandom();
//              A.c[i][j] = r*(1-t) + s*t;
//              tt = i2x2_det(&A);
//              mag = fabs(rr) + fabs(ss) + fabs(tt);
//              in_check_eps(tt,rr*(1.0 - t) + ss*t,000000001 * mag, NO, NO,
//                "i2x2_det error(1)"
//              );
//            }
//  
//          /* Row swap test: */
//          r = i2x2_det(&A);
//          for (j = 0; j < N; j++)
//            { double t = A.c[i][j]; A.c[i][j] = A.c[k][j]; A.c[k][j] = t; }
//          rr = i2x2_det(&A);
//          mag = fabs(r) + fabs(rr);
//          in_check_eps(r,-rr,000000001 * mag, NO, NO, "i2x2_det error(2)");
//  
//          /* Col swap test: */
//          r = i2x2_det(&A);
//          for (j = 0; j < N; j++)
//            { double t = A.c[j][i]; A.c[j][i] = A.c[j][k]; A.c[j][k] = t; }
//          rr = i2x2_det(&A);
//          mag = fabs(r) + fabs(rr);
//          in_check_eps(r,-rr,000000001 * mag, NO, NO, "i2x2_det error(3)");
//        }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_inv ---\n"); }
//      throw_matrix(&A);
//      i2x2_inv(&A, &B);
//      i2x2_mul(&A, &B, &C);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { double val = (i == j ? 1.0 : 0.0);
//              affirm(fabs(C.c[i][j] - val) < 000000001, "i2x2_inv error");
//            }
//        }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_sqrt ---\n"); }
//      do 
//        { throw_matrix(&A); }
//      while ((i2x2_det(&A) < 0) || (A.c[0][0] + A.c[1][1] + 2*sqrt(i2x2_det(&A)) <= 1.0e-10));
//      i2x2_sqrt(&A, &B);
//      i2x2_mul(&B, &B, &C);
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { affirm(fabs(C.c[i][j] - A.c[i][j]) < 000000001, "i2x2_sqrt error"); }
//        }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_norm,i2x2_norm_sqr,i2x2_mod_norm ---\n"); }
//      throw_matrix(&A);
//      s = i2x2_norm_sqr(&A);
//      r = i2x2_norm(&A);
//      t = i2x2_mod_norm_sqr(&A);
//      ss = 0; tt = 0;
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { double Aij = A.c[i][j];
//              double Dij = (i == j ? Aij - 1 : Aij);
//              ss += Aij*Aij;
//              tt += Dij*Dij;
//            }
//        }
//      affirm(ss >= 0, "i2x2_norm_sqr error");
//      affirm(fabs(ss - s) < 000000001, "i2x2_norm_sqr error");
//      rr = sqrt(ss);
//      affirm(fabs(rr - r) < 000000001, "i2x2_norm error");
//      affirm(tt >= 0, "i2x2_mod_norm_sqr error");
//      affirm(fabs(tt - t) < 000000001, "i2x2_mod_norm_sqr error");
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_sym_eigen ---\n"); }
//      if (drandom() < 0.2) 
//        { /* Test with a diagonal matrix: */
//          throw_diag_matrix(&A);
//        }
//      else
//        { /* Test with a general symmetric matrix: */
//          throw_symmetric_matrix(&A);
//        }
//      i2x2_sym_eigen(&A, &a, &B);
//      // if (verbose)
//      //   { fprintf(stderr, "testing i2x2_sym_eigen:\n");
//      //     fprintf(stderr, "A = ");
//      //     i2x2_print(stderr, &A);
//      //     fputc('\n', stderr);
//      //     fprintf(stderr, "eigenvalues = ");
//      //     i2_print(stderr, &a);
//      //     fputc('\n', stderr);
//      //     fprintf(stderr, "eigenvectors = ");
//      //     i2x2_print(stderr, &B);
//      //     fputc('\n', stderr);
//      //   }
//      /* Check order of eigenvalues: */
//      for (i = 1; i < N; i++)
//        { affirm(a.c[i-1] >= a.c[i], "i2x2_sym_eigen error: order"); }
//      /* Check whether {B} is orthonormal: */
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { double sum = 0.0;
//              for (k = 0; k < N; k++) 
//                { sum += B.c[k][i]*B.c[k][j]; }
//              double val = (i == j ? 1.0 : 0.0);
//              in_check_eps(val,sum,0.000000001 * fabs(sum), NO, NO, 
//                "i2x2_sym_eigen error: not orthormal"
//              );
//            }
//        }
//      /* Check whether {B} is right-handed: */
//      rr = i2x2_det(&B);
//      affirm(fabs(rr - 1.0) < 000000001, "i2x2_sym_eigen error: not right-handed");
//      /* Check whether {A = B'*DIAG(e)*B}: */
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { double sum = 0.0;
//              for (k = 0; k < N; k++) 
//                { sum += B.c[k][i]*a.c[k]*B.c[k][j]; }
//              in_check_eps(A.c[i][j],sum,0.000000001 * fabs(sum), NO, NO, 
//                "i2x2_sym_eigen error: decomp"
//              );
//            }
//        }
//  
//      if (verbose) { fprintf(stderr, "--- i2x2_print ---\n"); }
//      if (verbose)
//        { throw_matrix (&A);
//          fprintf(stderr, "A = ");
//          i2x2_print(stderr, &A);
//          fputc('\n', stderr);
//        }
//    }  
//  
//  void throw_matrix(i2x2_t *m)
//    {
//      int32_t i, j;
//      i2_t a;
//      for (i = 0; i < N; i++)
//        { i2_throw_cube(trad, &a);
//          for (j = 0; j < N; j++) { m->c[i][j] = a.c[j]; }
//        }
//    }
//  
//  void throw_diag_matrix(i2x2_t *m)
//    {
//      int32_t i, j;
//      for (i = 0; i < N; i++)
//        { for (j = 0; j < N; j++)
//            { m->c[i][j] = (i == j ? 2*drandom() - 1.0 : 0.0); }
//        }
//    }
//  
//  void throw_symmetric_matrix(i2x2_t *m)
//    {
//      int32_t i, j;
//      for (i = 0; i < N; i++)
//        { /* Note: {j} runs to {i} not {N-1}! */
//          for (j = 0; j <= i; j++)
//            { m->c[i][j] = 2*drandom() - 1.0;
//              if (j != i) { m->c[j][i] = m->c[i][j]; }
//            }
//        }
//    }
