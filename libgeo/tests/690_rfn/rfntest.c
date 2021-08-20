/* rfntest --- test program for rfn.h, rfmxn.h  */
/* Last edited on 2021-08-20 17:26:37 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include <affirm.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <flt.h>
#include <bool.h>

#include <rn.h>
#include <rn_test_tools.h>

#include <rfn.h>
#include <rfmxn.h>
#include <rfn_check.h>

#define NO NULL

/* Internal prototypes */

int32_t main (int32_t argc, char **argv);
void test_rfn(int32_t verbose);
void test_rfmxn(int32_t verbose);
void throw_matrix(int32_t m, int32_t n, float *Amn);
void throw_LT_matrix(int32_t m, float *Lmm);
void print_matrix(FILE *wr, int32_t m, int32_t n, float *Amn);

void check_ortho_matrix(int32_t n, float M[]);
  /* Checks whether the {n x n} matrix {M} orthonormal; that is
    whether the rows are pairwise orthogonal and have length 1. */

int32_t main (int32_t argc, char **argv)
  { int32_t i;
    srand(1993);
    srandom(1993);
    for (i = 0; i < 100; i++) test_rfn(i <= 3);
    for (i = 0; i < 100; i++) test_rfmxn(i <= 3);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }

void test_rfn (int32_t verbose)
  { float *a, *b, *c, *d, *e, *para, *perp;
    double r, s, t;
    float tf; /* rf, sf,  */
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;
    int32_t maxsize = (verbose ? 5 : 10);
    int32_t n = rand()/(RAND_MAX/maxsize) + 1;  /* Vector size. */

    fprintf(stderr, "test_rfn:  n = %d\n", n);
    
    a = rfn_alloc(n);
    b = rfn_alloc(n);
    c = rfn_alloc(n);
    d = rfn_alloc(n);
    e = rfn_alloc(n);
    para = rfn_alloc(n);
    perp = rfn_alloc(n);

    if (verbose) { fprintf(stderr, "--- rfn_zero ---\n"); }
    rfn_zero(n, a);
    for (i = 0; i < n; i++)
      { rn_check_eq(a[i],0.0, &i, NO, "rfn_zero error"); }

    if (verbose) { fprintf(stderr, "--- rfn_all ---\n"); }
    rfn_all(n, 3.14f, a);
    for (i = 0; i < n; i++)
      { rn_check_eq(a[i],3.14f, &i, NO, "rfn_all error"); }

    if (verbose) { fprintf(stderr, "--- rfn_axis ---\n"); }
    for (k = 0; k < n; k++)
      { rfn_axis(n, k, a);
        for (i = 0; i < n; i++)
          { rn_check_eq(a[i],(i == k ? 1.0 : 0.0), &i, NO, "rfn_axis error"); }
      }

    if (verbose) { fprintf(stderr, "--- rfn_throw_cube ---\n"); }
    /* Should check uniformity... */
    rfn_throw_cube(n, a);
    for (i = 0; i < n; i++)
      { affirm((n == 1) || (a[i] != a[(i+1)%n]), "rfn_throw_cube probable error(1)"); 
        /* Check whether there are more than 8 nonzero bits: */
        float vv = a[i]*256.0f;
        if (vv == (float)floor(vv))
          { 
            fprintf(stderr, "!! {rfn_throw_cube} gives less than 8 leading bits: %24.16e %24.16e\n", a[i], vv);
          }
        affirm((a[i] > -1.0) && (a[i] < 1.0), "rfn_throw_cube error(2)"); 
      }

    if (verbose) { fprintf(stderr, "--- rfn_throw_dir ---\n"); }
    /* Should check uniformity... */
    rfn_throw_dir(n, a);
    /* Check variation: */
    for (i = 0; i < n; i++) { affirm((n == 1) || (a[i] != a[(i+1)%n]), "rfn_throw_dir error(1)"); }
    /* Check whether the norm is 1: */
    rr = 0;
    for (i = 0; i < n; i++) { double ai = a[i]; rr += ai*ai; }
    rn_check_eps(1, rr, 0.000001*n*rr, NO, NO, "rfn_throw_dir error (2)");

    if (verbose) { fprintf(stderr, "--- rfn_throw_ball ---\n"); }
    /* Should check uniformity... */
    rfn_throw_ball(n, a);
    /* Check variation: */
    for (i = 0; i < n; i++) { affirm((n == 1) || (a[i] != a[(i+1)%n]), "rfn_throw_ball error(1)"); }
    /* Check whether the norm is at most 1: */
    rr = 0;
    for (i = 0; i < n; i++) { double ai = a[i]; rr += ai*ai; }
    demand(rr <= 1 + 0.000001*rr, "rfn_throw_ball error (2)");

    if (verbose) { fprintf(stderr, "--- rfn_add ---\n"); }
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    rfn_add(n, a, b, d);
    for (i = 0; i < n; i++)
      { rn_check_eq(d[i],a[i] + b[i], &i, NO, "rfn_add error"); }

    if (verbose) { fprintf(stderr, "--- rfn_sub ---\n"); }
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    rfn_sub(n, a, b, d);
    for (i = 0; i < n; i++)
      { rn_check_eq(d[i],a[i] - b[i], &i, NO, "rfn_sub error"); }

    if (verbose) { fprintf(stderr, "--- rfn_neg ---\n"); }
    rfn_throw_cube(n, a);
    rfn_neg(n, a, d);
    for (i = 0; i < n; i++)
      { double dxi = -a[i];
        rn_check_eq(d[i],dxi, &i, NO, "rfn_neg error"); }

    if (verbose) { fprintf(stderr, "--- rfn_scale ---\n"); }
    s = drandom();
    rfn_throw_cube(n, a);
    rfn_scale(n, s, a, d);
    for (i = 0; i < n; i++)
      { float dxi = (float)(s*a[i]);
        fprintf(stderr, "%5d %24.16e %24.16e %24.16e %24.16e\n", i, s, a[i], d[i], dxi);
        rn_check_eq(d[i], dxi, &i, NO, "rfn_scale error(1)");
      }

    if (verbose) { fprintf(stderr, "--- rfn_shift ---\n"); }
    s = drandom();
    rfn_throw_cube(n, a);
    rfn_shift(n, s, a, d);
    for (i = 0; i < n; i++)
      { float dxi = (float)(s + a[i]);
        fprintf(stderr, "%5d %24.16e %24.16e %24.16e %24.16e\n", i, s, a[i], d[i], dxi);
        rn_check_eq(d[i],dxi, &i, NO, "rfn_scale error(1)");
      }

    if (verbose) { fprintf(stderr, "--- rfn_mix ---\n"); }
    s = drandom();
    t = drandom();
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    rfn_mix(n, s, a, t, b, d);
    for (i = 0; i < n; i++)
      { float dxi = (float)(s * a[i] + t * b[i]);
        rn_check_eq(d[i],dxi, &i, NO, "rfn_mix error");
      }

    if (verbose) { fprintf(stderr, "--- rfn_mix_in ---\n"); }
    s = drandom();
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    for (i = 0; i < n; i++) { c[i] = b[i]; }
    rfn_mix_in(n, s, a, c);
    for (i = 0; i < n; i++)
      { float dxi = (float)(b[i] + s * a[i]);
        rn_check_eq(c[i],dxi, &i, NO, "rfn_mix_in error");
      }

    if (verbose) { fprintf(stderr, "--- rfn_weigh ---\n"); }
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    rfn_weigh(n, a, b, d);
    for (i = 0; i < n; i++)
      { float dxi = (float)(((double)a[i]) * b[i]);
        rn_check_eq(d[i],dxi, &i, NO, "rfn_weigh error");
      }

    if (verbose) { fprintf(stderr, "--- rfn_unweigh ---\n"); }
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    rfn_unweigh(n, a, b, d);
    for (i = 0; i < n; i++)
      { float dxi = (float)(((double)a[i]) / b[i]);
        rn_check_eq(d[i],dxi, &i, NO, "rfn_unweigh error");
      }

    if (n >= 2)
      { if (verbose) { fprintf(stderr, "--- rfn_rot_axis ---\n"); }
        { rfn_throw_cube(n, a);
          int32_t i = int32_abrandom(0, n-1);
          int32_t j = int32_abrandom(0, n-2); if (j >= i) { j++; }
          double ang = 2.1*M_PI*drandom();
          rfn_rot_axis(n, a, i, j, ang, d);
          rfn_check_rot_axis(n, a, i, j, ang, d, "rfn_rot_axis error");
        }
      }

    if (verbose) { fprintf(stderr, "--- rfn_sum ---\n"); }
    rfn_throw_cube(n, a);
    s = rfn_sum(n, a);
    ss = 0.0;
    for (i = 0; i < n; i++) { ss += a[i]; }
    rr = rfn_L_inf_norm(n, a);
    rn_check_eps(s,ss, 0.000001*rr, NO, NO, "rfn_sum error");

    if (verbose) { fprintf(stderr, "--- rfn_norm, rfn_norm_sqr, rfn_L_inf_norm ---\n"); }
    rfn_throw_cube(n, a);
    r = rfn_norm(n, a);
    s = rfn_norm_sqr(n, a);
    t = rfn_L_inf_norm(n, a);
    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < n; i++)
      { double ai = fabs(a[i]);
        ss += ai*ai; 
        if (ai > tt) { tt = ai; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000001*n*rr, NO, NO, "rfn_norm error");
    rn_check_eps(s,ss,0.000001*n*ss, NO, NO, "rfn_norm_sqr error");
    rn_check_eq(t,tt, NO, NO, "rfn_L_inf_norm error");


    if (verbose) { fprintf(stderr, "--- rfn_dist, rfn_dist_sqr, rfn_L_inf_dist ---\n"); }
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    r = rfn_dist(n, a, b);
    s = rfn_dist_sqr(n, a, b);
    t = rfn_L_inf_dist(n, a, b);

    ss = 0.0;
    tt = 0.0;
    for (i = 0; i < n; i++)
      { double di = fabs(((double)a[i]) - b[i]);
        ss += di*di; 
        if (di > tt) { tt = di; }
      }
    rr = sqrt(ss);
    rn_check_eps(r,rr,0.000001*n*rr, NO, NO, "rfn_dist error");
    rn_check_eps(s,ss,0.000001*n*ss, NO, NO, "rfn_dist_sqr error");
    rn_check_eq(t,tt, NO, NO, "rfn_L_inf_dist error");

    if (verbose) { fprintf(stderr, "--- rfn_dir, rfn_L_inf_dir ---\n"); }
    rfn_throw_cube(n, a);
    r = rfn_dir(n, a, b);
    s = rfn_L_inf_dir(n, a, d);
    ss = rfn_norm(n, a);
    tf = rfn_L_inf_norm(n, a);
    for (i = 0; i < n; i++)
      { float bi = (float)(a[i]/ss);
        rn_check_eps(b[i], bi, 0.000001*n*ss, NO, NO, "rfn_dir error");
        float di = (float)(a[i]/((double)tf));
        rn_check_eps(d[i], di, 0.000001*n*tt, NO, NO, "rfn_L_inf_dir error");
      }

    if (verbose) { fprintf(stderr, "--- rfn_dot, rfn_cos, rfn_sin, rfn_angle ---\n"); }
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    r = rfn_dot(n, a, b);
    double S = rfn_sin(n, a, b);
    double C = rfn_cos(n, a, b);
    double A = rfn_angle(n, a, b);
    mag = sqrt(rfn_dot(n, a,a)*rfn_dot(n, b,b));
    rr = 0.0;
    for (i = 0; i < n; i++) { rr += a[i]*b[i]; }
    double CC = rr/(rfn_norm(n, a)*rfn_norm(n, b));
    rn_check_eps(r,rr,0.000001*n*mag, NO, NO, "rfn_dot error(1)");
    rn_check_eps(C,CC,0.000001, NO, NO, "rfn_cos error(1)");
    for (i = 0; i < n; i++) { d[i] = a[i]; }
    rfn_mix_in(n, -rr/rfn_norm_sqr(n, b), b, d);
    double SS = rfn_norm(n, d)/rfn_norm(n, a);
    rn_check_eps(S,SS,0.000001, NO, NO, "rfn_sin error(1)");
    double AA = atan2(SS, CC);
    rn_check_eps(A,AA,0.000001, NO, NO, "rfn_angle error(1)");
    for (i = 0; i < n; i++)
      { rfn_axis(n, i, a);
        for (j = 0; j < n; j++)
          { rfn_axis(n, j, b);
            r = rfn_dot(n, a, b);
            s = rfn_sin(n, a, b);
            t = rfn_cos(n, a, b);
            rr = (i == j ? 1.0 : 0.0);
            rn_check_eq(r,rr, &i, &j, "rfn_dot error(2)");
            rn_check_eq(t,rr, &i, &j, "rfn_dot error(3)");
            rn_check_eq(s,1.0 - rr, &i, &j, "rfn_dot error(4)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rfn_cross ---\n"); }
    { float **z = (float **)malloc((n-1)*sizeof(float *));
      double *magz = rn_alloc(n);
      /* Test on basis vectors: */
      for (i = 0; i < n; i++)
        { float sign = ((n-1)*i % 2 == 0 ? +1.0 : -1.0);
          for (k = 0; k < n-1; k++)
            { int32_t ik = (i + k) % n; 
              z[k] = rfn_alloc(n);
              rfn_axis(n, ik, z[k]);
            }
          { int32_t in1 = (i + n-1) % n; rfn_axis(n, in1, b); }
          rfn_cross(n, z, c);
          for (j = 0; j < n; j++)
            { float cxj = sign*b[j];
              rn_check_eq(c[j],cxj, &i, &j, "rfn_cross error(x)");
            }
        }
      /* Test on random vectors: */
      mag = 1.0;
      for (k = 0; k < n-1; k++)
        { z[k] = rfn_alloc(n);
          rfn_throw_cube(n, z[k]);
          { magz[k] = rfn_norm(n, z[k]); mag *= magz[k]; }
        }
      rfn_cross(n, z, c);
      for (k = 0; k < n-1; k++)
        { r = rfn_dot(n, z[k], c);
          rn_check_eps(r,0.0,0.000001*n*mag*magz[k], NO, NO, "rfn_cross error(1)");
        }
      for (k = 0; k < n-1; k++) { free(z[k]); }
      free(z); free(magz);
    }
    
    if (verbose) { fprintf(stderr, "--- rfn_det ---\n"); }
    { float **z = (float **)malloc(n*sizeof(float *));
      double *magz = rn_alloc(n);
      /* Test on basis vectors: */
      for (i = 0; i < n; i++)
        { float sign = ((n-1)*i % 2 == 0 ? +1.0 : -1.0);
          for (k = 0; k < n; k++)
            { int32_t ik = (i + k) % n; 
              z[k] = rfn_alloc(n);
              rfn_axis(n, ik, z[k]);
            }
          r = rfn_det(n, z);
          rn_check_eq(r,sign, &i, NO, "rfn_det error(2)");
        }
      /* Test on random vectors: */
      mag = 1.0;
      for (k = 0; k < n; k++)
        { z[k] = rfn_alloc(n);
          rfn_throw_cube(n, z[k]);
          { magz[k] = rfn_norm(n, z[k]); mag *= magz[k]; }
        }
      r = rfn_det(n, z);
      rfn_cross(n, z, c);
      rr = rfn_dot(n, c, z[n-1]);
      rn_check_eps(r,rr,0.000001*n*mag, NO, NO, "rfn_det error(1)");

      for (k = 0; k < n; k++) { free(z[k]); }
      free(z); free(magz);
    }

    if (verbose) { fprintf(stderr, "--- rfn_decomp ---\n"); }
    rfn_throw_cube(n, a);
    rfn_throw_cube(n, b);
    r = rfn_decomp(n, a, b, para, perp);
    rr = rfn_dot(n, a, b)/rfn_norm_sqr(n, b);  
    rn_check_eps(r,rr,0.000001*n*(fabs(r) + fabs(rr)), NO, NO, "rfn_decomp error(1)");
    rfn_add(n, para, perp, c);
    s = rfn_dist(n, a, c);
    rn_check_eps(s,0.0,0.000001*n*rfn_norm(n, a), NO, NO, "rfn_decomp error(2)");
    s = rfn_dot(n, perp, b);
    rn_check_eps(s,0.0,0.000001*n*rfn_norm(n, b), NO, NO, "rfn_decomp error(3)");
    t = rfn_dot(n, para, perp);
    rn_check_eps(t,0.0,0.000001*n*rfn_norm(n, a), NO, NO, "rfn_decomp error(4)");

    if (verbose) { fprintf(stderr, "--- rfn_mirror ---\n"); }
    rfn_throw_cube(n, a);
    rfn_throw_dir(n, d);
    r = rfn_mirror(n, a, d, b);
    /* The dot products must be equal and opposite: */
    double tol = 1.0e-6*rfn_norm(n,a);
    r = rfn_dot(n, a, d);  
    s = rfn_dot(n, b, d);
    rn_check_eps(r,-s,tol, NO, NO, "rfn_mirror error(1)");
    /* Compute the average {c} of {a} and {b}: */
    rfn_mix(n, 0.5, a, 0.5, b, c);
    /* The average of {a} and {b} must be orthogonal to {d}: */
    t = rfn_dot(n, c, d);
    rn_check_eps(t,0,tol, NO, NO, "rfn_mirror error(2)");

    if (verbose) { fprintf(stderr, "--- rfn_print ---\n"); }
    if (verbose)
      { rfn_throw_cube (n, a);
        fprintf(stderr, "a = ");
        rfn_print(stderr, n, a);
        fputc('\n', stderr);
      }
      
    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(para);
    free(perp);
  }
  
void test_rfmxn(int32_t verbose)
  { float *Amn, *Bmn, *Cmn, *Amp, *Bpn, *Amm, *Bmm, *Cmm;
    float *am, *bm, *cm, *an, *bn, *cn;
    double r, s, t;
    float rf, sf;
    double rr, ss, tt;
    double mag;
    int32_t i, j, k;
    int32_t maxsize = (verbose ? 5 : 10);
    int32_t m = rand()/(RAND_MAX/maxsize) + 1;  /* Number of rows. */
    int32_t n = rand()/(RAND_MAX/maxsize) + 1;  /* Number of columns. */
    int32_t p = rand()/(RAND_MAX/maxsize) + 1;  /* Middle dimension for rfmxn_mul. */

    fprintf(stderr, "test_rfmxn:  m = %d  n = %d  p = %d\n", m, n, p);
    
    Amn = rfmxn_alloc(m, n);
    Bmn = rfmxn_alloc(m, n);
    Cmn = rfmxn_alloc(m, n);

    Amp = rfmxn_alloc(m, p);
    Bpn = rfmxn_alloc(p, n);

    Amm = rfmxn_alloc(m, m);
    Bmm = rfmxn_alloc(m, m);
    Cmm = rfmxn_alloc(m, m);

    am = rfn_alloc(m);
    bm = rfn_alloc(m);
    cm = rfn_alloc(m);
    
    an = rfn_alloc(n);
    bn = rfn_alloc(n);
    cn = rfn_alloc(n);

    if (verbose) { fprintf(stderr, "--- rfmxn_zero, rfmxn_ident ---\n"); }
    rfmxn_zero(m, n, Amn);
    rfmxn_ident(m, n, Bmn);
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { rn_check_eq(Amn[n*i + j],0.0, &i, &j, "rfmxn_zero error"); 
            rn_check_eq(Bmn[n*i + j],(i == j ? 1.0 : 0.0), &i, &j, "rfmxn_ident error");
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_copy ---\n"); }
    throw_matrix(m, n, Amn);
    rfmxn_copy(m, n, Amn, Bmn);
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { rn_check_eq(Amn[n*i + j],Bmn[n*i + j], &i, &j, "rfmxn_copy error"); }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_get_row, rfmxn_set_row ---\n"); }
    throw_matrix(m, n, Amn);
    rfn_throw_cube(n, an);
    rfmxn_copy(m, n, Amn, Bmn);
    rfn_copy(n, an, bn);
    k = int32_abrandom(0, m-1);
    rfmxn_get_row(m, n, Bmn, k, cn);
    rfmxn_set_row(m, n, Bmn, k, bn);
    for (j = 0; j < n; j++)
      { /* Check whether {rfmxn_get_row} copied the row correctly: */
        double cjObs = cn[j];
        double cjExp = Amn[k*n + j];
        rn_check_eq(cjObs,cjExp, &k, &j, "rfmxn_get_row error (1)"); 
        /* Check whether {rfmxn_set_row} modified the vector arg: */
        double bjObs = bn[j];
        double bjExp = an[j];
        rn_check_eq(bjObs,bjExp, &k, &j, "rfmxn_set_row error (1)");
        for (i = 0; i < m; i++)
          { double BijObs = Bmn[i*n + j];
            double BijExp = (i == k ? an[j] : Amn[i*n + j]);
            rn_check_eq(BijObs,BijExp, &i, &j, "rfmxn_set_row error (2)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_get_col, rfmxn_set_col ---\n"); }
    throw_matrix(m, n, Amn);
    rfn_throw_cube(m, am);
    rfmxn_copy(m, n, Amn, Bmn);
    rfn_copy(m, am, bm);
    k = int32_abrandom(0, n-1);
    rfmxn_get_col(m, n, Bmn, k, cm);
    rfmxn_set_col(m, n, Bmn, k, bm);
    for (i = 0; i < m; i++)
      { /* Check whether {rfmxn_get_col} copied the column correctly: */
        double ciObs = cm[i];
        double ciExp = Amn[i*n + k];
        rn_check_eq(ciObs,ciExp, &i, &k, "rfmxn_get_col error (1)"); 
        /* Check whether {rfmxn_set_col} modified the vector arg: */
        double biObs = bm[i];
        double biExp = am[i];
        rn_check_eq(biObs,biExp, &i, &k, "rfmxn_set_col error (1)");
        for (j = 0; j < n; j++)
          { double BijObs = Bmn[i*n + j];
            double BijExp = (j == k ? am[i] : Amn[i*n + j]);
            rn_check_eq(BijObs,BijExp, &i, &j, "rfmxn_set_col error (2)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_scale ---\n"); }
    s = drandom();
    throw_matrix(m, n, Amn);
    rfmxn_scale(m, n, s, Amn, Cmn);
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { float zij = (float)(s*Amn[n*i + j]);
            rn_check_eq(Cmn[n*i + j],zij, &i, &j, "rfmxn_scale error(1)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_mix ---\n"); }
    s = drandom();
    t = drandom();
    throw_matrix(m, n, Amn);
    throw_matrix(m, n, Bmn);
    rfmxn_mix(m, n, s, Amn, t, Bmn, Cmn);
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { float zij = (float)(s*Amn[n*i + j] + t*Bmn[n*i + j]);
            rn_check_eq(Cmn[n*i + j],zij, &i, &j, "rfmxn_mix error(1)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_rel_diff ---\n"); }
    s = drandom();
    t = drandom();
    throw_matrix(m, n, Amn);
    throw_matrix(m, n, Bmn);
    rfmxn_rel_diff(m, n, Amn, Bmn, Cmn);
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { float dij = (float)rel_diff(Amn[n*i + j], Bmn[n*i + j]);
            rn_check_eq(Cmn[n*i + j], dij, &i, &j, "rfmxn_rel_diff error(1)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_map_row, rfmxn_map_col ---\n"); }
    throw_matrix(m, n, Amn);
    rfn_throw_cube(m, am);
    rfn_throw_cube(n, an);
    rfmxn_map_row(m, n, am, Amn, bn);
    rfmxn_map_col(m, n, Amn, an, bm);
    rfn_zero(m, cm);
    rfn_zero(n, cn);
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { cn[j] += am[i] * Amn[n*i + j];
            cm[i] += Amn[n*i + j] * an[j];
          }
      }
    r = rfn_dist(m, bm, cm);
    rn_check_eps(r, 0.0, 0.000001*n*rfn_norm(m, cm), NO, NO, "rfn_map_row error");
    s = rfn_dist(n, bn, cn);
    rn_check_eps(s, 0.0, 0.000001*n*rfn_norm(n, cn), NO, NO, "rfn_map_col error");

    if (verbose) { fprintf(stderr, "--- rfmxn_mul ---\n"); }
    throw_matrix(m, p, Amp);
    throw_matrix(p, n, Bpn);
    rfmxn_mul(m, p, n, Amp, Bpn, Cmn);
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double sum = 0.0;
            for (k = 0; k < p; k++) { sum += Amp[p*i + k]*Bpn[n*k + j]; }
            rn_check_eps(Cmn[n*i + j], (float)sum, 0.000001*n*fabs(sum), NO, NO, 
              "rfmxn_mul error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_mul_tr ---\n"); }
    throw_matrix(m, n, Amn);
    throw_matrix(p, n, Bpn);
    rfmxn_mul_tr(m, p, n, Amn, Bpn, Amp);
    for (i = 0; i < m; i++)
      { for (j = 0; j < p; j++)
          { double sum = 0.0;
            for (k = 0; k < n; k++) { sum += ((double)Amn[n*i + k])*Bpn[n*j + k]; }
            rn_check_eps(Amp[p*i + j], (float)sum,0.0000003*(m+n) * fabs(sum), &i, &j, 
              "rfmxn_mul_tr error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_tr_mul ---\n"); }
    throw_matrix(m, p, Amp);
    throw_matrix(m, n, Amn);
    rfmxn_tr_mul(m, p, n, Amp, Amn, Bpn);
    for (i = 0; i < p; i++)
      { for (j = 0; j < n; j++)
          { double sum = 0.0;
            for (k = 0; k < m; k++) { sum += ((double)Amp[k*p + i])*Amn[k*n + j]; }
            rn_check_eps(Bpn[i*n + j], (float)sum, 0.00001 * fabs(sum), &i, &j, 
              "rfmxn_tr_mul error"
            );
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_det ---\n"); }
    throw_matrix(m, m, Amm);
    for (i = 0; i < m; i++)
      { int32_t k = (i + 1) % m;
        for (j = 0; j < m; j++)
          { /* Check for linearity */
            rf = frandom();
            Amm[m*i + j] = rf;
            rr = rfmxn_det(m, Amm);

            sf = frandom();
            Amm[m*i + j] = sf;
            ss = rfmxn_det(m, Amm);

            t = frandom();
            Amm[m*i + j] = (float)(rf*(1-t) + sf*t);
            tt = rfmxn_det(m, Amm);
            mag = fabs(rr) + fabs(ss) + fabs(tt);
            rn_check_eps(tt, (rr*(1-t) + ss*t),0.0000003*(m+n) * mag, NO, NO, 
              "rfmxn_det error(1)"
            );
          }

        if (m > 1) 
          { /* Row swap test: */
            k = (i + 1) % m;
            r = rfmxn_det(m, Amm);
            for (j = 0; j < m; j++)
              { float *Aij = &(Amm[m*i + j]);
                float *Akj = &(Amm[m*k + j]);
                { float t = *Aij; *Aij = *Akj; *Akj = t; }
              }
            rr = rfmxn_det(m, Amm);
            mag = fabs(r) + fabs(rr);
            rn_check_eps(r, (-rr), 0.000001*n*mag, NO, NO, "rfmxn_det error(2)");

            /* Col swap test: */
            r = rfmxn_det(m, Amm);
            for (j = 0; j < m; j++)
              { float *Aji = &(Amm[m*j + i]);
                float *Ajk = &(Amm[m*j + k]);
                { float t = *Aji; *Aji = *Ajk; *Ajk = t; }
              }
            rr = rfmxn_det(m, Amm);
            mag = fabs(r) + fabs(rr);
            rn_check_eps(r, (-rr), 0.000001*n*mag, NO, NO, "rfmxn_det error(3)");
          }
      }

    if (verbose) { fprintf(stderr, "--- rfmxn_cof ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rfmxn_cof not tested\n"); }
    
    if (verbose) { fprintf(stderr, "--- rfmxn_adj ---\n"); }
    if (verbose) { fprintf(stderr, "!! warning: rfmxn_adj not tested\n"); }
    
    if (verbose) { fprintf(stderr, "--- rfmxn_inv ---\n"); }
    throw_matrix(m, m, Amm);
    r = rfmxn_det(m, Amm);
    s = rfmxn_inv(m, Amm, Bmm);
    rfmxn_mul(m, m, m, Amm, Bmm, Cmm);
    for (i = 0; i < m; i++)
      { for (j = 0; j < m; j++)
          { double val = (i == j ? 1.0 : 0.0);
            rn_check_eps(Cmm[m*i + j], val, 0.000007*m, NO, NO, "rfmxn_inv error");
          }
      }
    rn_check_eps(r,s,0.000001, NO, NO, "rfmxn_inv/rfmxn_det error");
    
    if (verbose) { fprintf(stderr, "--- rfmxn_inv_full ---\n"); }
    s = rfmxn_inv_full(m, Amm, Bmm);
    rfmxn_mul(m, m, m, Amm, Bmm, Cmm);
    for (i = 0; i < m; i++)
      { for (j = 0; j < m; j++)
          { double val = (i == j ? 1.0 : 0.0);
            rn_check_eps(Cmm[m*i + j], val, 0.000007*m, NO, NO, "rfmxn_inv_full error");
          }
      }
    rn_check_eps(r,s,0.000001, NO, NO, "rfmxn_inv-full/rfmxn_det error");

    if (verbose) { fprintf(stderr, "--- rfmxn_cholesky ---\n"); }
    { float *Lmm = Amm, *Jmm = Cmm;
      throw_LT_matrix(m, Jmm);
      for(i = 0; i < m; i++) { Jmm[m*i + i] = fabsf(Jmm[m*i + i]) + 1.0f; }
      rfmxn_mul_tr(m, m, m, Jmm, Jmm, Bmm);
      rfmxn_cholesky(m, Bmm, Lmm);
      for (i = 0; i < m; i++)
        { for (j = 0; j < m; j++)
            { double eps = 0.000001*m;
              rn_check_eps(Lmm[m*i + j], Jmm[m*i + j], eps, NO, NO, 
                "rfmxn_cholesky error"
              );
            }
        }
    }
      
    if (verbose) { fprintf(stderr, "--- rfmxn_LT_inv_map_row, rfmxn_LT_inv_map_col ---\n"); }
    { float *Lmm = Amm;
      int32_t by_row; 
      throw_LT_matrix(m, Lmm);
      (void)rfmxn_inv(m, Lmm, Bmm);
      rfn_throw_cube(m, am);
      for (by_row = 0; by_row < 2; by_row++)
        { char *msg;
          if (by_row) 
            { rfmxn_LT_inv_map_row(m, am, Lmm, bm);
              rfmxn_map_row(m, m, am, Bmm, cm);
              msg = "rfmxn_LT_inv_map_row error";
            }
          else
            { rfmxn_LT_inv_map_col(m, Lmm, am, bm);
              rfmxn_map_col(m, m, Bmm, am, cm);
              msg = "rfmxn_LT_inv_map_col error";
            }
          r = rfn_dist(m, bm, cm);
          double eps =  0.000005*m * rfn_norm(m, cm);
          rn_check_eps(r, 0.0, eps, NO, NO, msg);
        }
    }
      
    if (verbose) { fprintf(stderr, "--- rfmxn_LT_pos_div ---\n"); }
    { float *Lmm = Amm, *Jmm = Bmm;
      /* We must use "Amn" as if it were "n x m" not "m x n". */
      float *Anm = Amn, *Bnm = Bmn, *Cnm = Cmn; 
      throw_LT_matrix(m, Lmm);
      (void)rfmxn_inv(m, Lmm, Jmm);
      throw_matrix(n, m, Anm);
      rfmxn_LT_pos_div(n, m, Anm, Lmm, Cnm);
      rfmxn_mul(n, m, m, Anm, Jmm, Bnm);
      for (i = 0; i < n; i++)
        { for (j = 0; j < m; j++)
            { double eps = 0.00002*m*(fabs(Bnm[m*i + j]) + fabs(Cnm[m*i + j]) + 1.0e-200);
              rn_check_eps(Cnm[m*i + j],Bnm[m*i + j],eps, NO, NO, 
                "rfmxn_LT_pos_div error"
              );
            }
        }
    }
    
    if (verbose) { fprintf(stderr, "--- rfmxn_LT_pre_div ---\n"); }
    { float *Lmm = Amm, *Jmm = Bmm;
      throw_LT_matrix(m, Lmm);
      (void)rfmxn_inv(m, Lmm, Jmm);
      throw_matrix(m, n, Amn);
      rfmxn_LT_pre_div(m, n, Lmm, Amn, Cmn);
      rfmxn_mul(m, m, n, Jmm, Amn, Bmn);
      for (i = 0; i < m; i++)
        { for (j = 0; j < n; j++)
            { double eps = 0.00002*m*(fabs(Bmn[n*i + j]) + fabs(Cmn[n*i + j]) + 1.0e-200);
              rn_check_eps(Cmn[n*i + j], Bmn[n*i + j], eps, NO, NO, 
                "rfmxn_LT_pre_div error"
              );
            }
        }
    }
     
    if (verbose) { fprintf(stderr, "--- rfmxn_norm,rfmxn_norm_sqr ---\n"); }
    throw_matrix(m,n,Amn);
    s = rfmxn_norm_sqr(m,n,Amn);
    r = rfmxn_norm(m,n,Amn);
    ss = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { double Aij = Amn[n*i + j]; 
            ss += Aij*Aij;
          }
      }
    affirm(ss >= 0, "rfmxn_norm_sqr error");
    affirm(fabs(ss - s) < 0.000001, "rfmxn_norm_sqr error");
    rr = sqrt(ss);
    affirm(fabs(rr - r) < 0.000001, "rfmxn_norm error");

    if (verbose) { fprintf(stderr, "--- rfmxn_norm,rfmxn_norm_sqr,rfmxn_mod_norm ---\n"); }
    throw_matrix(m,m,Amm);
    t = rfmxn_mod_norm_sqr(m,Amm);
    tt = 0;
    for (i = 0; i < m; i++)
      { for (j = 0; j < m; j++)
          { double Aij = Amm[m*i + j]; 
            double Dij = (i == j ? Aij - 1 : Aij);
            tt += Dij*Dij;
          }
      }
    affirm(tt >= 0, "rfmxn_mod_norm_sqr error");
    affirm(fabs(tt - t) < 0.000001, "rfmxn_mod_norm_sqr error");

    if (verbose) { fprintf(stderr, "--- rfmxn_print ---\n"); }
    if (verbose)
      { throw_matrix (m, n, Amn);
        fprintf(stderr, "Amn = ");
        rfmxn_print(stderr, m, n, Amn);
        fputc('\n', stderr);
      }
    
    free(Amn); free(Bmn); free(Cmn);
    free(Amp); free(Bpn);
    free(Amm); free(Bmm); free(Cmm);
    free(am); free(bm); free(cm);
    free(an); free(bn); free(cn);
  }  

void throw_matrix(int32_t m, int32_t n, float *Amn)
  { int32_t i, j;
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++) 
          { Amn[n*i + j] = 2.0f * frandom() - 1.0f; }
      }
  }

void throw_LT_matrix(int32_t m, float *Lmm)
  { int32_t i, j;
    for (i = 0; i < m; i++)
      { for (j = 0; j < m; j++) 
          { Lmm[m*i + j] = (j <= i ? 2.0f * frandom() - 1.0f : 0.0f); }
      }
  }

void check_ortho_matrix(int32_t n, float M[])
  { 
    double tol = 1.0e-12;
    int32_t i0, i1, j;
    for (i0 = 0; i0 < n; i0++)
      { for (i1 = 0; i1 <= i0; i1++)
          { /* Compute dot product {sObs} of rows {i0} and {i1}: */
            double sObs = 0;
            for (j = 0; j < n; j++) { sObs += ((double)M[i0*n + j])*M[i1*n + j]; }
            /* Check dot product against expected value {sExp}: */
            double sExp = (i0 == i1 ? 1 : 0);
            rn_check_eps(sExp, sObs, tol, &i0, &i1, "rfmxn_throw_ortho error");
          }
      }
  }

void print_matrix(FILE *wr, int32_t m, int32_t n, float *Amn)
  { int32_t i, j;
    fprintf(wr, "%d × %d matrix\n", m, n);
    for (i = 0; i < m; i++)
      { for (j = 0; j < n; j++)
          { fprintf(wr, " %22.16e", Amn[i*n + j]); }
        fprintf(wr, "\n");
      }
    fprintf(wr, "\n");
  }
