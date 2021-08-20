/* See {rfn.h}. */
/* Last edited on 2021-08-20 16:16:40 by stolfi */
/* Based on VectorN.mg, created  95-02-27 by J. Stolfi. */

/* !!! We should use Kahan's summation for all scalar products. */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include <rn.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_elim.h>

#include <rfn.h>

void rfn_zero (int32_t n, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = 0.0; }
  }

void rfn_all (int32_t n, float x, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = x; }
  }

void rfn_axis (int32_t n, int32_t i, float *r)
  { int32_t j;
    affirm((i >= 0) && (i < n), "rfn_axis: bad index");
    for (j = 0; j < n; j++) { r[j] = 0.0; }
    r[i] = 1.0;
  }

void rfn_copy (int32_t n, float *a, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = a[i]; }
  }

void rfn_add (int32_t n, float *a, float *b, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(((double)a[i]) + b[i]); }
  }

void rfn_sub (int32_t n, float *a, float *b, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = a[i] - b[i]; }
  }

void rfn_neg (int32_t n, float *a, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = - a[i]; }
  }

void rfn_scale (int32_t n, double s, float *a, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(s * a[i]); }
  }

void rfn_shift (int32_t n, double s, float *a, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(s + a[i]); }
  }

void rfn_mix (int32_t n, double s, float *a, double t, float *b, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(s * a[i] + t * b[i]); }
  }

void rfn_mix_in (int32_t n, double s, float *a, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(r[i] + s * a[i]); }
  }

void rfn_weigh (int32_t n, float *a, float *w, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(((double)a[i]) * w[i]); }
  }

void rfn_unweigh (int32_t n, float *a, float *w, float *r)
  { int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(((double)a[i]) / w[i]); }
  }

void rfn_rot_axis (int32_t n, float *a, int32_t i, int32_t j, double ang, float *r)
  {
    affirm((i >= 0) && (i < n), "rn_rot_axis: bad index {i}");
    affirm((j >= 0) && (j < n), "rn_rot_axis: bad index {j}");
    affirm(i != j, "rn_rot_axis: axes not distinct");
    double c = cos(ang);
    double s = sin(ang);
    float x = (float)(+ c*a[i] - s*a[j]);
    float y = (float)(+ s*a[i] + c*a[j]);
    for (int32_t k = 0; k < n; k++) { r[k] = (k == i ? x : (k == j ? y : a[k])); }
  }

double rfn_sum (int32_t n, float *a)
  { int32_t i;
    double sum = 0;
    for (i = 0; i < n; i++) { sum += a[i]; }
    return sum;
  }

double rfn_norm (int32_t n, float *a)
  { /* Don't worry about overflow. */
    /* Client should use {rfn_L_inf_dir} first if that is a problem. */
    double sum = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) { double ai = a[i]; sum += ai*ai; }
    return sqrt(sum);
  }

double rfn_norm_sqr (int32_t n, float *a)
  { double sum = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) { double ai = a[i]; sum += ai*ai; }
    return sum;
  }

float rfn_L_inf_norm (int32_t n, float *a)
  { float mag = 0.0;
    for (int32_t i = 0; i < n; i++) 
      { float mi = fabsf(a[i]); if (mi > mag) { mag = mi; } }
    return mag;
  }

double rfn_dist (int32_t n, float *a, float *b)
  { /* Don't worry about overflow. */
    /* Client should use {rfn_L_inf_dir} first if that is a problem. */
    double sum = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) { double di = a[i] - b[i]; sum += di*di; }
    return sqrt(sum);
  }

double rfn_dist_sqr (int32_t n, float *a, float *b)
  { double sum = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) { double di = (a[i] - b[i]); sum += di*di; }
    return sum;
  }

double rfn_L_inf_dist (int32_t n, float *a, float *b)
  { double mag = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) 
      { double mi = fabs(a[i] - b[i]); if (mi > mag) { mag = mi; } }
    return mag;
  }

double rfn_dir (int32_t n, float *a, float *r)
  { /* Don't worry about overflow. */
    /* Client should use {rfn_L_inf_dir} first if that is a problem. */
    double d = rfn_norm(n, a);
    int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(a[i]/d); }
    return d;
  }

float rfn_L_inf_dir (int32_t n, float *a, float *r)
  { float mag = rfn_L_inf_norm(n, a);
    int32_t i;
    for (i = 0; i < n; i++) { r[i] = (float)(((double)a[i])/mag); }
    return mag;
  }

double rfn_dot (int32_t n, float *a, float *b)
  { double sum = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) { sum += ((double)a[i])*b[i]; }
    return sum;
  }

double rfn_cos (int32_t n, float *a, float *b)
  { double aa = 0.0;
    double bb = 0.0;
    double ab = 0.0;
    int32_t i;
    for (i = 0; i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi; ab += ai*bi;
      }
    return ab/(sqrt(aa)*sqrt(bb));
  }

double rfn_sin (int32_t n, float *a, float *b)
  { int32_t i;
    /* Compute {aa = a*a, bb = b*b}: */ 
    double aa = 0.0;
    double bb = 0.0;
    for (i = 0; i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi;
      }
    /* Normalize {a, b}, compute {d = a-b}, {s = a+b}, {dd = |d|^2}, {ss = |s|^2}: */ 
    { double na = 1.0/sqrt(aa);
      double nb = 1.0/sqrt(bb);
      double dd = 0.0;
      double ss = 0.0;

      for (i = 0; i < n; i++)
        { double ai = na * a[i]; 
          double bi = nb * b[i]; 
          double di = ai - bi;
          double si = ai + bi;
          dd += di*di; ss += si*si;
        }
      /* Now we must have {sqrt(dd) = 2*sin(t/2)}, {sqrt(ss) = 2*cos(t/2)}. */
      /* We must have also {dd + ss = 4}, modulo roundoff errors. */
      /* Use formula {sin(t) = 2*sin(t/2)*cos(t/2) = sqrt(dd)*sqrt(ss)/2}. */
      /* Divide by {sqrt(dd + ss)/4} to compensate by roundoff (?). */ 
      return 2.0 * sqrt(dd*ss)/(dd + ss);
    }
  }

double rfn_angle (int32_t n, float *a, float *b)
  { int32_t i;
    /* Compute {aa = a*a, bb = b*b}: */ 
    double aa = 0.0;
    double bb = 0.0;
    for (i = 0; i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi;
      }
    /* Normalize {a, b}, compute {d = a-b}, {s = a+b}, {dd = |d|^2}, {ss = |s|^2}: */ 
    { double na = 1.0/sqrt(aa);
      double nb = 1.0/sqrt(bb);
      double dd = 0.0;
      double ss = 0.0;

      for (i = 0; i < n; i++)
        { double ai = na * a[i]; 
          double bi = nb * b[i]; 
          double di = ai - bi;
          double si = ai + bi;
          dd += di*di; ss += si*si;
        }
      /* Now we must have {sqrt(dd) = 2*sin(t/2)}, {sqrt(ss) = 2*cos(t/2)}. */
      return 2.0 * atan2(sqrt(dd), sqrt(ss));
    }
  }

void rfn_cross (int32_t n, float **a, float *r)
  { int32_t nn1 = (n-1)*n; int32_t t;
    double *C = (double *)notnull(malloc(nn1*sizeof(double)), "no mem for C");
    t = 0;
    for (int32_t i = 0; i < n-1; i++) 
      { float *ai = a[i]; 
        for (int32_t j = 0; j < n; j++) { C[t] = (double)ai[j]; t++; }
      }
    gsel_triangularize(n-1, n, C, TRUE, 0.0);
    gsel_diagonalize(n-1, n, C, TRUE);
    /* If {det(C)} is not zero, set {d = det(C)}, {izer = -1}.
      Else set {izer} to the first zero column in {C}, and 
      set {d} to the determinant of {C} excluding that column. */
    double d = 1.0; int32_t izer = -1;
    t = 0;
    for (int32_t i = 0; i < n-1; i++)
      { if (C[t] == 0.0)
          { if (izer < 0) { izer = i; t++; } else { d = 0.0; break; } }
        d *= C[t]; t += n+1;
      }
    if (izer < 0)
      { r[n-1] = (float)d;
        /* For {i=0..n-2}, set {r[i] = -d*(C[i,n-1]/C[i,i])}: */  
        int32_t s = nn1-1; 
        t = nn1-2;
        for (int32_t i = n-2; i >= 0; i--) 
          { r[i] = -(float)(d*C[s]/C[t]); s -= n; t -= n+1; }
      }
    else
      { /* Set {r[izer] = (-1)^(n-1-izer)*d}, all other elems to 0; */
        for (int32_t i = 0; i < n; i++) { r[i] = 0.0; }
        r[izer] = (float)((n - 1 - izer) % 2 == 0 ? d : -d);
      }
    free(C);
  }

double rfn_det (int32_t n, float **a)
  { int32_t n2 = n*n; int32_t t;
    t = 0;
    double *C = (double *)notnull(malloc(n2*sizeof(double)), "no mem for C");
    double d;
    for (int32_t i = 0; i < n; i++) 
      { float *ai = a[i]; 
        for (int32_t j = 0; j < n; j++) { C[t] = (double)ai[j]; t++; }
      }
    gsel_triangularize(n, n, C, FALSE, 0.0);
    d = 1.0;
    for (int32_t k = 0; k < n2; k += n+1) { d *= C[k]; }
    free(C);
    return d;
  }

double rfn_decomp (int32_t n, float *a, float *u, float *para, float *perp)
  { double sau = 0.0;
    double suu = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) 
      { double ai = a[i]; double ui = u[i];
        suu += ui*ui;
        sau += ai*ui;
      }
    if (sau == 0.0) 
      { for (i = 0; i < n; i++) 
          { if (para != NULL) { para[i] = 0.0; }
            if (perp != NULL) { perp[i] = a[i]; } 
          } 
        return 0.0;
      }
    else
      { double c = sau / suu;
        for (i = 0; i < n; i++) 
          { double pi = c * u[i]; 
            if (para != NULL) { para[i] = (float)pi; }
            if (perp != NULL) { perp[i] = (float)(a[i] - pi); }
          } 
        return c;
      }
  }

double rfn_mirror (int32_t n, float *a, float *u, float *r)
  { double sau = 0.0;
    int32_t i;
    for (i = 0; i < n; i++) 
      { double ai = a[i]; double ui = u[i];
        sau += ai*ui;
      }
    if (sau != 0.0) 
      { double c = 2*sau;
        for (i = 0; i < n; i++) { r[i] = (float)(a[i] - c*u[i]); }
      }
    return sau;
  }

void rfn_throw_cube (int32_t n, float *r)
  { int32_t i;
    for (i = 0; i < n; i++)
      { r[i] = (float)(2.0 * frandom() - 1.0); }
  }

void rfn_throw_normal (int32_t n, float *r)
  { for (int32_t i = 0; i < n; i++) { r[i] = fgaussrand(); }
  }

void rfn_throw_dir (int32_t n, float *r)
  { if (n == 0) 
      { return; }
    else if (n == 1) 
      { r[0] = (drandom() <= 0.5 ? +1 : -1); }
    else
      { /* Generate a nonzero Gaussian random vector: */
        int32_t i;
        double r2;
        do
          { rfn_throw_normal(n, r);
            /* Discard if too close to origin: */
            r2 = 0.0;
            for (i = 0; i < n; i++) { double ci = r[i]; r2 += ci*ci; }
          }
        while (r2 < 1.0e-5);
        /* Normalize to unit length: */
        double m = sqrt(r2);
        for (i = 0; i < n; i++) { r[i] = (float)(r[i]/m); }
      }
  }

void rfn_throw_ball (int32_t n, float *r)
  { int32_t i;
    if (n == 0) 
      { return; }
    else if (n == 1) 
      { r[0] = (float)(2.0 * drandom() - 1.0); }
    else
      { /* Generate a random unit vector: */
        rfn_throw_dir(n, r);
        /* Generate a random radius {z} with density {n*r^{n-1}}: */
        double z = drandom();
        if (z > 0) { z = exp(log(z)/n); }
        /* Scale {r[0..n-1]} by {z}: */
        for (i = 0; i < n; i++) { r[i] = (float)(r[i]*z); }
      }
  }

double rfn_abs_rel_diff(int32_t n, float *a, float *b, double abs_tol, double rel_tol)
  { double max_error = 0.0;
    int32_t i;
    for (i = 0; i < n; i++)
      { double error = abs_rel_diff(a[i], b[i], abs_tol, rel_tol);
        if (fabs(error) > max_error) { max_error = fabs(error); }
      }
    return max_error;
  }

void rfn_print (FILE *f, int32_t n, float *a)
  { rfn_gen_print(f, n, a, NULL, NULL, NULL, NULL); }

void rfn_gen_print 
  ( FILE *f, int32_t n, float *a, 
    char *fmt, 
    char *lp, char *sep, char *rp
  )
  { int32_t i;
    if (fmt == NULL) { fmt = "%16.8e"; }
    if (lp == NULL) { lp = "("; }
    if (sep == NULL) { sep = " "; }
    if (rp == NULL) { rp = ")"; }
    fputs(lp, f);
    for (i = 0; i < n; i++)
      { if (i > 0) { fputs(sep, f); }
        fprintf(f, fmt, a[i]);
      }
    fputs(rp, f);
  }

float *rfn_alloc(int32_t n)
  { void *p = malloc(n*sizeof(float));
    affirm(p != NULL, "not enough memory");
    return (float *)p;
  }

