/*Last edited on 2024-11-22 02:01:10 by stolfi */
/*
  Based on VectorN.mg, created  95-02-27 by J. Stolfi.
  Last edited by stolfi 
*/

/* !!! We should use Kahan's summation for all scalar products. */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <rn.h>
#include <rmxn.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_elim.h>

void rn_zero (uint32_t n, double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = 0.0; }
  }

void rn_all (uint32_t n, double x, double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = x; }
  }

void rn_axis (uint32_t n, uint32_t i, double r[])
  { 
    affirm(i < n, "rn_axis: bad index");
    for (int32_t j = 0; j < n; j++) { r[j] = 0.0; }
    r[i] = 1.0;
  }

void rn_copy (uint32_t n, double a[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = a[i]; }
  }

void rn_add (uint32_t n, double a[], double b[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = a[i] + b[i]; }
  }

void rn_sub (uint32_t n, double a[], double b[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = a[i] - b[i]; }
  }

void rn_neg (uint32_t n, double a[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = - a[i]; }
  }

void rn_scale (uint32_t n, double s, double a[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = s * a[i]; }
  }

void rn_shift (uint32_t n, double s, double a[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = s + a[i]; }
  }

void rn_mix (uint32_t n, double s, double a[], double t, double b[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = s * a[i] + t * b[i]; }
  }

void rn_mix_in (uint32_t n, double s, double a[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] += s * a[i]; }
  }

void rn_weigh (uint32_t n, double a[], double w[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = a[i] * w[i]; }
  }

void rn_unweigh (uint32_t n, double a[], double w[], double r[])
  { 
    for (int32_t i = 0; i < n; i++) { r[i] = a[i] / w[i]; }
  }

void rn_rot_axis (uint32_t n, double a[], uint32_t i, uint32_t j, double ang, double r[])
  {
    affirm(i < n, "rn_rot_axis: bad index {i}");
    affirm(j < n, "rn_rot_axis: bad index {j}");
    affirm(i != j, "rn_rot_axis: axes not distinct");
    double c = cos(ang);
    double s = sin(ang);
    double x = + c*a[i] - s*a[j];
    double y = + s*a[i] + c*a[j];
    for (int32_t k = 0; k < n; k++) { r[k] = (k == i ? x : (k == j ? y : a[k])); }
  }

double rn_sum (uint32_t n, double a[])
  {
    double sum = 0;
    for (int32_t i = 0; i < n; i++) { sum += a[i]; }
    return sum;
  }

double rn_norm (uint32_t n, double a[])
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double sum = 0.0;
    for (int32_t i = 0; i < n; i++) { double ai = a[i]; sum += ai*ai; }
    return sqrt(sum);
  }

double rn_norm_sqr (uint32_t n, double a[])
  { double sum = 0.0;
    for (int32_t i = 0; i < n; i++) { double ai = a[i]; sum += ai*ai; }
    return sum;
  }

double rn_L_inf_norm (uint32_t n, double a[])
  { double mag = 0.0;
    for (int32_t i = 0; i < n; i++) 
      { double mi = fabs(a[i]); if (mi > mag) { mag = mi; } }
    return mag;
  }

double rn_dist (uint32_t n, double a[], double b[])
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double sum = 0.0;
    for (int32_t i = 0; i < n; i++) { double di = a[i] - b[i]; sum += di*di; }
    return sqrt(sum);
  }

double rn_dist_sqr (uint32_t n, double a[], double b[])
  { double sum = 0.0;
    for (int32_t i = 0; i < n; i++) { double di = (a[i] - b[i]); sum += di*di; }
    return sum;
  }

double rn_L_inf_dist (uint32_t n, double a[], double b[])
  { double mag = 0.0;
    for (int32_t i = 0; i < n; i++) 
      { double mi = fabs(a[i] - b[i]); if (mi > mag) { mag = mi; } }
    return mag;
  }

double rn_dir (uint32_t n, double a[], double r[])
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double d = rn_norm(n, a);
    for (int32_t i = 0; i < n; i++) { r[i] = a[i]/d; }
    return d;
  }

double rn_L_inf_dir (uint32_t n, double a[], double r[])
  { double mag = rn_L_inf_norm(n, a);
    for (int32_t i = 0; i < n; i++) { r[i] = a[i]/mag; }
    return mag;
  }

double rn_dot (uint32_t n, double a[], double b[])
  { double sum = 0.0;
    for (int32_t i = 0; i < n; i++) { sum += a[i]*b[i]; }
    return sum;
  }

double rn_cos (uint32_t n, double a[], double b[])
  { double aa = 0.0;
    double bb = 0.0;
    double ab = 0.0;
    for (int32_t i = 0; i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi; ab += ai*bi;
      }
    return ab/(sqrt(aa)*sqrt(bb));
  }

double rn_sin (uint32_t n, double a[], double b[])
  { 
    /* Compute {aa = a*a, bb = b*b}: */ 
    double aa = 0.0;
    double bb = 0.0;
    for (int32_t i = 0; i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi;
      }
    /* Normalize {a, b}, compute {d = a-b}, {s = a+b}, {dd = |d|^2}, {ss = |s|^2}: */ 
    { double na = 1.0/sqrt(aa);
      double nb = 1.0/sqrt(bb);
      double dd = 0.0;
      double ss = 0.0;

      for (int32_t i = 0; i < n; i++)
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

double rn_angle (uint32_t n, double a[], double b[])
  { /* Compute {aa = a*a, bb = b*b}: */ 
    double aa = 0.0;
    double bb = 0.0;
    for (int32_t i = 0; i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi;
      }
    /* Normalize {a, b}, compute {d = a-b}, {s = a+b}, {dd = |d|^2}, {ss = |s|^2}: */ 
    { double na = 1.0/sqrt(aa);
      double nb = 1.0/sqrt(bb);
      double dd = 0.0;
      double ss = 0.0;

      for (int32_t i = 0; i < n; i++)
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

void rn_cross (uint32_t n, double *a[], double r[])
  { uint32_t nn1 = (n-1)*n;
    double *C = rmxn_alloc(n-1,n);
    { uint32_t t = 0;
      for (int32_t i = 0; i < n-1; i++) 
        { double *ai = a[i]; 
          for (int32_t j = 0; j < n; j++) { C[t] = ai[j]; t++; }
        }
    }
    gsel_triangularize(n-1, n, C, TRUE, 0.0);
    gsel_diagonalize(n-1, n, C);
    /* If {det(C)} is not zero, set {d = det(C)}, {izer = -1}.
      Else set {izer} to the first row {i} such that {C[i,i]} is zero, 
      and  {d} to the determinant of {C} excluding column {izer}. */
    double d = 1.0; 
    int32_t izer = -1; 
    { uint32_t t = 0;
      for (int32_t i = 0; i < n-1; i++)
        { if (C[t] == 0.0)
            { if (izer < 0) { izer = (int32_t)i; t++; } else { d = 0.0; break; } }
          d *= C[t]; t += n+1;
        }
    }
    if (izer < 0)
      { /* The main diagonal of {C} is non-zero, and possibly column {n-1}. */
        r[n-1] = d;
        /* For {i=0..n-2}, set {r[i] = -d*(C[i,n-1]/C[i,i])}: */  
        uint32_t s = nn1-1; 
        uint32_t t = nn1-2;
        for (int32_t i = (int32_t)n-2; i >= 0; i--) 
          { r[i] = -d*C[s]/C[t]; s -= n; t -= n+1; }
      }
    else
      { /* Set {r[izer] = (-1)^(n-1-izer)*d}, all other elems to 0; */
        assert((izer >= 0) && (izer < n-1)); 
        for (int32_t i = 0; i < n; i++) { r[i] = 0.0; }
        r[izer] = ((n - 1 - (uint32_t)izer) % 2 == 0 ? d : -d);
      }
    free(C);
  }

double rn_det (uint32_t n, double *a[])
  { uint32_t n2 = n*n;
    double *C = rmxn_alloc(n,n);
    double d;
    uint32_t t = 0;
    for (int32_t i = 0; i < n; i++) 
      { double *ai = a[i]; 
        for (int32_t j =  0; j < n; j++) { C[t] = ai[j]; t++; }
      }
    gsel_triangularize(n, n, C, FALSE, 0.0);
    d = 1.0;
    for (int32_t k = 0; k < n2; k += n+1) { d *= C[k]; }
    free(C);
    return d;
  }

double rn_decomp (uint32_t n, double a[], double u[], double para[], double perp[])
  { double sau = 0.0;
    double suu = 0.0;
    for (int32_t i = 0; i < n; i++) 
      { double ai = a[i]; double ui = u[i];
        suu += ui*ui;
        sau += ai*ui;
      }
    if (sau == 0.0) 
      { for (int32_t i = 0; i < n; i++) 
          { if (para != NULL) { para[i] = 0.0; }
            if (perp != NULL) { perp[i] = a[i]; } 
          } 
        return 0.0;
      }
    else
      { double c = sau / suu;
        for (int32_t i = 0; i < n; i++) 
          { double pi = c * u[i]; 
            if (para != NULL) { para[i] = pi; }
            if (perp != NULL) { perp[i] = a[i] - pi; }
          } 
        return c;
      }
  }

double rn_mirror (uint32_t n, double a[], double u[], double r[])
  { double sau = 0.0;
    for (int32_t i = 0; i < n; i++) 
      { double ai = a[i]; double ui = u[i];
        sau += ai*ui;
      }
    if (sau != 0.0) 
      { double c = 2*sau;
        for (int32_t i = 0; i < n; i++) { r[i] = a[i] - c*u[i]; }
      }
    return sau;
  }

uint32_t rn_remove_zeros(uint32_t n, double a[], double r[])
  { uint32_t m = 0;
    for (int32_t i = 0; i < n; i++)
      { if (a[i] != 0) { r[m] = a[i]; m++; } }
    return m;
  }
    
uint32_t rn_insert_zeros(uint32_t n, double a[], double b[], double r[])
  { uint32_t m = 0;
    for (int32_t i = 0; i < n; i++)
      { if (a[i] != 0) { r[i] = b[i]; m++; } else { r[i] = 0; } }
    return m;
  }

double rn_rad_rel_max_diff (uint32_t n, double a[], double b[], double rad[])
  {
    double dmax = 0.0;
    for (int32_t i = 0; i < n; i++)
      { double di = fabs(a[i] - b[i]);
        if (rad != NULL)
          { double ri = rad[i];
            demand(isfinite(ri) && (ri >= 0), "invalid search radius");
            if (ri != 0)
              { di /= ri; }
            else
              { demand(di == 0, "null {rad[i]} with unequal {a[i],b[i]}"); }
          }
        dmax = fmax(dmax, di);
      }
    return dmax;
  }

double rn_rad_rel_dist_sqr (uint32_t n, double a[], double b[], double rad[])
  {
    double d2 = 0.0;
    for (int32_t i = 0; i < n; i++)
      { double di = a[i] - b[i];
        if (rad != NULL)
          { double ri = rad[i];
            demand(isfinite(ri) && (ri >= 0), "invalid search radius");
            if (ri != 0)
              { di /= ri; }
            else
              { demand(di == 0, "null {rad[i]} with unequal {a[i],b[i]}"); }
          }
        d2 += di*di;
      }
    return d2;
  }


void rn_throw_cube (uint32_t n, double r[])
  { 
    for (int32_t i = 0; i < n; i++)
      { r[i] = 2.0 * drandom() - 1.0; }
  }

void rn_throw_normal (uint32_t n, double r[])
  { 
    for (int32_t i = 0; i < n; i++) 
      { r[i] = dgaussrand(); }
  }

void rn_throw_dir (uint32_t n, double r[])
  { if (n == 0) 
      { return; }
    else if (n == 1) 
      { r[0] = (drandom() <= 0.5 ? +1 : -1); }
    else
      { /* Generate a nonzero Gaussian random vector: */
        double r2;
        do
          { rn_throw_normal(n, r);
            /* Discard if too close to origin: */
            r2 = 0.0;
            for (int32_t i = 0; i < n; i++) { double ci = r[i]; r2 += ci*ci; }
          }
        while (r2 < 1.0e-5);
        /* Normalize to unit length: */
        double m = sqrt(r2);
        for (int32_t i = 0; i < n; i++) { r[i] /= m; }
      }
  }

void rn_throw_ball (uint32_t n, double r[])
  { if (n == 0) 
      { return; }
    else if (n == 1) 
      { r[0] = 2.0 * drandom() - 1.0; }
    else
      { /* Generate a random unit vector: */
        rn_throw_dir(n, r);
        /* Generate a random radius {z} with density {n*r^{n-1}}: */
        double z = drandom();
        if (z > 0) { z = exp(log(z)/n); }
        /* Scale {r[0..n-1]} by {z}: */
        for (int32_t i = 0; i < n; i++) { r[i] *= z; }
      }
  }

double rn_abs_rel_diff(uint32_t n, double a[], double b[], double abs_tol, double rel_tol)
  { double max_error = 0.0;
    for (int32_t i = 0; i < n; i++)
      { double error = abs_rel_diff(a[i], b[i], abs_tol, rel_tol);
        if (fabs(error) > max_error) { max_error = fabs(error); }
      }
    return max_error;
  }

void rn_print (FILE *wr, uint32_t n, double a[])
  { rn_gen_print(wr, n, a, NULL, NULL, NULL, NULL); }

void rn_gen_print 
  ( FILE *wr, uint32_t n, double a[], 
    char *fmt, 
    char *lp, char *sep, char *rp
  )
  { if (fmt == NULL) { fmt = "%16.8e"; }
    if (lp == NULL) { lp = "("; }
    if (sep == NULL) { sep = " "; }
    if (rp == NULL) { rp = ")"; }
    fputs(lp, wr);
    for (int32_t i = 0; i < n; i++)
      { if (i > 0) { fputs(sep, wr); }
        fprintf(wr, fmt, a[i]);
      }
    fputs(rp, wr);
  }

void rn_rad_rel_print(FILE *wr, uint32_t n, double a[], double rad[])
  {
    rn_gen_print(wr, n, a, "%12.7f", "[", " ", "]\n");
    rn_gen_print(wr, n, rad, "%12.7f", "[", " ", "]\n");
    double ar[n];
    rn_unweigh(n, a, rad, ar);
    rn_gen_print(wr, n, ar, "%12.7f", "[", " ", "]\n");
  }

double *rn_alloc(uint32_t n)
  { 
    double *p = talloc(n, double);
    return p;
  }

