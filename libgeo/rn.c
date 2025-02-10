/*Last edited on 2025-02-05 15:37:27 by stolfi */
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
#include <gausol_triang.h>
#include <gausol_solve.h>

void rn_zero (uint32_t n, double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = 0.0; }
  }

void rn_all (uint32_t n, double x, double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = x; }
  }

void rn_axis (uint32_t n, uint32_t i, double r[])
  { 
    affirm(i < n, "rn_axis: bad index");
    for (uint32_t j = 0;  j < n; j++) { r[j] = 0.0; }
    r[i] = 1.0;
  }

void rn_copy (uint32_t n, double a[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = a[i]; }
  }

void rn_add (uint32_t n, double a[], double b[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = a[i] + b[i]; }
  }

void rn_sub (uint32_t n, double a[], double b[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = a[i] - b[i]; }
  }

void rn_neg (uint32_t n, double a[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = - a[i]; }
  }

void rn_scale (uint32_t n, double s, double a[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = s * a[i]; }
  }

void rn_shift (uint32_t n, double s, double a[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = s + a[i]; }
  }

void rn_mix (uint32_t n, double s, double a[], double t, double b[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = s * a[i] + t * b[i]; }
  }

void rn_mix_in (uint32_t n, double s, double a[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] += s * a[i]; }
  }

void rn_weigh (uint32_t n, double a[], double w[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = a[i] * w[i]; }
  }

void rn_unweigh (uint32_t n, double a[], double w[], double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) { r[i] = a[i] / w[i]; }
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
    for (uint32_t k = 0;  k < n; k++) { r[k] = (k == i ? x : (k == j ? y : a[k])); }
  }

double rn_sum (uint32_t n, double a[])
  {
    double sum = 0;
    for (uint32_t i = 0;  i < n; i++) { sum += a[i]; }
    return sum;
  }

double rn_norm (uint32_t n, double a[])
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double sum = 0.0;
    for (uint32_t i = 0;  i < n; i++) { double ai = a[i]; sum += ai*ai; }
    return sqrt(sum);
  }

double rn_norm_sqr (uint32_t n, double a[])
  { double sum = 0.0;
    for (uint32_t i = 0;  i < n; i++) { double ai = a[i]; sum += ai*ai; }
    return sum;
  }

double rn_L_inf_norm (uint32_t n, double a[])
  { double mag = 0.0;
    for (uint32_t i = 0;  i < n; i++) 
      { double mi = fabs(a[i]); if (mi > mag) { mag = mi; } }
    return mag;
  }

double rn_dist (uint32_t n, double a[], double b[])
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double sum = 0.0;
    for (uint32_t i = 0;  i < n; i++) { double di = a[i] - b[i]; sum += di*di; }
    return sqrt(sum);
  }

double rn_dist_sqr (uint32_t n, double a[], double b[])
  { double sum = 0.0;
    for (uint32_t i = 0;  i < n; i++) { double di = (a[i] - b[i]); sum += di*di; }
    return sum;
  }

double rn_L_inf_dist (uint32_t n, double a[], double b[])
  { double mag = 0.0;
    for (uint32_t i = 0;  i < n; i++) 
      { double mi = fabs(a[i] - b[i]); if (mi > mag) { mag = mi; } }
    return mag;
  }

double rn_dir (uint32_t n, double a[], double r[])
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double d = rn_norm(n, a);
    for (uint32_t i = 0;  i < n; i++) { r[i] = (d == 0 ? NAN : a[i]/d); }
    return d;
  }

double rn_L_inf_dir (uint32_t n, double a[], double r[])
  { double mag = rn_L_inf_norm(n, a);
    for (uint32_t i = 0;  i < n; i++) { r[i] = (mag == 0 ? NAN : a[i]/mag); }
    return mag;
  }

double rn_dot (uint32_t n, double a[], double b[])
  { double sum = 0.0;
    for (uint32_t i = 0;  i < n; i++) { sum += a[i]*b[i]; }
    return sum;
  }

double rn_cos (uint32_t n, double a[], double b[])
  { double aa = 0.0;
    double bb = 0.0;
    double ab = 0.0;
    for (uint32_t i = 0;  i < n; i++)
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
    for (uint32_t i = 0;  i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi;
      }
    /* Normalize {a, b}, compute {d = a-b}, {s = a+b}, {dd = |d|^2}, {ss = |s|^2}: */ 
    { double na = 1.0/sqrt(aa);
      double nb = 1.0/sqrt(bb);
      double dd = 0.0;
      double ss = 0.0;

      for (uint32_t i = 0;  i < n; i++)
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
    for (uint32_t i = 0;  i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi;
      }
    /* Normalize {a, b}, compute {d = a-b}, {s = a+b}, {dd = |d|^2}, {ss = |s|^2}: */ 
    { double na = 1.0/sqrt(aa);
      double nb = 1.0/sqrt(bb);
      double dd = 0.0;
      double ss = 0.0;

      for (uint32_t i = 0;  i < n; i++)
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
  { demand(n >= 1, "invalid {n}");
    uint32_t m = n-1;
    double *C = rmxn_alloc(m, n);
    double *Cij = C;
    for (uint32_t i = 0; i < m; i++) 
      { double *ai = a[i];
        for (uint32_t j = 0; j < n; j++) { (*Cij) = ai[j]; Cij++; }
      }
    /* !!! Using full pivoting. Would partial pivoting be better? !!! */
    uint32_t prow[m], pcol[n]; 
    double tiny = 1.0e-180;
    double det; /* Det of {P00}. */
    uint32_t rank;
    gausol_triang_reduce(m,prow, n,pcol,C, 0,NULL, tiny, &det, &rank);
    gausol_triang_diagonalize(m,prow, n,pcol,C, 0,NULL, rank, tiny);
    if (rank == m)
      { /* Rows are linearly independent, good. */
        demand(isfinite(det), "invalid {det} from triangulation");
        demand(det != 0, "triangulation {det} should not be zero");
        /* Make {r[0..n-1]} orthogonal to {C[t,*]} for each {t}: */
        uint32_t pcol_m = pcol[m]; assert(pcol_m < n);
        r[pcol_m] = det;
        for (uint32_t t = 0; t < m; t++)
          { uint32_t prow_t = prow[t]; assert(prow_t < m);
            uint32_t pcol_t = pcol[t]; assert(pcol_t < n);
            double Ctt = C[prow_t*n + pcol_t];
            assert(isfinite(Ctt) && Ctt != 0);
            double Ctm = C[prow_t*n + pcol_m];
            double Rt = -det*Ctm/Ctt;
            demand(isfinite(Rt), "overflow");
            r[pcol_t] = Rt;
          }
      }
    else
      { /* Rows are linearly DEpendent. */
        for (uint32_t j = 0; j < n; j++) { r[j] = 0.0; }
      }
    free(C);
  }

double rn_det (uint32_t n, double *a[])
  { double *C = rmxn_alloc(n,n);
    for (uint32_t i = 0;  i < n; i++) 
      { double *ai = a[i]; 
        double *Ci = &(C[i*n]); 
        for (uint32_t j = 0;  j < n; j++) { Ci[j] = ai[j]; }
      }
    /* !!! Using full pivoting. Would partial pivoting be better? !!! */
    uint32_t prow[n], pcol[n]; /* Permutation matrices for Gaussian eleimination. */
    double det;
    uint32_t rank;
    double tiny = 1.0e-180;
    gausol_triang_reduce(n, prow, n, pcol, C, 0, NULL, tiny, &det, &rank);
    assert(! isnan(det));
    free(C);
    return (rank == n ? det : 0);
  }

double rn_decomp (uint32_t n, double a[], double u[], double para[], double perp[])
  { double sau = 0.0;
    double suu = 0.0;
    for (uint32_t i = 0;  i < n; i++) 
      { double ai = a[i]; double ui = u[i];
        suu += ui*ui;
        sau += ai*ui;
      }
    if (sau == 0.0) 
      { for (uint32_t i = 0;  i < n; i++) 
          { if (para != NULL) { para[i] = 0.0; }
            if (perp != NULL) { perp[i] = a[i]; } 
          } 
        return 0.0;
      }
    else
      { double c = sau / suu;
        for (uint32_t i = 0;  i < n; i++) 
          { double pi = c * u[i]; 
            if (para != NULL) { para[i] = pi; }
            if (perp != NULL) { perp[i] = a[i] - pi; }
          } 
        return c;
      }
  }

double rn_mirror (uint32_t n, double a[], double u[], double r[])
  { double sau = 0.0;
    for (uint32_t i = 0;  i < n; i++) 
      { double ai = a[i]; double ui = u[i];
        sau += ai*ui;
      }
    if (sau != 0.0) 
      { double c = 2*sau;
        for (uint32_t i = 0;  i < n; i++) { r[i] = a[i] - c*u[i]; }
      }
    return sau;
  }

uint32_t rn_remove_zeros(uint32_t n, double a[], double r[])
  { uint32_t m = 0;
    for (uint32_t i = 0;  i < n; i++)
      { if (a[i] != 0) { r[m] = a[i]; m++; } }
    return m;
  }
    
uint32_t rn_insert_zeros(uint32_t n, double a[], double b[], double r[])
  { uint32_t m = 0;
    for (uint32_t i = 0;  i < n; i++)
      { if (a[i] != 0) { r[i] = b[i]; m++; } else { r[i] = 0; } }
    return m;
  }

double rn_rad_rel_max_diff (uint32_t n, double a[], double b[], double rad[])
  {
    double dmax = 0.0;
    for (uint32_t i = 0;  i < n; i++)
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
    for (uint32_t i = 0;  i < n; i++)
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
    for (uint32_t i = 0;  i < n; i++)
      { r[i] = 2.0 * drandom() - 1.0; }
  }

void rn_throw_normal (uint32_t n, double r[])
  { 
    for (uint32_t i = 0;  i < n; i++) 
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
            for (uint32_t i = 0;  i < n; i++) { double ci = r[i]; r2 += ci*ci; }
          }
        while (r2 < 1.0e-5);
        /* Normalize to unit length: */
        double m = sqrt(r2);
        for (uint32_t i = 0;  i < n; i++) { r[i] /= m; }
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
        for (uint32_t i = 0;  i < n; i++) { r[i] *= z; }
      }
  }

double rn_abs_rel_diff(uint32_t n, double a[], double b[], double abs_tol, double rel_tol)
  { double max_error = 0.0;
    for (uint32_t i = 0;  i < n; i++)
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
    for (uint32_t i = 0;  i < n; i++)
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

