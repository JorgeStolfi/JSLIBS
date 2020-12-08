/*
  Last edited on 2014-03-24 23:32:23 by stolfilocal
  Based on VectorN.mg, created  95-02-27 by J. Stolfi.
  Last edited by stolfi 
*/

/* !!! We should use Kahan's summation for all scalar products. */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <rn.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <gauss_elim.h>

void rn_zero (int n, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = 0.0; }
  }

void rn_all (int n, double x, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = x; }
  }

void rn_axis (int n, int i, double *r)
  { int j;
    affirm((i >= 0) && (i < n), "rn_axis: bad index");
    for (j = 0; j < n; j++) { r[j] = 0.0; }
    r[i] = 1.0;
  }

void rn_copy (int n, double *a, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = a[i]; }
  }

void rn_add (int n, double *a, double *b, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = a[i] + b[i]; }
  }

void rn_sub (int n, double *a, double *b, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = a[i] - b[i]; }
  }

void rn_neg (int n, double *a, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = - a[i]; }
  }

void rn_scale (int n, double s, double *a, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = s * a[i]; }
  }

void rn_shift (int n, double s, double *a, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = s + a[i]; }
  }

void rn_mix (int n, double s, double *a, double t, double *b, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = s * a[i] + t * b[i]; }
  }

void rn_mix_in (int n, double s, double *a, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] += s * a[i]; }
  }

void rn_weigh (int n, double *a, double *w, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = a[i] * w[i]; }
  }

double rn_sum (int n, double *a)
  { int i;
    double sum = 0;
    for (i = 0; i < n; i++) { sum += a[i]; }
    return sum;
  }

double rn_norm (int n, double *a)
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double sum = 0.0;
    int i;
    for (i = 0; i < n; i++) { double ai = a[i]; sum += ai*ai; }
    return sqrt(sum);
  }

double rn_norm_sqr (int n, double *a)
  { double sum = 0.0;
    int i;
    for (i = 0; i < n; i++) { double ai = a[i]; sum += ai*ai; }
    return sum;
  }

double rn_L_inf_norm (int n, double *a)
  { double mag = 0.0;
    int i;
    for (i = 0; i < n; i++) 
      { double mi = fabs(a[i]); if (mi > mag) { mag = mi; } }
    return mag;
  }

double rn_dist (int n, double *a, double *b)
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double sum = 0.0;
    int i;
    for (i = 0; i < n; i++) { double di = a[i] - b[i]; sum += di*di; }
    return sqrt(sum);
  }

double rn_dist_sqr (int n, double *a, double *b)
  { double sum = 0.0;
    int i;
    for (i = 0; i < n; i++) { double di = (a[i] - b[i]); sum += di*di; }
    return sum;
  }

double rn_L_inf_dist (int n, double *a, double *b)
  { double mag = 0.0;
    int i;
    for (i = 0; i < n; i++) 
      { double mi = fabs(a[i] - b[i]); if (mi > mag) { mag = mi; } }
    return mag;
  }

double rn_dir (int n, double *a, double *r)
  { /* Don't worry about overflow. */
    /* Client should use {rn_L_inf_dir} first if that is a problem. */
    double d = rn_norm(n, a);
    int i;
    for (i = 0; i < n; i++) { r[i] = a[i]/d; }
    return d;
  }

double rn_L_inf_dir (int n, double *a, double *r)
  { double mag = rn_L_inf_norm(n, a);
    int i;
    for (i = 0; i < n; i++) { r[i] = a[i]/mag; }
    return mag;
  }

double rn_dot (int n, double *a, double *b)
  { double sum = 0.0;
    int i;
    for (i = 0; i < n; i++) { sum += a[i]*b[i]; }
    return sum;
  }

double rn_cos (int n, double *a, double *b)
  { double aa = 0.0;
    double bb = 0.0;
    double ab = 0.0;
    int i;
    for (i = 0; i < n; i++)
      { double ai = a[i]; 
        double bi = b[i]; 
        aa += ai*ai; bb += bi*bi; ab += ai*bi;
      }
    return ab/(sqrt(aa)*sqrt(bb));
  }

double rn_sin (int n, double *a, double *b)
  { int i;
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

double rn_angle (int n, double *a, double *b)
  { int i;
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

void rn_cross (int n, double **a, double *r)
  { int nn1 = (n-1)*n, t = 0, s, i, j, izer;
    double *C = (double *)notnull(malloc(nn1*sizeof(double)), "no mem for C");
    double d;
    for (i = 0; i < n-1; i++) 
      { double *ai = a[i]; 
        for (j = 0; j < n; j++) { C[t] = ai[j]; t++; }
      }
    gsel_triangularize(n-1, n, C, TRUE, 0.0);
    gsel_diagonalize(n-1, n, C, TRUE);
    /* If {det(C)} is not zero, set {d = det(C)}, {izer = -1}.
      Else set {izer} to the first zero column in {C}, and 
      set {d} to the determinant of {C} excluding that column. */
    d = 1.0; izer = -1; t = 0;
    for (i = 0; i < n-1; i++)
      { if (C[t] == 0.0)
          { if (izer < 0) { izer = i; t++; } else { d = 0.0; break; } }
        d *= C[t]; t += n+1;
      }
    if (izer < 0)
      { r[n-1] = d;
        /* For {i=0..n-2}, set {r[i] = -d*(C[i,n-1]/C[i,i])}: */  
        s = nn1-1; t = nn1-2;
        for (i = n-2; i >= 0; i--) 
          { r[i] = -d*C[s]/C[t]; s -= n; t -= n+1; }
      }
    else
      { /* Set {r[izer] = (-1)^(n-1-izer)*d}, all other elems to 0; */
        for (i = 0; i < n; i++) { r[i] = 0.0; }
        r[izer] = ((n - 1 - izer) % 2 == 0 ? d : -d);
      }
    free(C);
  }

double rn_det (int n, double **a)
  { int n2 = n*n, t = 0, i, j;
    double *C = (double *)notnull(malloc(n2*sizeof(double)), "no mem for C");
    double d;
    for (i = 0; i < n; i++) 
      { double *ai = a[i]; 
        for (j = 0; j < n; j++) { C[t] = ai[j]; t++; }
      }
    gsel_triangularize(n, n, C, FALSE, 0.0);
    d = 1.0;
    for (t = 0; t < n2; t += n+1) { d *= C[t]; }
    free(C);
    return d;
  }

double rn_decomp (int n, double *a, double *u, double *para, double *perp)
  { double sau = 0.0;
    double suu = 0.0;
    int i;
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
            if (para != NULL) { para[i] = pi; }
            if (perp != NULL) { perp[i] = a[i] - pi; }
          } 
        return c;
      }
  }

double rn_mirror (int n, double *a, double *u, double *r)
  { double sau = 0.0;
    int i;
    for (i = 0; i < n; i++) 
      { double ai = a[i]; double ui = u[i];
        sau += ai*ui;
      }
    if (sau != 0.0) 
      { double c = 2*sau;
        for (i = 0; i < n; i++) { r[i] = a[i] - c*u[i]; }
      }
    return sau;
  }

void rn_throw_cube (int n, double *r)
  { int i;
    for (i = 0; i < n; i++)
      { r[i] = 2.0 * drandom() - 1.0; }
  }

void rn_throw_normal (int n, double *r)
  { int i;
    for (i = 0; i < n; i++) { r[i] = dgaussrand(); }
  }

void rn_throw_dir (int n, double *r)
  { if (n == 0) 
      { return; }
    else if (n == 1) 
      { r[0] = (drandom() <= 0.5 ? +1 : -1); }
    else
      { /* Generate a nonzero Gaussian random vector: */
        int i;
        double r2;
        do
          { rn_throw_normal(n, r);
            /* Discard if too close to origin: */
            r2 = 0.0;
            for (i = 0; i < n; i++) { double ci = r[i]; r2 += ci*ci; }
          }
        while (r2 < 1.0e-5);
        /* Normalize to unit length: */
        double m = sqrt(r2);
        for (i = 0; i < n; i++) { r[i] /= m; }
      }
  }

void rn_throw_ball (int n, double *r)
  { int i;
    if (n == 0) 
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
        for (i = 0; i < n; i++) { r[i] *= z; }
      }
  }

double rn_abs_rel_diff(int n, double *a, double *b, double abs_tol, double rel_tol)
  { double max_error = 0.0;
    int i;
    for (i = 0; i < n; i++)
      { double error = abs_rel_diff(a[i], b[i], abs_tol, rel_tol);
        if (fabs(error) > max_error) { max_error = fabs(error); }
      }
    return max_error;
  }

void rn_print (FILE *f, int n, double *a)
  { rn_gen_print(f, n, a, NULL, NULL, NULL, NULL); }

void rn_gen_print 
  ( FILE *f, int n, double *a, 
    char *fmt, 
    char *lp, char *sep, char *rp
  )
  { int i;
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

double *rn_alloc(int n)
  { void *p = malloc(n*sizeof(double));
    affirm(p != NULL, "not enough memory");
    return (double *)p;
  }

double rn_ball_vol(double r, int d)
  {
    demand(d >= 0, "bad d");
    if (d == 0)
      { return 1; }
    else if (d == 1)
      { return 2*r; }
    else
      { return r*r * M_PI/(((double)d)/2) * rn_ball_vol(r, d-2); }
  }

double rn_ball_cap_vol_frac_pos(int d, double z)
  {
    /* fprintf(stderr, "    vol_frac_pos(%3d,%.10f)", d, z); */
    demand(d >= 0, "bad d");
    double f;
    if (z < -1.0) 
      { f = 0.0; }
    else if (z > 1.0) 
      { f = 1.0; }
    else if (d == 0)
      { /* The zero-dimensional sphere has just two points: */
        f = (z == -1.0 ? 0.25 : (z == 1 ? 0.75 : 0.50));
      }
    else if (z == -1.0) 
      { f = 0.0; }
    else if (z == 1.0) 
      { f = 1.0; }
    else if (d == 1)
      { /* The one-dimensional sphere is the interval [-1 _ +1]: */
        f = (z + 1)/2;
      }
    else
      { double w = asin(z);
        f = 0.5 + rn_ball_cap_vol_frac_ang(d,w);
      }
    /* fprintf(stderr, " = %.10f\n", f); */
    return f;
  }

double rn_ball_cap_vol_frac_ang(int d, double w)
  {
    /* fprintf(stderr, "    vol_frac_ang(%3d,%.10f)", d, w); */
    demand(d >= 0, "bad d");
    int p = (d % 2);
    double cw = cos(w);
    
    double K = (p == 0 ? cw : 1);
    double f = K;
    int i;
    for (i = 3 - p; i < d; i += 2)
      { K = K * cw*cw*((double)i-1)/((double) i); 
        f += K; 
      }
    f = f*sin(w);
    if (p == 0) { f = 2*(f + w)/M_PI; }
    f = 0.5*f;
    /* fprintf(stderr, " = %.10f\n", f); */
    return f;
  }
