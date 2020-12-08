/* See {neuromat_poly.h}. */
/* Last edited on 2014-03-26 16:11:37 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <rmxn.h>
#include <gauss_elim.h>

#include <neuromat_eeg.h>
#include <neuromat_poly.h>
   
void neuromat_poly_compute_lsq_matrix(int n, double x[], double w[], int g, double A[]);
  /* Assumes that {A} is a {g+1} by {g+1} matrix inearzed by rows.
    Fills it with the linear system's matrix for the weighted least squares 
    polynomial approximation problem at arguments {x[0..n-1]}
    with weight {w[0..n-1]}. 
    
    Namely, sets {A[i*(g+1)+j]} to {SUM{ w[k]*x[k]^i*x[k]^j : k \in 0..n-1 }}
    for all {i,j} in {0..g}. 
    
    If {x} is null, assumes {x[k] = 2*(k + 1/2)/n - 1} for {k} in {0..n-1}.
    If {w} is null, assumes unit weights for all arguments.*/

void neuromat_poly_compute_lsq_vector(int n, double x[], double y[], double w[], int g, double b[]);
  /* Assumes that {b} is a vector with {g+1} elements.
    Fills it with the linear system's right-hand-side vector for the weighted least squares 
    polynomial approximation problem at arguments {x[0..n-1]} of values {y[0..n-1]}
    with weight {w[0..n-1]}. 
    
    Namely, sets {b[i]} to {SUM{ w[k]*x[k]^i*y[k] : k \in 0..n-1 }}
    for all {i} in {0..g}. 
    
    If {x} is null, assumes {x[k] = 2*(k + 1/2)/n - 1} for {k} in {0..n-1}.
    If {w} is null, assumes unit weights for all arguments.*/

double neuromat_poly_bayes(double v, double pri_gud, double avg_gud, double dev_gud, double avg_bad, double dev_bad);
  /* Uses Bayes's formula to compute the probability of {v} being a good sample,
    assuming good and bad samples have normal distributions with averages
    {avg_gud,avg_bad} and deviations {dev_gud,dev_bad}, respectively, 
    and that, a priori, the sample is good with probability {pri_gud}. */

void neuromat_poly_compute_stats
  ( int n, 
    double v[], 
    double w[], 
    double p[], 
    bool_t good, 
    double *avgP, 
    double *devP, 
    double *priP
  );
  /* Estimates the average {*avgP} and deviation {*devP} of samples {v[0..n-1]}.
    If {good} is TRUE, estimates the parameters of the `good' samples (inliers),
    if {good} is FALSE estimates those of the `bad' ones (outliers).
    Also returns in {*priP} (if not null) the estimated `a priori' probability of a 
    sample being good or bad, respectively.

    Assumes that {w[k]} is the relevance weight of sample {v[k]},
    and {p[k]} is the current probability of that sample being a good sample. 
    If {w} is null, assumes {w[k]==1.0} for all {k}.
    If {p} is null, assumes {p[k]==0.5} for all {k}. */


double neuromat_poly_eval(int g, double P[], double x)
  {
    demand(g >= 0, "invalid power");
    double y = P[g];
    int i;
    for (i = g-1; i >= 0; i--) { y = P[i] + x*y; }
    return y;
  }

void neuromat_poly_eval_multi(int g, double P[], int n, double x[], double s[])
  {
    demand(g >= 0, "invalid power");
    int k;
    for (k = 0; k < n; k++) 
      { double xk = (x != NULL ? x[k] : 2*(k + 0.5)/n - 1);
        s[k] = neuromat_poly_eval(g, P, xk);
      }
  }

void neuromat_poly_compute_lsq_matrix(int n, double x[], double w[], int g, double A[])
  {
    demand(g >= 0, "invalid power");
    int g1 = g + 1;
    rmxn_zero(g1, g1, A);

    /* Fill the lower triangular half of {A}: */
    double *p = notnull(malloc(g1*sizeof(double)), "no mem"); /* Powers of each {x}. */
    int i, j, k;
    for (k = 0; k < n; k++)
      { double xk = (x != NULL ? x[k] : 2*(k + 0.5)/n - 1);
        double wk = (w != NULL ? w[k] : 1);
        for (i = 0; i <= g; i++)
          { p[i] = (i == 0 ? 1.0 : p[i-1]*xk);
            for (j = 0; j <= i; j++)
              { double pipj = p[i]*p[j];
                A[i*g1 + j] += wk*pipj;
              }
          }
      }
    free(p);

    /* Replicate the lower half of {A} into the upper half: */
    for (i = 1; i <= g; i++) { for (j = 0; j < i; j++) { A[j*g1 + i] = A[i*g1 + j]; } }
  }

void neuromat_poly_compute_lsq_vector(int n, double x[], double y[], double w[], int g, double b[])
  {
    demand(g >= 0, "invalid power");
    int g1 = g + 1;
    rn_zero(g1, b);
    int i, k;
    for (k = 0; k < n; k++)
      { double xk = (x != NULL ? x[k] : 2*(k + 0.5)/n - 1);
        double yk = y[k];
        double wk = (w != NULL ? w[k] : 1);
        double pi = 1.0; /* Power {xk^i}. */
        for (i = 0; i <= g; i++)
          { b[i] += wk*pi*yk;
            pi = pi*xk;
          }
      }
  }

void neuromat_poly_shift(int g, double P[], double a, double Q[])
  {
    demand(g >= 0, "invalid power");
    int k;
    for (k = g; k >= 0; k--)
      { /* Assume that {Q'=Q[k+1..g]} is the polynomial {P'=P[k+1..g]} shifted by {a},
          that is {Q'(x) = P'(x-a)} for all {x}. Now set {Q[k..g]} to the 
          coeffs of {Q''} that is {P''=P[k..g]} shifted by {a}.
          That is,  { Q''(x) = P''(x-a) = P[k]+(x-a)*P'(x-a) = P[k]+(x-a)*Q'(x) }
        */
        Q[k] = P[k];
        int i;
        for (i = k; i < g; i++) { Q[i] = Q[i] - a*Q[i+1]; }
      }
  }

void neuromat_poly_stretch(int g, double P[], double h, double Q[])
  {
    demand(g >= 0, "invalid power");
    double hk = 1.0;
    int k;
    for (k = 1; k <= g; k++)
      { hk *= h;  /* Value is {h^k}. */
        Q[k] = P[k]/hk;
      }
  }

void neuromat_poly_bezier(int g, int i, double a, double b, double P[])
  {
    demand((0 <= i) && (i <= g), "invalid Bezier exponent or index");
    double h = b - a;
    double ah = a/h, bh = b/h;
    P[0] = (double)comb(g,i); /* Hopefully there is no overflow or rounding. */
    /* Multiply {P[0..0]} by {((x-a)/h)^i} yielding {P[0..i]}: */
    int k, r;
    for (k = 0; k < i; k++)
      { /* Multiply {P[0..k]} by {(x-a)/h} yielding {P[0..k+1]}: */
        P[k+1] = P[k]/h;
        for (r = k; r > 0; r--) { P[r] = P[r-1]/h - P[r]*ah; }
        P[0] = -P[0]*ah;
      }
    /* Multiply {P[0..i]} by {((b-x)/h)^(g-i)} yielding {P[0..g]}: */
    for (k = i; k < g; k++)
      { /* Multiply {P[0..k]} by {(b-x)/h} yielding {P[0..k+1]}: */
        P[k+1] = -P[k]/h;
        for (r = k; r > 0; r--) { P[r] = -P[r-1]/h + P[r]*bh; }
        P[0] = P[0]*bh;
      }
  }

void neuromat_poly_fit(int n, double x[], double y[], double w[], int g, double P[])
  {
    demand(g >= 0, "invalid power");
    demand(n > g, "too few data points");
    int g1 = g + 1; /* Number of unnowns (coefficients). */
    double *A = notnull(malloc(g1*g1*sizeof(double)), "no mem"); /* Moment matrix. */
    double *b = notnull(malloc(g1*sizeof(double)), "no mem"); /* Right-hand side. */
    
    neuromat_poly_compute_lsq_matrix(n, x, w, g, A);
    rmxn_inv_full(g1, A, A);
    neuromat_poly_compute_lsq_vector(n, x, y, w, g, b);
    rmxn_map_col(g1, g1, A, b, P);
    free(A);
    free(b);
  }

void neuromat_poly_fit_robust
  ( int n, 
    double x[], 
    double y[], 
    double w[], 
    int maxiter, 
    int g, 
    double P[],
    neuromat_poly_report_proc_t *report
  )
  {
    /* bool_t debug = (report != NULL); */
    bool_t debug = TRUE;
    
    demand(g >= 0, "invalid power");
    demand(n > g, "too few data points");
    
    int g1 = g+1; /* Number of unnowns (coefficients). */
    double *A = notnull(malloc(g1*g1*sizeof(double)), "no mem"); /* Moment matrix. */
    double *b = notnull(malloc(g1*sizeof(double)), "no mem"); /* Right-hand side. */

    /* The matrix does not change: */
    neuromat_poly_compute_lsq_matrix(n, x, w, g, A);
    rmxn_inv_full(g1, A, A);

    /* Solve first without correction: */
    neuromat_poly_compute_lsq_vector(n, x, y, w, g, b);
    rmxn_map_col(g1, g1, A, b, P);
    if (report != NULL) { report(0, P, y); }

    if (maxiter > 0)
      { 
        /* Outlier suppression method: */
        double *ya = notnull(malloc(n*sizeof(double)), "no mem"); /* Current approximation. */
        double *yc = notnull(malloc(n*sizeof(double)), "no mem"); /* Deviations or fake values. */

        double *p = notnull(malloc(n*sizeof(double)), "no mem"); /* Probability of each datum to be inlier. */
        
        /* Assumed parameters of inlier residuals and outlier values: */
        double avg_gud, avg_bad;   /* Average. */
        double dev_gud, dev_bad;   /* Deviation. */
        double pri_gud, pri_bad;   /* A priori probability of being inlier or outlier. */
        
        int iter;
        rn_all(n, 0.5, p); /* A priori, inliers and outliers are equally likely: */
        for (iter = 1; iter <= maxiter; iter++)
          { 
            if (debug) { fprintf(stderr, "    iteration %d\n", iter); }

            /* Set {ya[0..n-1]} to current approximation: */
            neuromat_poly_eval_multi(g, P, n, x, ya);
            
            /* Recompute parameters of inlier deviations: */
            rn_sub(n, y, ya, yc); /* Temporary: residuals of current spproximation. */
            neuromat_poly_compute_stats(n, yc, w, p, TRUE,  &avg_gud, &dev_gud, &pri_gud);
            if (debug) { fprintf(stderr, "      gud:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_gud, dev_gud, pri_gud); }
            
            /* Recompute parameters of outlier values: */
            neuromat_poly_compute_stats(n, y, w, p, FALSE, &avg_bad, &dev_bad, &pri_bad);
            if (debug) { fprintf(stderr, "      bad:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_bad, dev_bad, pri_bad); }

            assert(fabs(pri_bad + pri_gud - 1.0) < 0.0001);
            
            /* Widen the inlier deviation distribution to account for the shrink effect: */
            double beta = 1.1; /* !!! Figure out the right factor !!! */
            dev_gud = beta*dev_gud + 1.0e-10*dev_bad + 1.0e-200;
            if (debug) { fprintf(stderr, "      gud:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_gud, dev_gud, pri_gud); }

            /* Make sure that the outlier distrib is broader than the inlier one: */
            double alpha = 2.0; /* Difference factor. */
            if (dev_bad < dev_gud*alpha) { dev_bad = dev_gud*alpha; }
            if (debug) { fprintf(stderr, "      bad:  avg = %+9.4f  dev = %9.4f  pri = %9.4f\n", avg_bad, dev_bad, pri_bad); }

            /* Recompute data point inlier/outlier probabilities and adjusted data {yc[0..n-1]}: */
            int k;
            for (k = 0; k < n; k++) 
              { /* Decide the probability {p[k]} of each sample {y[k]} being an inlier: */
                double pk = neuromat_poly_bayes(y[k], pri_gud, ya[k] + avg_gud, dev_gud, avg_bad, dev_bad);
                p[k] = pk;
                /* Set {yc[k]} to the expected value of the inlier part of {y[k]}: */
                yc[k] = pk*y[k] + (1 - pk)*ya[k];
              }
              
            /* Recompute coeffs using {yc[0..n-1]} instead of {y[0..n-1]}: */
            neuromat_poly_compute_lsq_vector(n, x, yc, w, g, b);
            rmxn_map_col(g1, g1, A, b, P);

            if (report != NULL) { report(iter, P, yc); }
          }
        free(ya);
        free(yc);
        free(p);
      }
  }

void neuromat_poly_compute_stats
  ( int n, 
    double v[], 
    double w[], 
    double p[], 
    bool_t good, 
    double *avgP, 
    double *devP, 
    double *priP
  )
  {
    /* Compute the average and overall probability: */
    double sum_wpv = 0;
    double sum_wp = 0;
    double sum_w = 0;
    int k;
    for (k = 0; k < n; k++) 
      { double vk = v[k];
        double wk = (w == NULL ? 1.0 : w[k]);
        double pk = (p == NULL ? 0.5 : (good ? p[k] : 1.0 - p[k]));
        sum_wpv += wk*pk*vk;
        sum_wp += wk*pk;
        sum_w += wk;
      }
    double avg = sum_wpv/(sum_wp + 1.0e-300);
    double pri = sum_wp/(sum_w + 1.0e-300);
    assert(! isnan(avg));
    /* Compute the standard deviation: */
    double sum_wpd2 = 0;
    for (k = 0; k < n; k++) 
      { double dk = v[k] - avg;
        double wk = (w == NULL ? 1.0 : w[k]);
        double pk = (p == NULL ? 0.5 : (good ? p[k] : 1.0 - p[k]));
        sum_wpd2 += wk*pk*dk*dk;
      }
    double var = sum_wpd2/(sum_wp + 1.0e-300);
    double dev = sqrt(var); /* Should compensate for bias... */
    /* Return: */
    assert(! isnan(dev));
    assert(! isnan(pri));
    (*avgP) = avg;
    (*devP) = dev;
    if (priP != NULL) { (*priP) = pri; }
  }

double neuromat_poly_bayes(double v, double pri_gud, double avg_gud, double dev_gud, double avg_bad, double dev_bad)
  {
    double dk_gud = (v - avg_gud)/dev_gud;
    double dk_bad = (v - avg_bad)/dev_bad;
    if (dk_gud > 6.0)
      { return 0.0; }
    else if (dk_bad > 6.0)
      { return 1.0; }
    else
      { double pri_bad = 1 - pri_gud;
        double P_gud = exp(-dk_gud*dk_gud/2)/dev_gud; /* Prob of sample times {sqrt(2*PI)}, assuming good. */
        double P_bad = exp(-dk_bad*dk_bad/2)/dev_bad; /* Prob of sample times {sqrt(2*PI)}, assuming bad. */ 
        double PP_gud = pri_gud*P_gud;
        double PP_bad = pri_bad*P_bad;
        return PP_gud/(PP_gud + PP_bad);
      }
  }
