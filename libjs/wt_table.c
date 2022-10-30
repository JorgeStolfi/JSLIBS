/* See wt_table.h */
/* Last edited on 2022-10-30 19:28:47 by stolfi */

#define wt_table_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <gauss_distr.h>

#include <wt_table.h>

/* IMPLEMENTATIONS */

double_vec_t wt_table_make_gaussian(double sigma, double maxLoss, bool_t norm)
  { /* Find {r} so that the omitted weight is at most {maxLoss}: */
    int32_t r = 0;
    while (wt_table_gaussian_loss(2*r + 1, sigma) > maxLoss) { r++; }
    assert(r > 0);
    /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(2*r + 1);
    wt_table_fill_gaussian(sigma, wt.ne, wt.e, norm);
    return wt;
  }
  
double wt_table_gaussian_loss(int32_t n, double sigma)
  { double r1 = 0.5*n, t1 = r1/sigma; 
    double ws = erfc(t1/M_SQRT2); 
    return ws;
  }

double_vec_t wt_table_make_binomial(int32_t r, bool_t norm)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(2*r + 1);
    wt_table_fill_binomial(wt.ne, wt.e, norm);
    return wt;
  }

double_vec_t wt_table_make_triangular(int32_t r, bool_t norm)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(2*r + 1);
    wt_table_fill_triangular(wt.ne, wt.e, norm);
    return wt;
  }
   
double_vec_t wt_table_make_hann(int32_t r, bool_t norm)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(2*r + 1);
    wt_table_fill_hann(wt.ne, wt.e, norm);
    return wt;
  }
   
double wt_table_avg(int32_t n, double wt[])
  { double sum_w = 0;   /* Sum of {wt[k]}. */
    double sum_iw = 0;  /* Sum of {k*wt[k]}. */
    for (int32_t k = 0; k < n; k++)
      { double w = wt[k];
        sum_w += w;
        sum_iw += k*w;
      }
    double avg = sum_iw/sum_w;
    return avg;
  }
   
double wt_table_var(int32_t n, double wt[], double avg)
  { double sum_w = 0;   /* Sum of {wt[k]}. */
    double sum_d2w = 0; /* Sum of {(k-avg)^2*wt[k]}. */
    for (int32_t k = 0; k < n; k++)
      { double w = wt[k];
        double d = k - avg;
        sum_w += w;
        sum_d2w += d*d*w; 
      }
    double var = sum_d2w/sum_w;
    return var;
  }

void wt_table_normalize_sum(int32_t n, double wt[])
  { double sum = 0.0;
    for (int32_t k = 0; k < n; k++) { sum += wt[k]; }
    if (sum != 0.0)
      { for (int32_t k = 0; k < n; k++) { wt[k] /= sum; } }
  }

void wt_table_shifted_sum(int32_t n, double wt[], int32_t stride, double ws[])
  {
    demand(stride >= 1, "invalid {stride}");
    for (int32_t ka = 0; ka < n; ka++)
      { int32_t imin = -(ka/stride);
        int32_t imax = (n-1-ka)/stride;
        double sum = 0.0; /* Sum of all weights in a {stride} train. */
        for (int32_t i = imin; i <= imax; i++)
          { int32_t kb = ka + i*stride;
            assert((0 <=kb) && (kb < n));
            sum += wt[kb];
          }
        ws[ka] = sum;
      }
  }

void wt_table_print(FILE *wr, char *wtname, int32_t n, double wt[], int32_t stride)
  { double ctr = ((double)n-1)/2;
    double radius = ctr;
    fprintf(wr, "weight table\n");
    if (wtname != NULL) { fprintf(wr, "name = %s\n", wtname); }
    fprintf(wr, "length = %d\n", n);
    double ws[n];
    if (stride != 0) 
      { fprintf(wr, "stride = %d\n", stride); 
        wt_table_shifted_sum(n, wt, stride, ws);
      }
    double sum_w = 0;   /* Sum of {wt[k]}. */
    double sum_bw[2] = { 0, 0 }; /* Sum of even and odd elements. */
    double wsExp = 1.0/stride; /* Expected value of overlapped windows. */
    for (int32_t k = 0; k < n; k++)
      { double wa = wt[k];
        fprintf(wr, "  w[%03d] (%+6.1f) = %10.8f", k, k - ctr, wa);
        sum_w += wa;
        sum_bw[k%2] += wa;
        if (stride != 0)
          { fprintf(wr, " shifted sum = %18.16f", ws[k]);
            double err = ws[k] - wsExp;
            fprintf(wr, " err = %24.16e", err);
          }
        fprintf(wr, "\n");
      }

    fprintf(wr, "\n");
    fprintf(wr, "width =     %4d\n", n);
    fprintf(wr, "radius =    %6.1f\n", radius);
    fprintf(wr, "sum =       %13.8f\n", sum_w);
    if (sum_w == 0) { sum_w = 1.0e-200; }
    fprintf(wr, "unbalance = %13.8f\n", (sum_bw[1] - sum_bw[0])/sum_w);
    double avg = wt_table_avg(n, wt);
    fprintf(wr, "mean =      %13.8f\n", avg);
    double sum_d2w = 0;
    for (int32_t k = 0; k < n; k++)
      { double w = wt[k];
        double d = k - avg;
        sum_d2w += d*d*w; 
      }
    double var = wt_table_var(n, wt, avg);
    fprintf(wr, "variance =  %13.8f\n", var);
    fprintf(wr, "deviation = %13.8f\n", sqrt(var));
    fprintf(wr, "\n");
  }
  
void wt_table_fill_gaussian(double sigma, int32_t n, double wt[], bool_t norm)
  { /* Compute entries {wt[k]} and their sum {sumw}: */
    for (int32_t k = 0, j = n-1; k < n; k++, j--) 
      { wt[k] = (k <= j ? wt_table_gaussian_entry(n, k, sigma) : wt[j]); }
    if (norm) { wt_table_normalize_sum(n, wt); }
  }
  
double wt_table_gaussian_entry(int32_t n, int32_t k, double sigma)
  { /* Integral of Gaussian within interval {k} of {n} intervals centered at origin: */
    double r0 = ((double)k) - 0.5*n;
    double r1 = r0 + 1.0; 
    double wr = gauss_distr_integral(r0, r1, 0.0, sigma);
    return wr;
  }

void wt_table_fill_binomial(int32_t n, double wt[], bool_t norm)
  { /* Repeated convolutions (Pascal's triangle): */
    wt[0] = 1;
    for (int32_t k = 1; k < n; k++)
      { wt[k] = wt[k-1];
        for (int32_t j = k-1; j > 0; j--) { wt[j] = wt[j] + wt[j-1]; }
        if (norm) { for (int32_t j = 0; j <= k; j++) { wt[j] /= 2; } }
      }
  }

void wt_table_fill_triangular(int32_t n, double wt[], bool_t norm)
  { /* Build triangular table and compute its sum: */
    double c = (n-1)/2.0;
    double h = (n+1)/2.0;
    for (int32_t k = 0; k < n; k++) { wt[k] = 1.0 - fabs((k-c)/h); }
    if (norm) { wt_table_normalize_sum(n, wt); }
  }
   
void wt_table_fill_hann(int32_t n, double wt[], bool_t norm)
  { /* Compute parameters: */
    double c = ((double)n-1)/2.0;  /* Center of {[0 _ n]} (may be half-integer). */
    double h = ((double)n+1)/2.0;  /* To make {wt[-1]=wt[n]=0}. */
    for (int32_t k = 0; k < n; k++)
      { wt[k] = (k <= n-1-k ? (1 + cos(M_PI*(k-c)/h))/2 : wt[n-1-k]); }
    if (norm) { wt_table_normalize_sum(n, wt); }
  }

bool_t wt_table_check_normalization(int32_t n, double wt[], double tol,bool_t die)
  { /* Check unit sum property: */
    double sumw = 0;
    for (int32_t k = 0; k < n; k++) { sumw += wt[k]; }
    double err = sumw - 1;
    if (isnan(sumw) || (fabs(err) > tol))
      { if (die)
          { fprintf(stderr, "normalization failure: sumw = %+24.15e  err = %+24.15e\n", sumw, err);
            assert(FALSE);
          }
        else
          { return FALSE; }
      }
    else
      { return TRUE; }
  }
  
bool_t wt_table_check_partition_of_unity(int32_t n, double wt[], int32_t stride, double tol, bool_t die)
  {
    demand(stride >= 1, "invalid {stride}");
    double ws[n];
    wt_table_shifted_sum(n, wt, stride, ws);
    for (int32_t k = 0; k < n; k++) 
      { if (fabs(ws[k] - 1.0) > tol)
          { if (die)
              { fprintf(stderr, "table is not  partition of unity");
                fprintf(stderr, " {ws[%d] = %18.16f", k, ws[k]);
                double err = ws[k] - 1.0;
                fprintf(stderr, " err = %24.16e\n", err);
                assert(FALSE);
              }
            else
              { return FALSE; }
          }
      }
    return TRUE;
  }

double_vec_t wt_table_args_parse(argparser_t *pp, bool_t unitNorm)
  { /* Allocate the weight vector, with size unknown: */
    double_vec_t w = double_vec_new(50);
    /* Parse the follwoing numeric arguments, save them in {w[0..nw-1]}: */
    int32_t nw = 0;
    while (argparser_next_is_number(pp))
      { double wt = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        double_vec_expand(&w, nw);
        w.e[nw] = wt; nw++;
      }
    /* Parse the optional "/ {DENOM}" args: */
    double denom = 1.0;
    if (argparser_keyword_present_next(pp, "/"))
      { double den = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        if (! unitNorm) { denom = den; }
      }
    if (unitNorm)
      { /* Set [denom} to the sum of all weights: */
        denom = 0;
        for (int32_t k = 0; k < nw; k++) { denom += w.e[k]; }
      }
    if (denom != 1)
      { /* Divide all weights by {denom}: */
        for (int32_t k = 0; k < nw; k++) { w.e[k] /= denom; }
      }
    double_vec_trim(&w, nw);
    return w;
  }

char *wt_table_make_descr(int32_t n, double wt[], char *fmt)
  { char_vec_t d = char_vec_new(n*8); /* Buffer for description. */
    int32_t nc = 0; /* Number of characters in description. */
    /* Start with open bracket: */
    char_vec_expand(&d, nc); d.e[nc] = '['; nc++;
    /* Append the elements: */
    for (int32_t k = 0; k < n; k++) 
      { /* Append a space to {d}: */
        char_vec_expand(&d, nc); d.e[nc] = ' '; nc++;
        /* Convert the element {wt[k]} to string: */
        char *wk = NULL;
        asprintf(&wk, fmt, wt[k]);
        /* Append it to {d}: */
        char *p = wk;
        while ((*p) != 0) { char_vec_expand(&d, nc); d.e[nc] = (*p); p++; nc++; }
        /* Cleanup: */
        free(wk);
      }
    /* Append a space, a close-bracket, and zero char: */
    char_vec_expand(&d, nc); d.e[nc] = ' '; nc++;
    char_vec_expand(&d, nc); d.e[nc] = ']'; nc++;
    char_vec_expand(&d, nc); d.e[nc] = 0; nc++;
    
    /* Trim and return: */
    char_vec_trim(&d, nc);
    return d.e;
  }
