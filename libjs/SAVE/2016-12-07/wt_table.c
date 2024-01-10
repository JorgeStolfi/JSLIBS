/* See wt_table.h */
/* Last edited on 2013-05-04 14:32:40 by stolfilocal */

#define wt_table_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>

#include <wt_table.h>

/* INTERNAL PROTOTYPES */

bool_t wt_table_check_normalization(int n, double wt[], double tol, bool_t die);
  /* Checks whether {wtb.e[0..wtb.ne-1]} add to 1. If not, either
    dies with error (if {die} is TRUE) or returns FALSE silently 
    (if {die} is FALSE). */

double wt_table_gaussian_loss(int n, double sigma);
  /* Returns the weight that is lost if a table with {n} entries is
    used to represent a Gaussian weight distribution with standard
    deviation {sigma}. */

double wt_table_gaussian_entry(int k, double sigma, bool_t odd);
  /* Computes the entry with position {k} in a table that represents a
    Gaussian weight distribution with standard deviation {sigma}. The
    position {k} starts from 0, and is counted from the central entry,
    if the table has odd length (indicated by {odd=TRUE}); or from the
    entry just after the center, if the table has even length. 
    
    !!! Inefficient -- should reuse the shared {erf} calls !!! */

/* IMPLEMENTATIONS */

double_vec_t wt_table_make_gaussian(double sigma, double maxLoss)
  { /* Find {r} so that the omitted weight is at most {maxLoss}: */
    int r = 0;
    while (wt_table_gaussian_loss(2*r + 1, sigma) > maxLoss) { r++; }
    assert(r > 0);
    /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(2*r + 1);
    wt_table_fill_gaussian(sigma, wt.ne, wt.e);
    return wt;
  }
  
double wt_table_gaussian_loss(int n, double sigma)
  { double r1 = 0.5*n, t1 = r1/sigma; 
    double ws = erfc(t1/M_SQRT2); 
    return ws;
  }

double_vec_t wt_table_make_binomial(int r)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(2*r + 1);
    wt_table_fill_binomial(wt.ne, wt.e);
    return wt;
  }

double_vec_t wt_table_make_triangular(int r)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(2*r + 1);
    wt_table_fill_triangular(wt.ne, wt.e);
    return wt;
  }
   
double_vec_t wt_table_make_hann(int r)
  { /* Allocate and fill the table: */
    double_vec_t wt = double_vec_new(2*r + 1);
    wt_table_fill_hann(wt.ne, wt.e);
    return wt;
  }
   
double wt_table_avg(int n, double wt[])
  { double sum_w = 0;   /* Sum of {wt[i]}. */
    double sum_iw = 0;  /* Sum of {i*wt[i]}. */
    int i;
    for (i = 0; i < n; i++)
      { double w = wt[i];
        sum_w += w;
        sum_iw += i*w;
      }
    double avg = sum_iw/sum_w;
    return avg;
  }
   
double wt_table_var(int n, double wt[], double avg)
  { double sum_w = 0;   /* Sum of {wt[i]}. */
    double sum_d2w = 0;
    int i;
    for (i = 0; i < n; i++)
      { double w = wt[i];
        double d = i - avg;
        sum_w += w;
        sum_d2w += d*d*w; 
      }
    double var = sum_d2w/sum_w;
    return var;
  }

void wt_table_print(FILE *wr, char *wtname, int n, double wt[])
  { double ctr = ((double)n-1)/2;
    double radius = ctr;
    fprintf(wr, "weight table\n");
    if (wtname != NULL) { fprintf(wr, "name = %s\n", wtname); }
    double sum_w = 0;   /* Sum of {wt[i]}. */
    double sum_bw[2] = { 0, 0 }; /* Sum of even and odd elements. */
    int i;
    for (i = 0; i < n; i++)
      { double w = wt[i];
        fprintf(wr, "  w[%3d] (%+6.1f) = %12.8f\n", i, i - ctr, w);
        sum_w += w;
        sum_bw[i%2] += w;
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
    for (i = 0; i < n; i++)
      { double w = wt[i];
        double d = i - avg;
        sum_d2w += d*d*w; 
      }
    double var = wt_table_var(n, wt, avg);
    fprintf(wr, "variance =  %13.8f\n", var);
    fprintf(wr, "deviation = %13.8f\n", sqrt(var));
    fprintf(wr, "\n");
  }
  
void wt_table_fill_gaussian(double sigma, int n, double wt[])
  { /* Paranoia - check functions: */
    double win = wt_table_gaussian_entry(1, sigma, TRUE);
    double wot = wt_table_gaussian_loss(1, sigma);
    assert(fabs(1 - (win+wot)) < 1.0e-10);
    /* Compute entries {wt[i]} and their sum {sumw}: */
    int r = (n-1)/2;
    bool_t odd = ((n % 2) != 0);
    double sumw = 0;
    int i;
    for (i = 0; i < n; i++) 
      { double w = (i > r ? wt[n-1-i] : wt_table_gaussian_entry(r-i, sigma, odd));
        wt[i] = w; sumw += w;
      }
    /* Normalize table to unit sum: */
    for (i = 0; i < n; i++) { wt[i] /= sumw; }
  }
  
double wt_table_gaussian_entry(int k, double sigma, bool_t odd)
  { /* Integral of Gaussian within each pixel: */
    double r0 = ((double)k) - (odd ? 0.5 : 0.0), t0 = r0/sigma; 
    double r1 = r0 + 1, t1 = r1/sigma; 
    double wr = (erf(t1/M_SQRT2) - erf(t0/M_SQRT2))/2;
    return wr;
  }

void wt_table_fill_binomial(int n, double wt[])
  { /* Repeated convolutions (Pascal's triangle): */
    wt[0] = 1;
    int k;
    for (k = 1; k < n; k++)
      { int j;
        wt[k] = wt[k-1]/2;
        for (j = k-1; j > 0; j--)
          { wt[j] = wt[j]/2 + wt[j-1]/2; }
        wt[0] /= 2;
      }
    if (n < 40)
      { /* The unit sum property should hold exactly for {n<50}: */
        wt_table_check_normalization(n, wt, 0.0, TRUE);
      }
  }

void wt_table_fill_triangular(int n, double wt[])
  { /* Compute distribution: */
    int r = (n-1)/2;
    demand(n == 2*r + 1, "n must be odd");
    double W = (r+1)*(r+1);
    int k;
    for (k = 0; k <= r; k++)
      { wt[k] = wt[2*r-k] = ((double)k+1)/W; }
    /* Unit sum property should hold except for roundoff: */
    wt_table_check_normalization(n, wt, 1.0e-10, TRUE);
  }
   
void wt_table_fill_hann(int n, double wt[])
  { /* Compute parameters: */
    double c = ((double)n-1)/2.0;
    double h = ((double)n)/2.0;
    double sumw = 0;
    int k;
    for (k = 0; k < n/2; k++)
      { double wk = (1 + cos(M_PI*(k-c)/h))/2; 
        wt[k] = wt[n-1-k] = wk;
        sumw += wk;
      }
    /* Normalize table to unit sum: */
    for (k = 0; k < n; k++) { wt[k] /= sumw; }
  }

bool_t wt_table_check_normalization(int n, double wt[], double tol,bool_t die)
  { /* Check unit sum property: */
    double sumw = 0;
    int i;
    for (i = 0; i < n; i++) { sumw += wt[i]; }
    double err = sumw - 1;
    if (fabs(err) > tol)
      { if (die)
          { fprintf(stderr, "normalization failure: sumw = 1 %+12.4e\n", err);
            assert(FALSE);
          }
        else
          { return FALSE; }
      }
    else
      { return TRUE; }
  }

double_vec_t wt_table_args_parse(argparser_t *pp, bool_t unitNorm)
  { /* Allocate the weight vector, with size unknown: */
    double_vec_t w = double_vec_new(50);
    /* Parse the follwoing numeric arguments, save them in {w[0..nw-1]}: */
    int nw = 0;
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
        int k;
        for (k = 0; k < nw; k++) { denom += w.e[k]; }
      }
    if (denom != 1)
      { /* Divide all weights by {denom}: */
        int k;
        for (k = 0; k < nw; k++) { w.e[k] /= denom; }
      }
    double_vec_trim(&w, nw);
    return w;
  }

char *wt_table_make_descr(int n, double wt[], char *fmt)
  { char_vec_t d = char_vec_new(n*8); /* Buffer for description. */
    int nc = 0; /* Number of characters in description. */
    /* Start with open bracket: */
    char_vec_expand(&d, nc); d.e[nc] = '['; nc++;
    /* Append the elements: */
    int i;
    for (i = 0; i < n; i++) 
      { /* Append a space to {d}: */
        char_vec_expand(&d, nc); d.e[nc] = ' '; nc++;
        /* Convert the element {wt[i]} to string: */
        char *wi = NULL;
        asprintf(&wi, fmt, wt[i]);
        /* Append it to {d}: */
        char *p = wi;
        while ((*p) != 0) { char_vec_expand(&d, nc); d.e[nc] = (*p); p++; nc++; }
        /* Cleanup: */
        free(wi);
      }
    /* Append a space, a close-bracket, and zero char: */
    char_vec_expand(&d, nc); d.e[nc] = ' '; nc++;
    char_vec_expand(&d, nc); d.e[nc] = ']'; nc++;
    char_vec_expand(&d, nc); d.e[nc] = 0; nc++;
    
    /* Trim and return: */
    char_vec_trim(&d, nc);
    return d.e;
  }
