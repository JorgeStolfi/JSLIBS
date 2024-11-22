/* See dnae_weight.h */
/* Last edited on 2024-11-20 06:00:11 by stolfi */

#define dnae_weight_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>

#include <vec.h>

#include <msm_basic.h>

#include <dnae_weight.h>

/* !!! Check if we can use the {weight_table.h} in JSLIBS/libjs. !!! */

/* INTERNAL PROTOTYPES */

bool_t check_normalization(double_vec_t *wtb, double tol, bool_t die);
  /* Checks whether {wtb.e[0..wtb.ne-1]} add to 1. If not, either
    dies with error (if {die} is TRUE) or returns FALSE silently 
    (if {die} is FALSE). */

double dnae_weight_table_gaussian_loss(int r, double sigma, bool_t odd);
  /* Returns the weight that is lost if a table with {N} entries is
    used to represent a Gaussian weight distribution with standard
    deviation {sigma}; where {N} is {2*r} if {odd} is FALSE, and
    {2*r+1} if {odd} is TRUE. */

double dnae_weight_table_gaussian_entry(int k, double sigma, bool_t odd);
  /* Returns entry {r+k} of a table with {N} entries that
    represents a Gaussian weight distribution with standard
    deviation {sigma}; where {N} is {2*r} if {odd} is FALSE, and
    {2*r+1} if {odd} is TRUE. */

/* IMPLEMENTATIONS */

void dnae_weight_table_print(FILE *wr, char *wtname, double_vec_t *wt)
  { int n = wt->ne;
    double ctr = ((double)n-1)/2;
    double radius = ctr;
    fprintf(wr, "weight table\n");
    fprintf(wr, "name = %s\n", wtname);
    double sum_w = 0;   /* Sum of {wt[i]}. */
    double sum_iw = 0;  /* Sum of {i*wt[i]}. */
    double sum_bw[2] = { 0, 0 }; /* Sum of even and odd elements. */
    int i;
    for (i = 0; i < n; i++)
      { double w = wt->e[i];
        fprintf(wr, "  w[%3d] (%+6.1f) = %12.8f\n", i, i - ctr, w);
        sum_w += w;
        sum_iw += i*w;
        sum_bw[i%2] += w;
      }
    fprintf(wr, "\n");
    fprintf(wr, "width =     %4d\n", n);
    fprintf(wr, "radius =    %6.1f\n", radius);
    fprintf(wr, "sum =       %13.8f\n", sum_w);
    if (sum_w == 0) { sum_w = 1.0e-200; }
    fprintf(wr, "unbalance = %13.8f\n", (sum_bw[1] - sum_bw[0])/sum_w);
    double avg = sum_iw/sum_w;
    fprintf(wr, "mean =      %13.8f\n", avg);
    double sum_d2w = 0;
    for (i = 0; i < n; i++)
      { double w = wt->e[i];
        double d = i - avg;
        sum_d2w += d*d*w; 
      }
    double var = sum_d2w/sum_w;
    fprintf(wr, "variance =  %13.8f\n", var);
    fprintf(wr, "deviation = %13.8f\n", sqrt(var));
    fprintf(wr, "\n");
  }
double_vec_t dnae_weight_table_make_gaussian(double sigma, double maxLoss)
  { /* Paranoia - check functions: */
    double win = dnae_weight_table_gaussian_entry(0, sigma, TRUE);
    double wot = dnae_weight_table_gaussian_loss(0, sigma, TRUE);
    assert(fabs(1 - (win+wot)) < 1.0e-10);
    /* Find {r} so that the omitted weight is at most {maxLoss}: */
    int r = 0;
    while (dnae_weight_table_gaussian_loss(r, sigma, TRUE) > maxLoss) { r++; }
    assert(r > 0);
    /* Allocate table: */
    int N = 2*r + 1;
    double_vec_t wtb = double_vec_new(N);
    /* Compute entries {wtb[i]} and their sum {sumw}: */
    double sumw = 0;
    int k;
    for (k = 0; k <= r; k++) 
      { double w = dnae_weight_table_gaussian_entry(k, sigma, TRUE);
        wtb.e[r+k] = wtb.e[r-k] = w;
        sumw += (k == 0 ? w : 2*w);
      }
    /* Normalize table to unit sum: */
    int i;
    for (i = 0; i < N; i++) { wtb.e[i] /= sumw; }
    /* Return {wtb}: */
    assert(wtb.ne == N);
    return wtb;
  }
  
double dnae_weight_table_gaussian_loss(int r, double sigma, bool_t odd)
  { double r1 = ((double)r) + (odd ? 0.5 : 0.0); 
    double t1 = r1/sigma; 
    double ws = erfc(t1/M_SQRT2); 
    fprintf(stderr, "  sigma = %12.8f r = %2d t = %12.8f loss = %12.8f\n", sigma, r, t1, ws); 
    return ws;
  }

double dnae_weight_table_gaussian_entry(int k, double sigma, bool_t odd)
  { double r0 = ((double)k) - (odd ? 0.5 : 0.0), t0 = r0/sigma; 
    double r1 = r0 + 1, t1 = r1/sigma; 
    double wr = (erf(t1/M_SQRT2) - erf(t0/M_SQRT2))/2;
    return wr;
  }

double_vec_t dnae_weight_table_make_binomial(int r)
  { /* Compute table size {N} and allocate table: */
    int N = 2*r + 1;
    double_vec_t wtb = double_vec_new(N);
    /* Compute distribution by repeated convolutions (Pascal's triangle): */
    wtb.e[0] = 1;
    int k;
    for (k = 1; k < N; k++)
      { int j;
        wtb.e[k] = wtb.e[k-1]/2;
        for (j = k-1; j > 0; j--)
          { wtb.e[j] = wtb.e[j]/2 + wtb.e[j-1]/2; }
        wtb.e[0] /= 2;
      }
    /* Unit sum property should hold exactly for {N<50}: */
    check_normalization(&wtb, 0.0, TRUE);
    /* Return {wtb}: */
    assert(wtb.ne == N);
    return wtb;
  }

double_vec_t dnae_weight_table_make_triangular(int r)
  { /* Compute table size {N} and allocate table: */
    int N = 2*r + 1;
    double_vec_t wtb = double_vec_new(N);
    /* Compute distribution: */
    double W = (r+1)*(r+1);
    int k;
    for (k = 0; k <= r; k++)
      { wtb.e[k] = wtb.e[2*r-k] = ((double)k+1)/W; }
    /* Unit sum property should hold except for roundoff: */
    check_normalization(&wtb, 1.0e-10, TRUE);
    /* Return {wtb}: */
    assert(wtb.ne == N);
    return wtb;
  }
   
bool_t check_normalization(double_vec_t *wtb, double tol,bool_t die)
  { /* Check unit sum property (should hold exactly for {N<50} or so): */
    double sumw = 0;
    int i;
    for (i = 0; i < wtb->ne; i++) { sumw += wtb->e[i]; }
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
   
void dnae_weight_table_pair_make_gaussian
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  )   
  { double maxLoss = 1.0e-3;
    double sigma0 = sqrt(var0);
    (*wtb0) = dnae_weight_table_make_gaussian(sigma0, maxLoss);
    (*wname0) = jsprintf("gaussian(sigma = %8.6f)", sigma0);
    if (verbose) { dnae_weight_table_print(stderr, (*wname0), wtb0); }
    
    double sigma1 = sqrt(var1);
    (*wtb1) = dnae_weight_table_make_gaussian(sigma1, maxLoss);
    (*wname1) = jsprintf("gaussian(sigma = %8.6f)", sigma1);
    if (verbose) { dnae_weight_table_print(stderr, (*wname1), wtb1); }
  }
   
void dnae_weight_table_pair_make_binomial
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  )   
  { /* Variance of the distribution {w[k] = choose(n,k)} is {n/4 = r/2}, so: */
    int r0 = (int)msm_round(2*var0);
    (*wtb0) = dnae_weight_table_make_binomial(r0);
    (*wname0) = jsprintf("binomial(n = %d, sigma = %8.6f)", 2*r0, sqrt(0.5*r0));
    if (verbose) { dnae_weight_table_print(stderr, (*wname0), wtb0); }
    
    int r1 = (int)msm_round(2*var1);
    (*wtb1) = dnae_weight_table_make_binomial(r1);
    (*wname1) = jsprintf("binomial(n = %d, sigma = %8.6f)", 2*r1, sqrt(0.5*r1));
    if (verbose) { dnae_weight_table_print(stderr, (*wname1), wtb1); }
  }
   
void dnae_weight_table_pair_make_triangular
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  )   
  { /* Variance of the distribution {w[k] = r + 1 - |r-k|} is {~r^2/4}, so: */
    int r0 = (int)msm_round(sqrt(4*var0));
    (*wtb0) = dnae_weight_table_make_triangular(r0);
    (*wname0) = jsprintf("triangular(r = %d)", r0);
    if (verbose) { dnae_weight_table_print(stderr, (*wname0), wtb0); }
    
    int r1 = (int)msm_round(sqrt(4*var1));
    (*wtb1) = dnae_weight_table_make_triangular(r1);
    (*wname1) = jsprintf("triangular(r = %d)", r1);
    if (verbose) { dnae_weight_table_print(stderr, (*wname1), wtb1); }
  }

double_vec_t dnae_weight_args_parse(argparser_t *pp, bool_t unitNorm)
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

char *dnae_weight_make_descr(double_vec_t *wt, char *fmt)
  { int nw = wt->ne;
    char_vec_t d = char_vec_new(nw*8); /* Buffer for description. */
    int nc = 0; /* Number of characters in description. */
    /* Start with open bracket: */
    char_vec_expand(&d, nc); d.e[nc] = '['; nc++;
    /* Append the elements: */
    int i;
    for (i = 0; i < nw; i++) 
      { /* Append a space to {d}: */
        char_vec_expand(&d, nc); d.e[nc] = ' '; nc++;
        /* Convert the element {wt[i]} to string: */
        char *wi = jsprintf(fmt, wt->e[i]);
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
