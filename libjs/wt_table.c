/* See wt_table.h */
/* Last edited on 2024-11-06 01:30:28 by stolfi */

#define wt_table_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <jsstring.h>

#include <wt_table.h>
  
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

double_vec_t wt_table_convolution(int32_t n1, double wt1[], int32_t n2, double wt2[], int32_t stride)
  {
    demand(stride >= 1, "invalid {stride}");
    int32_t ns = n1 + (n2-1)*stride;
    double_vec_t ws = double_vec_new(ns);
    for (int32_t i = 0; i < ns; i++)
      { int32_t k2min = (i < n1 ? 0 : (i - n1 + stride)/stride);
        int32_t k2max = (i >= n2*stride ? n2 - 1 : i/stride);
        assert(k2min >= 0);
        assert(k2max < n2);
        double sum = 0;
        for (int32_t k2 = k2min; k2 <= k2max; k2++)
          { int32_t k1 = i - k2*stride;
            assert(k1 >= 0);
            assert(k1 < n1);
            sum += wt1[k1]*wt2[k2];
          }
        ws.e[i] = sum;
      }
    return ws;
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
        fprintf(wr, "  w[%03d] (%+6.1f) = %18.16", k, k - ctr, wa);
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
    fprintf(wr, "sum =       %21.16f\n", sum_w);
    if (sum_w == 0) { sum_w = 1.0e-200; }
    fprintf(wr, "unbalance = %21.16f\n", (sum_bw[1] - sum_bw[0])/sum_w);
    double avg = wt_table_avg(n, wt);
    fprintf(wr, "mean =      %21.16f\n", avg);
    double sum_d2w = 0;
    for (int32_t k = 0; k < n; k++)
      { double w = wt[k];
        double d = k - avg;
        sum_d2w += d*d*w; 
      }
    double var = wt_table_var(n, wt, avg);
    fprintf(wr, "variance =  %21.16f\n", var);
    fprintf(wr, "deviation = %21.16f\n", sqrt(var));
    fprintf(wr, "\n");
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
  
bool_t wt_table_check_partition_of_constant
  ( int32_t n, 
    double wt[], 
    int32_t stride,
    double tol, 
    bool_t die
  )
  {
    demand(n >= 1, "invalid table length {n}");
    demand(stride >= 1, "invalid {stride}");
    double ws[n];
    wt_table_shifted_sum(n, wt, stride, ws);
    double wsExp = ws[0]; /* Expected value of overlapped windows. */
    for (int32_t k = 1; k < n; k++) 
      { if (fabs(ws[k] - wsExp) > tol)
          { if (die)
              { fprintf(stderr, "table is not  partition of constant");
                fprintf(stderr, " {ws[%d] = %18.16f", k, ws[k]);
                double err = ws[k] - wsExp;
                fprintf(stderr, " err = %24.16e\n", err);
                assert(FALSE);
              }
            else
              { return FALSE; }
          }
      }
    return TRUE;
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
