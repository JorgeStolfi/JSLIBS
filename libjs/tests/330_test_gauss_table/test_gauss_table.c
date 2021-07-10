#define PROG_NAME "test_gauss_table"
#define PROG_DESC "test of {gauss_table.h}"
#define PROG_VERS "1.0"

/* Last edited on 2021-07-04 05:11:37 by jstolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_hermite3_COPYRIGHT \
  "Copyright © 2017  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>
#include <wt_table.h>

#include <gauss_table.h>

int32_t main(int32_t argn, char **argv);

void do_test_print(int32_t nw, bool_t normSum, bool_t folded);
  /* Generates table with {nw} entries and options { normSum,folded} 
    and prints it. */

void compute_folded(int32_t nw, double avg, double dev, double wc[]);
  /* Computes a folded table {wc[0..nw-1]} by building a big enough table
    with the lib function and then folding it. Also checks consistency
    with {gauss_table_folded_bell}. */
    
void compute_unfolded(int32_t nw, double avg, double dev, double wc[]);
  /* Computes a basic table {wc[0..nw-1]} by the definition, without shortcuts. */

int32_t main (int32_t argc, char **argv)
  {
    for (int32_t nw = 1; nw <= 13; nw = 3*nw/2+1)
      { do_test_print(nw, FALSE, FALSE);
        do_test_print(nw, TRUE,  FALSE);
        do_test_print(nw, FALSE, TRUE);
        do_test_print(nw, TRUE,  TRUE);
      }
        
    return 0;
  }

void do_test_print(int32_t nw, bool_t normSum, bool_t folded)
  {
    fprintf(stderr, "=== testing nw = %d =============================================\n", nw);
    
    double avg = 0.37 * nw;
    double dev = 0.21 * nw;
    double *wt = gauss_table_make(nw, avg, dev, normSum, folded);         

    wt_table_print(stderr, "gauss_table_make", nw, wt, 0);
    if (normSum) { wt_table_check_normalization(nw, wt, 1.0e-8, TRUE); }
    
    /* Build a table {wc[0..nw-1]} independently of the lib: */
    double *wc = notnull(malloc(nw*sizeof(double)), "no mem");
    if (folded) 
      { compute_folded(nw, avg, dev, wc); }
    else
      { compute_unfolded(nw, avg, dev, wc); }
    /* Normalize {wc} if so requested: */
    if (normSum) 
      { double sum = 1.0e-100;
        for (int32_t i = 0; i < nw; i++) { sum += wc[i]; }
        for (int32_t i = 0; i < nw; i++) { wc[i] /= sum; }
      }
    /* Compare {wc} with library table  {wt}: */
    bool_t ok = TRUE;
    for (int32_t i = 0; i < nw; i++) 
      { if (fabs(wt[i] - wc[i]) > 1.0e-12)
          { fprintf(stderr, "** {gauss_table_make} error:");
            fprintf(stderr, " wt[i] = %.15f  wc[i] = %.15f",  wt[i], wc[i]);
            fprintf(stderr, " diff = %.15f",  wt[i] - wc[i]);
            fprintf(stderr, "\n");
            ok = FALSE;
          }
      }
    assert(ok);
  }
    
void compute_folded(int32_t nw, double avg, double dev, double wc[])
  { /* Allocate a table {wf} big enough to hold the full distribution: */
    int32_t mw = 7*nw;
    int32_t nf = mw + nw + mw;
    double *wf = gauss_table_make(nf, avg + (double)mw, dev, FALSE, FALSE);
    assert (wf[0] < 1.0e-15);
    assert (wf[nf-1] < 1.0e-15);
    /* Accumulate the entries congruent modulo {nw} into central segment, smallest first: */
    for (int32_t i = 0; i < nw; i++) { wc[i] = wf[mw +i]; }
    for (int32_t i = 0; i < mw; i++)
      { int32_t j = nf - 1 - i;
        wc[i] += wf[i] + wf[j];
      }
    /* Consistency with {gauss_table_folded_bell}: */
    bool_t ok = TRUE;
    for (int32_t i = 0; i < nw; i++) 
      { double z = ((double)i) - avg;
        double wgi = gauss_table_folded_bell(z, dev, nw);
        if (fabs(wgi - wc[i]) > 1.0e-12)
          { fprintf(stderr, "** {gauss_table_folded_bell} error:");
            fprintf(stderr, " wgi = %.15f  wc[i] = %.15f",  wgi, wc[i]);
            fprintf(stderr, " diff = %.15f",  wgi - wc[i]);
            fprintf(stderr, "\n");
            ok = FALSE;
          }
      }
    assert(ok);
  }
    
void compute_unfolded(int32_t nw, double avg, double dev, double wc[])
  { for (int32_t i = 0; i < nw; i++)
      { double z = (((double)i) - avg)/dev;
        wc[i] = exp(-z*z/2);
      }
  }

