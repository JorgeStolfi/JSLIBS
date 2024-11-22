#define PROG_NAME "test_gauss_table"
#define PROG_DESC "test of {gauss_table.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-19 06:14:51 by stolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_gauss_table_COPYRIGHT \
  "Copyright © 2017  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>
#include <wt_table.h>
#include <wt_table_gaussian.h>

#include <gauss_table.h>

int32_t main(int32_t argn, char **argv);

void test_gauss_table_make__gauss_table_folded_bell(uint32_t nw, bool_t norm, bool_t fold);
  /* Generates table with {nw} entries and options {norm,fold} 
    and prints it. */

void tgat_compute_folded(uint32_t nw, double avg, double dev, double wc[]);
  /* Computes a folded table {wc[0..nw-1]} by building a big enough table
    with the lib function and then folding it. Also checks consistency
    with {gauss_table_folded_bell}. */
    
void tgat_compute_unfolded(uint32_t nw, double avg, double dev, double wc[]);
  /* Computes a basic table {wc[0..nw-1]} by the definition, without shortcuts. */

int32_t main (int32_t argc, char **argv)
  {
    for (int32_t nw = 1; nw <= 13; nw = 3*nw/2+1)
      { test_gauss_table_make__gauss_table_folded_bell(nw, FALSE, FALSE);
        test_gauss_table_make__gauss_table_folded_bell(nw, TRUE,  FALSE);
        test_gauss_table_make__gauss_table_folded_bell(nw, FALSE, TRUE);
        test_gauss_table_make__gauss_table_folded_bell(nw, TRUE,  TRUE);
      }
        
    return 0;
  }

void test_gauss_table_make__gauss_table_folded_bell(uint32_t nw, bool_t norm, bool_t fold)
  {
    fprintf(stderr, "=== testing nw = %d norm = %c fold = %c ===\n", nw, "FT"[norm], "FT"[fold]);
    
    double avg = 0.37 * nw;
    double dev = 0.21 * nw;
    fprintf(stderr, "  avg = %+24.16e  dev = %24.16e\n", avg, dev);
    double *wt = gauss_table_make(nw, avg, dev, norm, fold);

    wt_table_print(stderr, "wt_table_gaussian_make", nw, wt, 0);
    if (norm) { wt_table_check_normalization(nw, wt, 1.0e-8, TRUE); }
    
    /* Build a table {wc[0..nw-1]} independently of the lib: */
    double *wc = talloc(nw, double);
    if (fold) 
      { tgat_compute_folded(nw, avg, dev, wc); }
    else
      { tgat_compute_unfolded(nw, avg, dev, wc); }
    /* Normalize {wc} if so requested: */
    if (norm) 
      { double sum = 1.0e-100;
        for (int32_t i = 0; i < nw; i++) { sum += wc[i]; }
        for (int32_t i = 0; i < nw; i++) { wc[i] /= sum; }
      }
    else if (fold)
      { /* Normalize to unity at {avg}: */
        double wmax = gauss_table_folded_bell(0.0, dev, nw);
        assert(wmax > 0.0);
        for (int32_t i = 0; i < nw; i++) { wc[i] /= wmax; }
      }
      
    /* Compare {wc} with library table  {wt}: */
    bool_t ok = TRUE;
    for (int32_t i = 0; i < nw; i++) 
      { if (fabs(wt[i] - wc[i]) > 1.0e-12)
          { fprintf(stderr, "** {gauss_table_make} error:");
            fprintf(stderr, " i = %d wt[i] = %.15f  wc[i] = %.15f", i, wt[i], wc[i]);
            fprintf(stderr, " diff = %.15f",  wt[i] - wc[i]);
            fprintf(stderr, "\n");
            ok = FALSE;
          }
      }
    assert(ok);
  }
    
void tgat_compute_folded(uint32_t nw, double avg, double dev, double wc[])
  { /* Reduce {avg} modulo {nw}: */
    avg = avg - nw*floor(avg/nw);
    while (avg >= nw) { avg -= nw; }
    while (avg < 0) { avg += nw; }
    assert(avg >= 0);
    assert(avg < nw);
    /* Allocate a table {wf} big enough to hold the full distribution: */
    uint32_t mw = 10*nw;
    uint32_t nf = mw + nw + mw;
    /* Build unfolded table with average at {mw+avg}: */
    double *wf = gauss_table_make(nf, avg + (double)mw, dev, FALSE, FALSE);
    assert (wf[0] < 1.0e-16);
    assert (wf[nf-1] < 1.0e-16);
    /* Accumulate the entries congruent modulo {nw} into central segment, smallest first: */
    assert((mw % nw) == 0); 
    for (int32_t i = 0; i < mw; i++)
      { uint32_t ki = (i % nw);
        wc[ki] += wf[i];
        uint32_t j = nf - 1 - i;
        uint32_t kj = (j % nw);
        wc[kj] += wf[j];
      }
    for (int32_t i = 0; i < nw; i++) { wc[i] += wf[mw + i]; }
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
    
void tgat_compute_unfolded(uint32_t nw, double avg, double dev, double wc[])
  { for (int32_t i = 0; i < nw; i++)
      { double z = (((double)i) - avg)/dev;
        wc[i] = exp(-z*z/2);
      }
  }

