#define PROG_NAME "test_wt_table"
#define PROG_DESC "test of {wt_table.h}"
#define PROG_VERS "1.0"

/* Last edited on 2017-06-12 00:10:09 by stolfilocal */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_hermite3_COPYRIGHT \
  "Copyright © 2017  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <wt_table.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>

int main(int argn, char **argv);

void do_test_basics(void);
  /* Basic consistency tests. */

void do_test_print(int nw, char *tname);
  /* Generates tables and prints them. */

int main (int argc, char **argv)
  {
    do_test_basics();
    for (int nw = 1; nw <= 13; nw = 3*nw/2+1)
      { do_test_print(nw, "gaussian");
        do_test_print(nw, "binomial");
        do_test_print(nw, "triangular");
        do_test_print(nw, "hann");
      }
        
    return 0;
  }

void do_test_basics(void)
  {
    double sigma = 3.0;
    /* Paranoia - check functions: */
    for (int nw = 1; nw <= 11; nw++)
      { /* Compute sum of all entries in table: */
        double win = 0.0;
        for (int k = 0; k < nw; k++) { win += wt_table_gaussian_entry(nw, k, sigma); }
        /* Compute total mass outside the table: */
        double wot = wt_table_gaussian_loss(nw, sigma);
        fprintf(stderr, "gaussian sigma =  %20.18f nw = %d", sigma, nw);
        fprintf(stderr, "  inside = %20.18f  outside =  %20.18f  sum = %20.18f\n", win, wot, win+wot);
        assert(fabs(1 - (win+wot)) < 1.0e-10);
      }
  }

void do_test_print(int nw, char *tname)
  {
    fprintf(stderr, "=== %s =============================================\n", tname);
    
    double wt[nw];         
    
    if (strcmp(tname, "gaussian") == 0)
      { double sigma = nw/5.0;
        wt_table_fill_gaussian(sigma, nw, wt); 
      }
    else if (strcmp(tname, "binomial") == 0)
      { wt_table_fill_binomial(nw, wt); 
      }
    else if (strcmp(tname, "triangular") == 0)
      { wt_table_fill_triangular(nw, wt); }
    else if (strcmp(tname, "hann") == 0)
      { wt_table_fill_hann(nw, wt); }
    else 
      { assert(FALSE); }

    wt_table_print(stderr, tname, nw, wt);
    wt_table_check_normalization(nw, wt, 1.0e-8, TRUE);
    
    /* !!! To be expanded !!! */
   
  }

