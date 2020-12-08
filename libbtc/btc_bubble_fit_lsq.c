/* See {btc_bubble_fit_lsq.h} */
/* Last edited on 2015-04-21 22:11:22 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <lsq_robust.h>

#include <btc_bubble_t.h>

#include <btc_bubble_fit_lsq.h>

void btc_bubble_fit_lsq
  ( int nd, 
    char* dt[], 
    double ap[], 
    double wt[],
    int nb, 
    btc_bubble_t bp[], 
    double bval[], 
    int maxIters,
    char* outPrefix
  )
  {
    bool_t verbose = FALSE;
    
    if (verbose)
      { fprintf(stderr, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
        fprintf(stderr, "robust least squares fitting\n");
      }

    double* Pr = notnull(malloc(nd*sizeof(double)), "no mem"); /* Probability of data point being good. */
    double* coef = notnull(malloc(nb*sizeof(double)), "no mem"); /* Coefficients in vector form. */
    
    /* Compute the coefficients by robust least squares: */
    bool_t verbose_lsq = FALSE;
    lsq_robust_fit
      ( nd, nb,
        bval, ap, wt,
        maxIters,
        coef, Pr,
        NULL,
        verbose_lsq
      );
        
    /* Save the coefficients: */
    int jb;
    for (jb = 0; jb < nb; jb++)
      { bp[jb].coef = coef[jb];
        if (verbose) { fprintf(stderr, "  coef[%02d] = %25.16e\n", jb, bp[jb].coef); }
      }
    if (verbose)
      { fprintf(stderr, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
        fprintf(stderr, "\n");
      }
    
    free(coef);
    free(Pr);
  }
  
