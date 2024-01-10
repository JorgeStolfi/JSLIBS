/* Plots function susing AA and IA. */
/* Last edited on 2016-12-26 22:19:40 by stolfilocal */

#define _GNU_SOURCE

#include <allgraphs.h>
#include <fn1_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <affirm.h>
#include <bool.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/*** INTERNAL PROTOTYPES ***/

int main(int argc, char **argv);

/*** MAIN PROGRAM ***/

int main(int argc, char **argv)
  { 
    aa_init();
    fprintf(stderr, "aa_stack_top = %p (%lu)\n", aa_top(), (uintptr_t)aa_top());

    /* Get function tag and plot steps: */
    affirm(argc == 3, "wrong number of arguments");
    char *ftag = argv[1];
    int psteps = atoi(argv[2]);
    fprintf(stderr, "function = %s plot steps = %d\n", ftag, psteps);
    
    /* Get function data: */
    fn1_data_t f = fn1_from_tag(ftag);

    /* Construct the output filename prefix: */
    char *fileprefix = NULL;
    asprintf(&fileprefix, "pf-%s-", ftag);
    
    bool_t epsformat;
    for (epsformat = FALSE; epsformat <= TRUE; epsformat++)
      { allgraphs_plot
          ( fileprefix,
            epsformat,
            f.descr,
            f.eval_fp,
            f.eval_ia,
            f.diff_ia,
            f.eval_aa,
            f.xd.lo, f.xd.hi,
            f.yd.lo, f.yd.hi,
            f.nsub,
            psteps
          );
      }
    return 0;
  }

