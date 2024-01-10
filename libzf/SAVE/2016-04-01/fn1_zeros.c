/* Tests univariate zero-finders with AA and IA test. */
/* Last edited on 2012-07-21 12:29:27 by stolfi */

#define _GNU_SOURCE

#include <fn1_zf.h>
#include <fn1_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <affirm.h>
#include <bool.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*** INTERNAL PROTOTYPES ***/

int main(int argc, char **argv);

/*** MAIN PROGRAM ***/

int main(int argc, char **argv)
  { 
    aa_init();
    fprintf(stderr, "aa_stack_top = %p (%lud)\n", aa_top(), (long unsigned)aa_top());

    /* Get function tag, and subd psteps: */
    affirm(argc == 3, "wrong number of arguments");
    char *ftag = argv[1];
    int psteps = atoi(argv[2]);
    fprintf(stderr, "function = %s plot steps = %d\n", ftag, psteps);
    
    /* Get function data: */
    fn1_data_t f = fn1_from_tag(ftag);

    /* Construct the output filename prefix: */
    char *fileprefix = NULL;
    asprintf(&fileprefix, "z1-%s-", ftag);
    
    fn1_zf_find_and_plot_zeros
      ( fileprefix,
        FALSE, 
        f.eval_fp,
        f.eval_ia,
        f.diff_ia,
        f.eval_aa,
        f.descr,
        f.xd, f.yd,
        f.epsilon, f.delta,
        psteps
      );
    return (0);
  }

