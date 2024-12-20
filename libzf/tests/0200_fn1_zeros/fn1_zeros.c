/* Tests univariate zero-finders with AA and IA test. */
/* Last edited on 2024-12-05 10:41:07 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <flt.h>
#include <ia.h>
#include <aa.h>

#include <fn1_functions.h>
#include <fn1_zf.h>

/*** INTERNAL PROTOTYPES ***/

int32_t main(int32_t argc, char **argv);

/*** MAIN PROGRAM ***/

int32_t main(int32_t argc, char **argv)
  { 
    aa_init();
    fprintf(stderr, "aa_stack_top = %p (%lud)\n", aa_top(), (long unsigned)aa_top());

    /* Get function tag, and subd psteps: */
    affirm(argc == 3, "wrong number of arguments");
    char *ftag = argv[1];
    int32_t psteps = atoi(argv[2]);
    fprintf(stderr, "function = %s plot steps = %d\n", ftag, psteps);
    
    /* Get function data: */
    fn1_data_t f = fn1_from_tag(ftag);

    /* Construct the output filename prefix: */
    char *prefix = jsprintf("zf1_%s", ftag);
    
    fn1_zf_find_and_plot_zeros
      ( prefix,
        f.eval_fp,
        f.eval_ia,
        f.diff_ia,
        f.eval_aa,
        f.descr,
        f.xd, f.yd,
        (float)f.epsilon, 
        (float)f.delta,
        psteps
      );
    return (0);
  }

