/* Tracing zeros of a function with various grids and validated models. */
/* Last edited on 2024-12-05 10:41:12 by stolfi */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <affirm.h>
#include <bool.h>

#include <fn2_zf_quad.h>
#include <fn2_zf_grid.h>
#include <fn2_functions.h>

/*** INTERNAL PROTOTYPES ***/

int32_t main(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
  {
    aa_init();
    fprintf(stderr, "aa_stack_top = %p (%lud)\n", aa_top(), (long unsigned)aa_top());
    
    /* Get plot method, function tag, and subd order: */
    affirm(argc == 4, "wrong number of arguments");
    char *ftag = argv[1];
    char *meth = argv[2];
    int32_t order = atoi(argv[3]);
    fprintf(stderr, "method = %s function = %s order = %d\n", meth, ftag, order);
    
    /* Get function data: */
    fn2_data_t f = fn2_from_tag(ftag);

    /* Construct the output filename prefix: */
    char *prefix = jsprintf("zf2_%s_%s_%03d", ftag, meth, order);
    
    if (strcmp(meth, "quad") == 0)
      { fn2_zf_quad_plots(
          prefix,
          f.descr,
          f.eval_fp,
          f.eval_ia,
          f.eval_aa,
          f.xd, f.yd,
          order, 
          128
        );
      }
    else
      { fn2_zf_grid_plots(
          prefix,
          f.descr,
          f.eval_fp,
          f.eval_ia,
          f.eval_aa,
          f.xd, f.yd,
          order, 
          128
        );
      }

    return (0);
  }
 
