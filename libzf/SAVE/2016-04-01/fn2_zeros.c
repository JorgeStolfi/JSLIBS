/* Tracing zeros of a function with various grids and validated models. */
/* Last edited on 2012-07-21 12:29:39 by stolfi */

#define _GNU_SOURCE

#include <fn2_zf_quad.h>
#include <fn2_zf_grid.h>
#include <fn2_functions.h>

#include <aa.h>
#include <ia.h>
#include <flt.h>
#include <affirm.h>
#include <bool.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*** INTERNAL PROTOTYPES ***/

int main(int argc, char **argv);

int main(int argc, char **argv)
  {
    aa_init();
    fprintf(stderr, "aa_stack_top = %p (%lud)\n", aa_top(), (long unsigned)aa_top());
    
    /* Get plot method, function tag, and subd order: */
    affirm(argc == 4, "wrong number of arguments");
    char *ftag = argv[1];
    char *meth = argv[2];
    int order = atoi(argv[3]);
    fprintf(stderr, "method = %s function = %s order = %d\n", meth, ftag, order);
    
    /* Get function data: */
    fn2_data_t f = fn2_from_tag(ftag);

    /* Construct the output filename prefix: */
    char *fileprefix = NULL;
    asprintf(&fileprefix, "z2-%s-%s-%03d-", ftag, meth, order);
    
    bool_t epsf;
    for (epsf = FALSE; epsf <= TRUE; epsf++)
      { 
        if (strcmp(meth, "quad") == 0)
          { fn2_zf_quad_plots(
              fileprefix,
              epsf,
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
              fileprefix,
              epsf,
              f.descr,
              f.eval_fp,
              f.eval_ia,
              f.eval_aa,
              f.xd, f.yd,
              order, 
              128
            );
          }
      }

    return (0);
  }
 
