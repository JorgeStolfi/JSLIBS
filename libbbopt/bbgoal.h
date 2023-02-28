/* A generic goal function for branch-and-bound optimization programs. */
/* Last edited on 2023-02-20 06:49:07 by stolfi */

#ifndef bbgoal_H
#define bbgoal_H

#define _GNU_SOURCE
#include <stdint.h>

#include <flt.h>
#include <ia.h>

typedef struct bbgoal_data_t
  {
    int32_t dim;      /* Dimension {d} of domain. */
    
    Float (*eval_fp)(Float *x);
      /* Returns an approximate value of the goal function
        for {xr[0..d-1]}. */

    Interval (*eval_ia)(Interval *xr);
      /* Returns a guaranteed enclosure for the exact value of the goal
        function, for any argument vector in the box {xr[0..d-1]}. Note
        that this enclosure may not contain the approximate floating-point
        value returned by {F} above. */

    void (*true_opt)(Interval *xr, Interval *sr);
      /* Returns in {sr[0..d-1]} the smallest possible box that contains
        some true global minimum of {F1} in the domain {xr[0..d-1]}. */

    Interval (*plot_range)(Interval *xr);
      /* Returns an interval of values of the goal function that is
        suitable as the function-axis plot range of graphs of {F}, when
        the arguments is restricted to {xr[0..d-1]}. May be larger or
        smaller than {F_ia(xr)}. */

    char *tag;    /* Function tag, e.g. "F1-IA". */
    char *descr;  /* Function description, possibly multi-line. */
  } bbgoal_data_t;

bbgoal_data_t bbgoal_from_tag(char *tag);
  /* Returns a canned goal function given its tag. */

#endif
