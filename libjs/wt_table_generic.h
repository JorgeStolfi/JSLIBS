#ifndef wt_table_generic_H
#define wt_table_generic_H

/* Generic weight tables with variable kind */
/* Last edited on 2023-11-26 13:21:08 by stolfi */

#define wt_table_generic_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>

typedef enum 
  { wt_table_kind_GAUSSIAN,
    wt_table_kind_HANN,
    wt_table_kind_BINOMIAL,
    wt_table_kind_TRIANGULAR,
    wt_table_kind_UNIFORM,
    wt_table_kind_INVALID
  } wt_table_kind_t;
  /* A code identifying a discrete windowing function type. See t_table_gaussian.h},
    {wt_table_hann.h}, etc.
    
    The {wt_table_kind_INVALID} value is not a valid kind. 
    It should be the last in the {enum}. */
    
wt_table_kind_t wt_table_kind_from_string(char *name);
  /* Returns a {wt_table_kind_t} from the given name, as explained in
    {wt_table_kind_from_string_INFO} below.  If the name is not recognized, returns {wt_table_kind_INVALID}. */
    
#define wt_table_kind_from_string_INFO \
  "        \"gaussian\"    see {wt_table_gaussian_fill} ({sigma})." \
  "        \"triangular\"  see {wt_table_triangular_fill}.\n" \
  "        \"binomial\"    see {wt_table_binomial_fill}.\n" \
  "        \"uniform\"     see {wt_table_uniform_fill} ({val}).\n" \
  "        \"hann\"        see {wt_table_hann_fill} ({flat}).\n" \
  "    The name can be any non-empty prefix of the above that is distinctive."

char *wt_table_kind_to_string(wt_table_kind_t kind);
  /* Returns a constant (non-heap) string with the name corresponding to the 
    given {kind}. Returns "invalid" for {wt_table_kind_INVALID}.*/ 

void wt_table_fill(wt_table_kind_t kind, int32_t n, double parm, double wt[], bool_t norm, int32_t *stride_P);
  /* Stores into {wt[0..n-1]} a window weight table of the given {kind}
    (which must not be {wt_table_kind_INVALID}.
    If the table has a parameter (such as the deviation for Gaussian,
    the flat fraction for Hann, etc.), that parameter is set to {parm},
    otherwise {parm} is ignored. 
    
    If {norm} is true, the table is normalized so that
    the sum of all entries is 1.  
    
    The table is a /partition of constant/ (PoC) with /stride/ {m}
    if the sum of infinitely many copies of the table, shifted by
    integer multiples of {m}, is a constant function.  Note that if {d} is a
    positive divisor of {m}, then the table will be a PoC also with
    stride {d}. If the table is a PoC, the largest value of {m} is
    returned in {*stride_P}. Otherwise {*stride_P} is set to 0. If
    {*stride_P} is {NULL}, it is ignored. */

double_vec_t wt_table_make(wt_table_kind_t kind, int32_t n, double parm, bool_t norm, int32_t *stride_P);
  /* Allocates a new {double_vec_t} {wt} with {wt.ne = n} elements and fills {wt.e[0..n-1]} 
    with {wt_table_fill(kind, n, parm, wt.e,norm,stride_P)}. */

#endif
