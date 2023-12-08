/* See wt_table_args_parse.h */
/* Last edited on 2023-11-25 18:59:22 by stolfi */

#define wt_table_args_parse_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <jsstring.h>
#include <argparser.h>

#include <wt_table.h>

#include <wt_table_args_parse.h>
  
wt_table_kind_t wt_table_args_parse_kind(argparser_t *pp)
  { 
    char *name = argparser_get_next_non_keyword(pp);
    wt_table_kind_t kind = wt_table_kind_from_string(name);
    if (kind == wt_table_kind_INVALID) { argparser_error(pp, "unrecognized weight table kind"); }
    return kind;
  }

double_vec_t wt_table_args_parse_weights(argparser_t *pp, bool_t unitNorm)
  { /* Allocate the weight vector, with size unknown: */
    double_vec_t w = double_vec_new(50);
    /* Parse the follwoing numeric arguments, save them in {w[0..nw-1]}: */
    int32_t nw = 0;
    while (argparser_next_is_number(pp))
      { double wt = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        double_vec_expand(&w, nw);
        w.e[nw] = wt; nw++;
      }
    /* Parse the optional "/ {DENOM}" args: */
    double denom = 1.0;
    if (argparser_keyword_present_next(pp, "/"))
      { double den = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        if (! unitNorm) { denom = den; }
      }
    if (unitNorm)
      { /* Set [denom} to the sum of all weights: */
        denom = 0;
        for (int32_t k = 0; k < nw; k++) { denom += w.e[k]; }
      }
    if (denom != 1)
      { /* Divide all weights by {denom}: */
        for (int32_t k = 0; k < nw; k++) { w.e[k] /= denom; }
      }
    double_vec_trim(&w, nw);
    return w;
  }
